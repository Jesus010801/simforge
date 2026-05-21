# workflows/workflow_graph.py
"""
WorkflowGraph — representación DAG del plan de simulación.

Responsabilidades:
    - Construir el grafo de dependencias desde un SimulationPlan
    - Validar que sea un DAG (sin ciclos, sin dependencias rotas)
    - Ordenar topológicamente (Kahn's algorithm)
    - Detectar el camino crítico (ruta más larga sin paralelismo)
    - Identificar grupos de steps que pueden ejecutarse en paralelo
    - Generar vistas para distintos consumidores:
        to_execution_view()  → list[SimulationStep] ordenado topológicamente
        to_user_view()       → list[str] simplificado para el usuario final
        to_stage_view()      → dict[StepStage, list[SimulationStep]] agrupado
        to_parallel_waves()  → list[list[SimulationStep]] ondas de ejecución paralela
        render_mermaid()     → str en formato Mermaid con agrupación por stage

Principios:
    - NO modifica el SimulationPlan
    - NO accede a archivos ni al SystemState
    - Solo consume SimulationStep.step_id y SimulationStep.depends_on
    - Toda la información semántica (stage, title, blocking) se preserva
      pero no dirige la lógica del grafo

Invariantes garantizados después de validate():
    - Todos los step_id son únicos
    - Todas las dependencias apuntan a step_id existentes
    - El grafo es acíclico (DAG)
    - sorted_steps contiene exactamente todos los nodos
"""

from __future__ import annotations

from collections import defaultdict, deque
from typing import Iterator

from core.execution_models import (
    SimulationPlan,
    SimulationStep,
    StepStage,
    StepType,
    PlanStatus,
)


# ═══════════════════════════════════════════════════════════════════════════════
# WorkflowGraph
# ═══════════════════════════════════════════════════════════════════════════════

class WorkflowGraph:
    """
    DAG del workflow de simulación.

    Construcción:
        graph = WorkflowGraph(plan)
        graph.validate()               # lanza ValueError si hay problemas
        ordered = graph.to_execution_view()

    El grafo es inmutable después de __init__.
    validate() puede llamarse múltiples veces sin efectos secundarios.
    """

    def __init__(self, plan: SimulationPlan) -> None:
        self.plan = plan

        # ── Nodos ─────────────────────────────────────────────────────────────
        # step_id → SimulationStep
        self.nodes: dict[str, SimulationStep] = {}

        # ── Aristas dirigidas ─────────────────────────────────────────────────
        # dep_id → [step_id, ...]   (dep_id debe ejecutarse ANTES de step_id)
        self.edges: dict[str, list[str]] = defaultdict(list)

        # step_id → [dep_id, ...]   (predecesores directos de step_id)
        self.reverse_edges: dict[str, list[str]] = defaultdict(list)

        # ── Cache de sort topológico ──────────────────────────────────────────
        self._sorted_steps: list[SimulationStep] = []
        self._validated: bool = False

        self._build_graph()

    # ═══════════════════════════════════════════════════════════════════════════
    # Construcción interna
    # ═══════════════════════════════════════════════════════════════════════════

    def _build_graph(self) -> None:
        """
        Registra nodos y aristas desde plan.steps.
        Detecta step_id duplicados y dependencias rotas en este paso.
        """

        # ── Registrar nodos ────────────────────────────────────────────────────
        for step in self.plan.steps:
            if step.step_id in self.nodes:
                raise ValueError(
                    f"[WorkflowGraph] step_id duplicado: '{step.step_id}'. "
                    "Cada step debe tener un identificador único."
                )
            self.nodes[step.step_id] = step

        # ── Registrar dependencias ─────────────────────────────────────────────
        for step in self.plan.steps:
            for dep in step.depends_on:
                if dep not in self.nodes:
                    raise ValueError(
                        f"[WorkflowGraph] Dependencia rota: "
                        f"'{step.step_id}' depende de '{dep}', "
                        f"pero '{dep}' no existe en el plan."
                    )
                self.edges[dep].append(step.step_id)
                self.reverse_edges[step.step_id].append(dep)

    # ═══════════════════════════════════════════════════════════════════════════
    # Validación
    # ═══════════════════════════════════════════════════════════════════════════

    def validate(self) -> None:
        """
        Verifica que el grafo sea un DAG válido.

        Detecta:
            - Ciclos (Kahn's algorithm: si sobran nodos al final → ciclo)
            - Grafo vacío (warning, no error)

        Lanza ValueError con mensaje descriptivo si falla.
        Idempotente: puede llamarse múltiples veces.
        """
        if not self.nodes:
            # Plan vacío — válido pero inútil
            self._validated = True
            return

        self._topological_sort()

        n_sorted = len(self._sorted_steps)
        n_nodes  = len(self.nodes)

        if n_sorted != n_nodes:
            # Kahn's: si no se procesaron todos los nodos → ciclo
            sorted_ids = {s.step_id for s in self._sorted_steps}
            remaining  = [sid for sid in self.nodes if sid not in sorted_ids]
            raise ValueError(
                f"[WorkflowGraph] El workflow contiene ciclos.\n"
                f"  Nodos procesados: {n_sorted}/{n_nodes}\n"
                f"  Steps involucrados en el ciclo: {remaining}"
            )

        self._validated = True

    # ═══════════════════════════════════════════════════════════════════════════
    # Ordenamiento topológico — Kahn's Algorithm
    # ═══════════════════════════════════════════════════════════════════════════

    def _topological_sort(self) -> list[SimulationStep]:
        """
        Kahn's algorithm con desempate por stage para producir un orden
        estable y semánticamente coherente.

        Desempate: dentro de los nodos con in_degree == 0, procesar primero
        los de stages anteriores en el pipeline (preparation < parametrization
        < validation < assembly < minimization < equilibration < ...).
        Esto garantiza que el orden refleje el flujo científico esperado.
        """
        _STAGE_ORDER = {
            StepStage.PREPARATION:       0,
            StepStage.PARAMETRIZATION:   1,
            StepStage.VALIDATION:        2,
            StepStage.ASSEMBLY:          3,
            StepStage.MINIMIZATION:      4,
            StepStage.EQUILIBRATION:     5,
            StepStage.ENHANCED_SAMPLING: 6,
            StepStage.PRODUCTION:        7,
            StepStage.ANALYSIS:          8,
        }

        in_degree: dict[str, int] = {
            node: len(self.reverse_edges[node])
            for node in self.nodes
        }

        # Cola priorizada: (stage_order, step_id) para desempate determinista
        # Usamos lista + sort manual (evita importar heapq para tan pocos nodos)
        ready: list[str] = [
            sid for sid, deg in in_degree.items() if deg == 0
        ]

        def _priority(sid: str) -> tuple[int, str]:
            stage = self.nodes[sid].stage
            return (_STAGE_ORDER.get(stage, 99), sid)

        ordered: list[SimulationStep] = []

        while ready:
            ready.sort(key=_priority)
            current = ready.pop(0)
            ordered.append(self.nodes[current])

            for neighbor in self.edges[current]:
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    ready.append(neighbor)

        self._sorted_steps = ordered
        return ordered

    # ═══════════════════════════════════════════════════════════════════════════
    # Camino crítico
    # ═══════════════════════════════════════════════════════════════════════════

    def critical_path(self) -> list[SimulationStep]:
        """
        Detecta el camino crítico del DAG: la secuencia más larga de
        dependencias que determina el tiempo mínimo de ejecución si se
        paraleliza todo lo que se pueda.

        Algoritmo: DP sobre el orden topológico.
            dist[v] = max(dist[u] + 1 for u in predecesores de v)
        El camino crítico es la cadena de nodos con mayor dist.

        Sin estimaciones de tiempo real (el executor no existe todavía),
        cada step cuenta como peso 1. En el futuro, SimulationStep puede
        tener un campo `estimated_duration_h` que se usa aquí.

        Retorna lista de SimulationStep en orden desde la raíz.
        """
        if not self._sorted_steps:
            self._topological_sort()

        # dist[step_id] = longitud del camino más largo que termina en ese nodo
        dist: dict[str, int]     = {sid: 1 for sid in self.nodes}
        prev: dict[str, str | None] = {sid: None for sid in self.nodes}

        for step in self._sorted_steps:
            for dep_id in self.reverse_edges[step.step_id]:
                candidate = dist[dep_id] + 1
                if candidate > dist[step.step_id]:
                    dist[step.step_id]  = candidate
                    prev[step.step_id] = dep_id

        # Nodo final del camino crítico = el de mayor dist
        end_node = max(dist, key=lambda k: dist[k])

        # Reconstruir camino hacia atrás
        path: list[str] = []
        node: str | None = end_node
        while node is not None:
            path.append(node)
            node = prev[node]

        path.reverse()
        return [self.nodes[sid] for sid in path]

    # ═══════════════════════════════════════════════════════════════════════════
    # Ondas de ejecución paralela
    # ═══════════════════════════════════════════════════════════════════════════

    def to_parallel_waves(self) -> list[list[SimulationStep]]:
        """
        Agrupa los steps en "ondas" (waves) de ejecución paralela.

        Una onda contiene todos los steps cuyos predecesores ya están
        en ondas anteriores. Los steps dentro de una misma onda pueden
        ejecutarse simultáneamente.

        Útil para el executor cuando tiene múltiples workers disponibles
        (ej: cluster HPC con varios nodos) o para mostrar al usuario
        qué puede correr en paralelo.

        Algoritmo: BFS por niveles desde los nodos raíz.
        """
        if not self._sorted_steps:
            self._topological_sort()

        in_degree: dict[str, int] = {
            sid: len(self.reverse_edges[sid])
            for sid in self.nodes
        }

        waves:   list[list[SimulationStep]] = []
        visited: set[str]                   = set()

        current_wave_ids: list[str] = [
            sid for sid, deg in in_degree.items() if deg == 0
        ]

        while current_wave_ids:
            # Ordenar dentro de la onda por stage para coherencia visual
            current_wave_ids.sort(
                key=lambda sid: (self.nodes[sid].stage.value, sid)
            )
            wave = [self.nodes[sid] for sid in current_wave_ids]
            waves.append(wave)
            visited.update(current_wave_ids)

            next_wave_ids: list[str] = []
            for sid in current_wave_ids:
                for neighbor in self.edges[sid]:
                    if neighbor in visited:
                        continue
                    # ¿Todos los predecesores del neighbor ya fueron visitados?
                    preds = self.reverse_edges[neighbor]
                    if all(p in visited for p in preds):
                        if neighbor not in next_wave_ids:
                            next_wave_ids.append(neighbor)

            current_wave_ids = next_wave_ids

        return waves

    # ═══════════════════════════════════════════════════════════════════════════
    # Vistas públicas
    # ═══════════════════════════════════════════════════════════════════════════

    def to_execution_view(self) -> list[SimulationStep]:
        """
        Lista completa de steps en orden topológico.
        Es el contrato principal que consume WorkspaceBuilder.
        """
        if not self._sorted_steps:
            self._topological_sort()
        return list(self._sorted_steps)

    def to_stage_view(self) -> dict[StepStage, list[SimulationStep]]:
        """
        Steps agrupados por stage, en orden topológico dentro de cada grupo.
        Útil para reportes y debugging.
        """
        if not self._sorted_steps:
            self._topological_sort()

        result: dict[StepStage, list[SimulationStep]] = defaultdict(list)
        for step in self._sorted_steps:
            result[step.stage].append(step)
        return dict(result)

    def to_user_view(self) -> list[str]:
        """
        Vista simplificada para el usuario final.

        Filtra steps internos (prepare_*) y reemplaza steps manuales
        con mensajes legibles. El resultado es una lista numerada de
        acciones que el usuario entiende sin conocer los internos.
        """
        if not self._sorted_steps:
            self._topological_sort()

        _HIDDEN_PREFIXES = (
            "prepare_",
            "validate_pose_",
        )

        _MANUAL_MESSAGES = {
            "review_parametrization": "⚠  Revisión manual de parametrización requerida",
        }

        user_steps: list[str] = []

        for step in self._sorted_steps:

            # ── Ocultar steps internos de preparación ─────────────────────────
            if any(step.step_id.startswith(p) for p in _HIDDEN_PREFIXES):
                continue

            # ── Steps manuales con mensaje descriptivo ────────────────────────
            matched = False
            for prefix, msg in _MANUAL_MESSAGES.items():
                if step.step_id.startswith(prefix):
                    user_steps.append(msg)
                    matched = True
                    break
            if matched:
                continue

            # ── Steps externos ────────────────────────────────────────────────
            if step.step_type == StepType.EXTERNAL:
                user_steps.append(f"[externo] {step.title}")
                continue

            user_steps.append(step.title)

        return user_steps

    # ═══════════════════════════════════════════════════════════════════════════
    # Helpers de consulta
    # ═══════════════════════════════════════════════════════════════════════════

    def predecessors(self, step_id: str) -> list[SimulationStep]:
        """Todos los predecesores directos de un step."""
        return [self.nodes[sid] for sid in self.reverse_edges.get(step_id, [])]

    def successors(self, step_id: str) -> list[SimulationStep]:
        """Todos los sucesores directos de un step."""
        return [self.nodes[sid] for sid in self.edges.get(step_id, [])]

    def ancestors(self, step_id: str) -> list[SimulationStep]:
        """Todos los ancestros transitivos de un step (BFS hacia atrás)."""
        visited: set[str] = set()
        queue   = list(self.reverse_edges.get(step_id, []))
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue
            visited.add(current)
            queue.extend(self.reverse_edges.get(current, []))
        return [self.nodes[sid] for sid in visited]

    def descendants(self, step_id: str) -> list[SimulationStep]:
        """Todos los descendientes transitivos de un step (BFS hacia adelante)."""
        visited: set[str] = set()
        queue   = list(self.edges.get(step_id, []))
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue
            visited.add(current)
            queue.extend(self.edges.get(current, []))
        return [self.nodes[sid] for sid in visited]

    def blocking_steps(self) -> list[SimulationStep]:
        """Steps marcados como blocking=True."""
        return [s for s in self.nodes.values() if s.blocking]

    def root_steps(self) -> list[SimulationStep]:
        """Steps sin dependencias (entrada del DAG)."""
        return [
            self.nodes[sid]
            for sid in self.nodes
            if not self.reverse_edges[sid]
        ]

    def leaf_steps(self) -> list[SimulationStep]:
        """Steps sin sucesores (salida del DAG)."""
        return [
            self.nodes[sid]
            for sid in self.nodes
            if not self.edges[sid]
        ]

    # ═══════════════════════════════════════════════════════════════════════════
    # Iteradores
    # ═══════════════════════════════════════════════════════════════════════════

    def __iter__(self) -> Iterator[SimulationStep]:
        """Itera en orden topológico."""
        return iter(self.to_execution_view())

    def __len__(self) -> int:
        return len(self.nodes)

    def __contains__(self, step_id: str) -> bool:
        return step_id in self.nodes

    # ═══════════════════════════════════════════════════════════════════════════
    # Estadísticas
    # ═══════════════════════════════════════════════════════════════════════════

    def stats(self) -> dict:
        """
        Resumen estadístico del grafo. Útil para logging y debugging.
        """
        if not self._sorted_steps:
            self._topological_sort()

        stage_counts: dict[str, int] = defaultdict(int)
        for step in self.nodes.values():
            stage_counts[step.stage.value] += 1

        type_counts: dict[str, int] = defaultdict(int)
        for step in self.nodes.values():
            type_counts[step.step_type.value] += 1

        cp = self.critical_path()

        return {
            "n_steps":            len(self.nodes),
            "n_edges":            sum(len(v) for v in self.edges.values()),
            "n_root_steps":       len(self.root_steps()),
            "n_leaf_steps":       len(self.leaf_steps()),
            "n_blocking_steps":   len(self.blocking_steps()),
            "n_parallel_waves":   len(self.to_parallel_waves()),
            "critical_path_len":  len(cp),
            "critical_path_ids":  [s.step_id for s in cp],
            "stage_counts":       dict(stage_counts),
            "type_counts":        dict(type_counts),
            "plan_status":        self.plan.status.value,
            "n_blocking_issues":  len(self.plan.blocking_issues),
        }

    # ═══════════════════════════════════════════════════════════════════════════
    # Mermaid export
    # ═══════════════════════════════════════════════════════════════════════════

    def render_mermaid(self, group_by_stage: bool = True) -> str:
        """
        Exporta el workflow como Mermaid flowchart.

        group_by_stage=True:  agrupa nodos en subgraphs por stage.
                              Más útil para visualización.
        group_by_stage=False: grafo plano, más limpio para DAGs pequeños.

        Formato del nodo:
            step_id["title"]           → step automático
            step_id{{"title"}}         → step manual (rombo)
            step_id(["title"])         → step externo (estadio)

        Los steps blocking se marcan en rojo vía classDef.
        """
        if not self._sorted_steps:
            self._topological_sort()

        lines: list[str] = ["graph TD"]

        # ── Estilos ───────────────────────────────────────────────────────────
        lines.append("")
        lines.append("    classDef blocking fill:#ff6b6b,stroke:#c0392b,color:#fff")
        lines.append("    classDef manual   fill:#f39c12,stroke:#e67e22,color:#fff")
        lines.append("    classDef external fill:#8e44ad,stroke:#6c3483,color:#fff")
        lines.append("    classDef default  fill:#2ecc71,stroke:#27ae60,color:#fff")
        lines.append("")

        def _node_shape(step: SimulationStep) -> str:
            title = step.title.replace('"', "'")
            if step.step_type == StepType.MANUAL:
                return f'{step.step_id}{{{{"{title}"}}}}'
            elif step.step_type == StepType.EXTERNAL:
                return f'{step.step_id}(["{title}"])'
            else:
                return f'{step.step_id}["{title}"]'

        def _node_class(step: SimulationStep) -> str:
            if step.blocking:
                return "blocking"
            if step.step_type == StepType.MANUAL:
                return "manual"
            if step.step_type == StepType.EXTERNAL:
                return "external"
            return "default"

        if group_by_stage:
            # ── Subgraphs por stage ───────────────────────────────────────────
            stage_view = self.to_stage_view()

            _STAGE_LABELS = {
                StepStage.PREPARATION:       "🔧 Preparation",
                StepStage.PARAMETRIZATION:   "⚗️  Parametrization",
                StepStage.VALIDATION:        "✅ Validation",
                StepStage.ASSEMBLY:          "🏗️  Assembly",
                StepStage.MINIMIZATION:      "⚡ Minimization",
                StepStage.EQUILIBRATION:     "🌡️  Equilibration",
                StepStage.ENHANCED_SAMPLING: "🔬 Enhanced Sampling",
                StepStage.PRODUCTION:        "🚀 Production MD",
                StepStage.ANALYSIS:          "📊 Analysis",
            }

            for stage, steps in stage_view.items():
                label = _STAGE_LABELS.get(stage, stage.value)
                # Sanitizar: Mermaid no acepta algunos caracteres en subgraph id
                subgraph_id = f"sg_{stage.value}"
                lines.append(f'    subgraph {subgraph_id}["{label}"]')
                for step in steps:
                    lines.append(f"        {_node_shape(step)}")
                lines.append("    end")
                lines.append("")
        else:
            # ── Grafo plano ───────────────────────────────────────────────────
            for step in self._sorted_steps:
                lines.append(f"    {_node_shape(step)}")
            lines.append("")

        # ── Aristas ───────────────────────────────────────────────────────────
        for step in self._sorted_steps:
            for dep in step.depends_on:
                lines.append(f"    {dep} --> {step.step_id}")

        lines.append("")

        # ── Clases por nodo ───────────────────────────────────────────────────
        class_groups: dict[str, list[str]] = defaultdict(list)
        for step in self._sorted_steps:
            cls = _node_class(step)
            class_groups[cls].append(step.step_id)

        for cls, ids in class_groups.items():
            lines.append(f"    class {','.join(ids)} {cls}")

        return "\n".join(lines)

    # ═══════════════════════════════════════════════════════════════════════════
    # __repr__
    # ═══════════════════════════════════════════════════════════════════════════

    def __repr__(self) -> str:
        validated = "✓" if self._validated else "?"
        return (
            f"WorkflowGraph("
            f"steps={len(self.nodes)}, "
            f"edges={sum(len(v) for v in self.edges.values())}, "
            f"validated={validated})"
        )