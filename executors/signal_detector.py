# executors/signal_detector.py
"""
AdaptiveReasoner — motor de diagnóstico y remediación de SimForge.

Responsabilidades:
    diagnose()          → lee stdout/stderr/exit_code → DiagnosisResult
    plan_remediation()  → DiagnosisResult → RemediationPlan con acciones concretas

Principios arquitectónicos:
    - NO ejecuta comandos
    - NO modifica archivos
    - NO conoce el SystemState (solo conoce el output del executor)
    - Solo lee texto y produce planes estructurados

El razonamiento está organizado en dos capas:

    Capa 1 — Detección de señales (pattern matching sobre stdout/stderr)
        Cada señal es un patrón de texto que identifica una categoría de error.
        Las señales tienen prioridad: las más específicas se evalúan primero.
        La primera señal que hace match determina la categoría del diagnóstico.

    Capa 2 — Generación de plan (tabla de estrategias por categoría)
        Cada categoría de error tiene una función de planificación asociada.
        Las funciones de planificación generan acciones concretas:
        patch_mdp, scale_timestep, inject_restraints, reset_step, etc.

Extensibilidad:
    - Agregar nueva señal → añadir entrada a _SIGNAL_PATTERNS
    - Agregar nueva estrategia → añadir función _plan_* y registrarla en _PLANNERS
    - El resto del sistema no cambia

Heurísticas implementadas:

    LINCS_WARNING / LINCS_FATAL
        → Reducir dt (0.002 → 0.001), agregar position_restraints,
          reducir nstlist si está alto

    NAN_ENERGY / EXPLODING_SYSTEM
        → Reducir dt más agresivamente (0.001), agregar emtol más bajo
          en minimización previa, verificar carga del sistema

    FMAX_NOT_CONVERGED
        → Aumentar nsteps de minimización (×2), bajar emtol,
          cambiar integrador steep → l-bfgs si ya intentó steep

    POOR_EQUILIBRATION
        → Reducir tau_t, agregar restraints de posición,
          extender nsteps de equilibración

    MISSING_PARAMETER / ATOM_TYPE_MISMATCH
        → Fatal: requiere intervención humana (reparametrización)

    TIMEOUT
        → Reducir nsteps a 50% para estimar tiempo, notificar

    NONZERO_EXIT / MISSING_OUTPUT / UNKNOWN
        → Reset + retry limpio (hasta max_retries)
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Callable

from executors.execution_state import StepExecutionRecord
from executors.remediation_models import (
    DiagnosisResult,
    RemediationPlan,
    RemediationAction,
    ErrorCategory,
    ErrorSeverity,
    ActionType,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Capa 1 — Señales de diagnóstico
# ═══════════════════════════════════════════════════════════════════════════════

class _Signal:
    """Patrón de detección de error con metadatos."""

    def __init__(
        self,
        category:   ErrorCategory,
        severity:   ErrorSeverity,
        patterns:   list[str],          # regex sobre stdout+stderr
        confidence: float = 0.90,
        explanation: str = "",
    ):
        self.category    = category
        self.severity    = severity
        self.patterns    = [re.compile(p, re.IGNORECASE | re.MULTILINE) for p in patterns]
        self.confidence  = confidence
        self.explanation = explanation

    def match(self, text: str) -> tuple[bool, str]:
        """
        Retorna (matched, primary_signal).
        primary_signal es la primera línea que hizo match.
        """
        for pat in self.patterns:
            m = pat.search(text)
            if m:
                # Extraer la línea completa que contiene el match
                start = text.rfind("\n", 0, m.start()) + 1
                end   = text.find("\n", m.end())
                line  = text[start:end if end != -1 else None].strip()
                return True, line
        return False, ""


# Orden de prioridad: señales más específicas primero.
# La primera que hace match gana.
_SIGNAL_PATTERNS: list[_Signal] = [

    # ── LINCS fatal (antes que warning para no confundir) ─────────────────────
    _Signal(
        category    = ErrorCategory.LINCS_FATAL,
        severity    = ErrorSeverity.RECOVERABLE,
        patterns    = [
            r"LINCS warning.*too many times",
            r"step\s+\d+,\s+LINCS\s+LINCS WARNING",
            r"Segmentation fault.*lincs",
        ],
        confidence  = 0.92,
        explanation = "LINCS no pudo corregir la geometría de los enlaces. "
                      "El sistema está inestable — timestep demasiado grande "
                      "o estructura inicial con clashes severos.",
    ),

    # ── LINCS warning (recuperable) ────────────────────────────────────────────
    _Signal(
        category    = ErrorCategory.LINCS_WARNING,
        severity    = ErrorSeverity.RECOVERABLE,
        patterns    = [
            r"LINCS warning",
            r"angle.*LINCS",
            r"Relative constraint deviation",
        ],
        confidence  = 0.88,
        explanation = "Advertencias de LINCS detectadas. El timestep puede ser "
                      "demasiado grande para este sistema.",
    ),

    # ── NaN/Inf en energías ────────────────────────────────────────────────────
    _Signal(
        category    = ErrorCategory.NAN_ENERGY,
        severity    = ErrorSeverity.RECOVERABLE,
        patterns    = [
            r"Potential Energy\s*=\s*[-]?nan",
            r"Potential Energy\s*=\s*[-]?inf",
            r"NaN.*Potential",
            r"Epot\s*=\s*nan",
            r"non-finite.*energy",
        ],
        confidence  = 0.95,
        explanation = "Energía potencial no finita (NaN/Inf). El sistema explotó — "
                      "átomos solapados o parámetros incompatibles.",
    ),

    # ── Sistema que explota ────────────────────────────────────────────────────
    _Signal(
        category    = ErrorCategory.EXPLODING_SYSTEM,
        severity    = ErrorSeverity.RECOVERABLE,
        patterns    = [
            r"atoms are not in the expected.*molecule",
            r"One or more water molecules can not be settled",
            r"coordinates.*out of range",
            r"Atom.*moved too far",
            r"step.*Atoms moved out of the box",
        ],
        confidence  = 0.90,
        explanation = "Átomos con coordenadas fuera de rango — sistema explosivo. "
                      "La minimización previa fue insuficiente o el timestep "
                      "es demasiado grande.",
    ),

    # ── Fmax no convergió ─────────────────────────────────────────────────────
    _Signal(
        category    = ErrorCategory.FMAX_NOT_CONVERGED,
        severity    = ErrorSeverity.RECOVERABLE,
        patterns    = [
            r"Norm of force\s*=\s*[\d.]+\s*\(exceeds",
            r"Fmax\s*=\s*[\d.eE+\-]+\s*>",
            r"Energy minimization has stopped.*machine precision",
            r"did not converge to Fmax",
        ],
        confidence  = 0.93,
        explanation = "La minimización energética no convergió al Fmax objetivo. "
                      "Se necesitan más pasos o un algoritmo más robusto.",
    ),

    # ── Mala equilibración ────────────────────────────────────────────────────
    _Signal(
        category    = ErrorCategory.POOR_EQUILIBRATION,
        severity    = ErrorSeverity.RECOVERABLE,
        patterns    = [
            r"Temperature.*diverged",
            r"Pressure.*oscillat",
            r"tau_t.*too small",
            r"velocity.*too large",
            r"System temperature.*\d{4,}",   # temperatura > 999K
        ],
        confidence  = 0.82,
        explanation = "Temperatura o presión inestable durante equilibración. "
                      "El sistema no está termostatizado correctamente.",
    ),

    # ── Parámetro de FF faltante ───────────────────────────────────────────────
    _Signal(
        category    = ErrorCategory.MISSING_PARAMETER,
        severity    = ErrorSeverity.FATAL,
        patterns    = [
            r"No default.*dihedral",
            r"No default.*angle",
            r"No default.*bond",
            r"atom type.*not found",
            r"missing parameter",
            r"No parameters for",
        ],
        confidence  = 0.95,
        explanation = "Parámetro de forcefield no encontrado. "
                      "La parametrización del ligando está incompleta.",
    ),

    # ── Tipo de átomo no reconocido ────────────────────────────────────────────
    _Signal(
        category    = ErrorCategory.ATOM_TYPE_MISMATCH,
        severity    = ErrorSeverity.FATAL,
        patterns    = [
            r"atomtype.*not found",
            r"Unknown atom type",
            r"can not find atom type",
            r"atom.*has no.*parameters",
        ],
        confidence  = 0.93,
        explanation = "Tipo de átomo no reconocido por el forcefield. "
                      "Requiere reparametrización manual.",
    ),

    # ── Desbalance de carga ────────────────────────────────────────────────────
    _Signal(
        category    = ErrorCategory.CHARGE_IMBALANCE,
        severity    = ErrorSeverity.NEEDS_REVIEW,
        patterns    = [
            r"total charge.*not an integer",
            r"System has non-zero total charge",
            r"charge.*imbalance",
        ],
        confidence  = 0.88,
        explanation = "La carga total del sistema no es entera. "
                      "Verificar neutralización con iones.",
    ),

    # ── Archivo de entrada faltante ───────────────────────────────────────────
    _Signal(
        category    = ErrorCategory.MISSING_INPUT_FILE,
        severity    = ErrorSeverity.NEEDS_REVIEW,
        patterns    = [
            r"No such file or directory",
            r"cannot open.*for reading",
            r"File.*not found",
            r"Error opening file",
        ],
        confidence  = 0.90,
        explanation = "Archivo de entrada no encontrado. "
                      "Un step previo no generó el output esperado.",
    ),

    # ── Checkpoint corrupto ────────────────────────────────────────────────────
    _Signal(
        category    = ErrorCategory.CORRUPT_CHECKPOINT,
        severity    = ErrorSeverity.RECOVERABLE,
        patterns    = [
            r"corrupt.*checkpoint",
            r"checkpoint.*invalid",
            r"Error reading checkpoint",
            r"cpt.*corrupted",
        ],
        confidence  = 0.85,
        explanation = "Checkpoint corrupto — probablemente por interrupción previa.",
    ),

    # ── Proceso bloqueado esperando input interactivo (TTY) ───────────────────
    _Signal(
        category    = ErrorCategory.INTERACTIVE_BLOCK,
        severity    = ErrorSeverity.FATAL,
        patterns    = [
            r"Select group for determining the orientation",
            r"Select a group:",
            r"Enter selection:",
            r"Choice\s*[\?:]",
            r"^Group\s+\d+\s+\(",
        ],
        confidence  = 0.98,
        explanation = "GROMACS intentó leer del terminal (TTY) para una selección interactiva. "
                      "El proceso fue terminado porque stdin no tiene terminal. "
                      "Causa más común: -princ en gmx editconf sin pipe de grupo.",
    ),

    # ── Timeout ───────────────────────────────────────────────────────────────
    _Signal(
        category    = ErrorCategory.TIMEOUT,
        severity    = ErrorSeverity.NEEDS_REVIEW,
        patterns    = [
            r"Timeout",
            r"time limit exceeded",
            r"wall.*time.*exceeded",
        ],
        confidence  = 0.98,
        explanation = "El step excedió el límite de tiempo configurado.",
    ),

    # ── Position restraint atom index out of bounds ──────────────────────────
    _Signal(
        category    = ErrorCategory.TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR,
        severity    = ErrorSeverity.FATAL,
        patterns    = [
            r"TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR",
            r"Atom index.*in position_restraints out of bounds",
            r"position_restraints.*out of bounds",
        ],
        confidence  = 0.99,
        explanation = (
            "A [ position_restraints ] file contains an atom index that exceeds the local atom count\n"
            "of its moleculetype. This happens when global system atom indices are written instead of\n"
            "local moleculetype indices (1-based, starting from 1 for each chain).\n"
            "Fix: split strong_posre.itp into per-chain files (strong_posre_Protein_chain_A.itp, …)\n"
            "and write only local indices parsed from each chain's [ atoms ] section.\n"
            "Retrying will not fix this — the restraint generator must be corrected."
        ),
    ),

    # ── Undefined moleculetype in topology [ molecules ] ─────────────────────
    _Signal(
        category    = ErrorCategory.TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR,
        severity    = ErrorSeverity.FATAL,
        patterns    = [
            r"TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR",
            r"No such moleculetype",
            r"\[ molecules \] contains .* but no included \[ moleculetype \]",
        ],
        confidence  = 0.99,
        explanation = (
            "GROMACS could not find a [ moleculetype ] definition for an entry in [ molecules ].\n"
            "Common causes:\n"
            "  1. Multichain protein: pdb2gmx generates 'Protein_chain_A 1', 'Protein_chain_B 1', etc.,\n"
            "     but SimForge wrote 'Protein 1' (stale hardcoded entry).\n"
            "  2. Lipid topology: pdb2gmx generated 'DPP 1' (one combined DPP moleculetype) but\n"
            "     SimForge wrote 'DPP 478' (GRO residue count — wrong for pdb2gmx-assembled ITPs).\n"
            "Fix: use molecule_entries from protein_topology_manifest.json for protein entries;\n"
            "extract [ molecules ] from lipid_topol.top for lipid entries.\n"
            "Retrying will not fix this — the topology assembler code must be corrected."
        ),
    ),

    # ── Duplicate [ defaults ] in expanded topology ───────────────────────────
    _Signal(
        category    = ErrorCategory.TOPOLOGY_DUPLICATE_DEFAULTS_ERROR,
        severity    = ErrorSeverity.FATAL,
        patterns    = [
            r"TOPOLOGY_DUPLICATE_DEFAULTS_ERROR",
            r"Found a second defaults directive",
            r"expanded topology has \d+ \[ defaults \] section",
        ],
        confidence  = 0.99,
        explanation = (
            "GROMACS found more than one [ defaults ] section in the expanded topology.\n"
            "Cause: forcefield.itp is included more than once. Most common source:\n"
            "pdb2gmx writes an absolute path for the FF include "
            "(e.g. #include \"/abs/.../oplsaa_membrane.ff/forcefield.itp\"), which\n"
            "_walk_includes resolves locally and stores in molecule_itps. The\n"
            "assembler then adds a second canonical FF include on top of it.\n"
            "Fix: _walk_includes now skips .ff/ paths; assembler guards the loop;\n"
            "_clean_topology strips [ defaults ] from lipid sub-includes. Retrying\n"
            "will not help — this is a deterministic code-level bug."
        ),
    ),

    # ── Topology include resolution failure ───────────────────────────────────
    _Signal(
        category    = ErrorCategory.TOPOLOGY_INCLUDE_RESOLUTION_ERROR,
        severity    = ErrorSeverity.FATAL,
        patterns    = [
            r"broken topology includes",
            r"Cannot resolve .* \(in .*\.top\)",
            r"Cannot resolve .* \(in .*\.itp\)",
            r"TOPOLOGY INCLUDE ERROR",
        ],
        confidence  = 0.99,
        explanation = (
            "Topology include path resolution failed.\n"
            "Cause: the assembled topol.top contains a #include path that cannot be\n"
            "resolved — typically caused by an absolute path from the manifest being\n"
            "naively concatenated with a relative prefix.\n"
            "Fix: the topology assembler must normalize all include paths using\n"
            "os.path.relpath(abs_path, topology_dir) regardless of whether the\n"
            "manifest entry is absolute or relative."
        ),
    ),

    # ── OXT/terminus mismatch — pdb2gmx run on embedded/mixed GRO ───────────
    _Signal(
        category    = ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR,
        severity    = ErrorSeverity.FATAL,
        patterns    = [
            r"Atom OXT in residue .* was not found in rtp entry",
            r"was not found in rtp entry .* with \d+ atoms",
        ],
        confidence  = 0.97,
        explanation = (
            "OXT/terminus mismatch detected in topology preprocessing.\n"
            "Cause: pdb2gmx was run on an embedded GRO file where chain/TER "
            "semantics are lost — the C-terminal residue cannot be identified correctly.\n"
            "Fix: run pdb2gmx only on the original protein PDB "
            "(generate_protein_topology step), then build the combined topology "
            "with assemble_system_topology."
        ),
    ),

    # ── Salida con error sin categoría ────────────────────────────────────────
    _Signal(
        category    = ErrorCategory.NONZERO_EXIT,
        severity    = ErrorSeverity.RECOVERABLE,
        patterns    = [
            r"Fatal error",
            r"gmx.*error",
            r"Error in user input",
        ],
        confidence  = 0.70,
        explanation = "GROMACS terminó con error. Ver stderr para detalles.",
    ),
]


def _extract_evidence_lines(text: str, patterns: list[re.Pattern], n: int = 5) -> list[str]:
    """Extrae líneas de contexto alrededor de los matches para el log."""
    lines = text.splitlines()
    evidence: list[str] = []
    for i, line in enumerate(lines):
        for pat in patterns:
            if pat.search(line):
                start = max(0, i - 1)
                end   = min(len(lines), i + n)
                chunk = lines[start:end]
                evidence.extend(chunk)
                break
    # Deduplicar preservando orden
    seen: set[str] = set()
    return [l for l in evidence if not (l in seen or seen.add(l))]


# ═══════════════════════════════════════════════════════════════════════════════
# Capa 2 — Funciones de planificación por categoría
# ═══════════════════════════════════════════════════════════════════════════════

def _plan_lincs(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """LINCS warning/fatal → reducir dt, agregar posres si hace falta."""
    actions: list[RemediationAction] = []

    # 1. Detectar dt actual en el MDP
    mdp_files = list(step_dir.glob("*.mdp"))
    current_dt = "0.002"
    new_dt     = "0.001"

    for mdp in mdp_files:
        content = mdp.read_text()
        m = re.search(r"^\s*dt\s*=\s*([\d.]+)", content, re.MULTILINE)
        if m:
            current_dt = m.group(1)
            try:
                new_dt = str(round(float(current_dt) / 2, 4))
            except ValueError:
                pass
            break

    if mdp_files:
        actions.append(RemediationAction(
            action_type = ActionType.PATCH_MDP,
            description = f"Reducir timestep de {current_dt} a {new_dt} ps",
            target_file = mdp_files[0].name,
            patch_key   = "dt",
            patch_value = new_dt,
            patch_old   = current_dt,
            rationale   = "LINCS falla cuando los átomos se mueven demasiado en un paso. "
                          "Reducir dt a la mitad es la corrección estándar.",
            confidence  = 0.90,
        ))

        # 2. Reducir nstlist si es > 10
        for mdp in mdp_files:
            content = mdp.read_text()
            m = re.search(r"^\s*nstlist\s*=\s*(\d+)", content, re.MULTILINE)
            if m and int(m.group(1)) > 10:
                actions.append(RemediationAction(
                    action_type = ActionType.PATCH_MDP,
                    description = "Reducir nstlist a 10 para mayor estabilidad",
                    target_file = mdp.name,
                    patch_key   = "nstlist",
                    patch_value = "10",
                    patch_old   = m.group(1),
                    rationale   = "nstlist alto con dt pequeño puede causar inestabilidades.",
                    confidence  = 0.75,
                ))
            break

    # 3. Reset del step (limpiar outputs del intento fallido)
    actions.append(RemediationAction(
        action_type = ActionType.RESET_STEP,
        description = "Limpiar outputs del intento fallido antes de retry",
        rationale   = "Los archivos parciales del intento fallido pueden "
                      "interferir con el retry.",
        confidence  = 1.0,
    ))

    return RemediationPlan(
        step_id          = diag.step_id,
        diagnosis        = diag,
        actions          = actions,
        is_applicable    = True,
        requires_human   = False,
        max_retries      = 3,
        strategy         = "Reducir timestep + limpiar estado + retry",
        expected_outcome = "LINCS debería resolverse con dt más pequeño. "
                           "Si persiste con dt=0.0005, puede indicar clashes severos.",
    )


def _plan_nan_energy(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """NaN en energía / sistema explosivo → reducción agresiva de dt + minimización extra."""
    actions: list[RemediationAction] = []

    mdp_files = list(step_dir.glob("*.mdp"))

    if mdp_files:
        # Reducir dt agresivamente
        actions.append(RemediationAction(
            action_type = ActionType.PATCH_MDP,
            description = "Reducir dt a 0.0005 ps (modo ultra-conservador)",
            target_file = mdp_files[0].name,
            patch_key   = "dt",
            patch_value = "0.0005",
            rationale   = "NaN indica sistema explosivo. "
                          "Se necesita timestep muy pequeño para estabilizar.",
            confidence  = 0.88,
        ))

        # Si es equilibración: agregar restraints de posición temporales
        if diag.stage in ("equilibration", "production"):
            actions.append(RemediationAction(
                action_type = ActionType.INJECT_RESTRAINTS,
                description = "Agregar define=-DPOSRES para restraints de posición",
                target_file = mdp_files[0].name,
                patch_key   = "define",
                patch_value = "-DPOSRES",
                rationale   = "Los restraints de posición estabilizan el sistema "
                              "durante la fase inicial de equilibración.",
                confidence  = 0.82,
            ))

    actions.append(RemediationAction(
        action_type = ActionType.RESET_STEP,
        description = "Limpiar outputs del intento fallido",
        confidence  = 1.0,
        rationale   = "Partir limpio después de NaN.",
    ))

    return RemediationPlan(
        step_id          = diag.step_id,
        diagnosis        = diag,
        actions          = actions,
        is_applicable    = True,
        requires_human   = False,
        max_retries      = 2,
        strategy         = "Reducción agresiva de dt + restraints de posición + retry",
        expected_outcome = "Con dt=0.0005 y posres el sistema debería estabilizarse. "
                           "Si persiste: revisar la minimización energética.",
    )


def _plan_fmax(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """Fmax no convergió → más pasos, emtol más bajo, cambio de integrador."""
    actions: list[RemediationAction] = []

    mdp_files = list(step_dir.glob("*.mdp"))

    if mdp_files:
        mdp = mdp_files[0]
        content = mdp.read_text()

        # Detectar nsteps actual y doblar
        m_nsteps = re.search(r"^\s*nsteps\s*=\s*(\d+)", content, re.MULTILINE)
        if m_nsteps:
            current = int(m_nsteps.group(1))
            new_val = str(current * 2)
            actions.append(RemediationAction(
                action_type = ActionType.PATCH_MDP,
                description = f"Doblar nsteps de minimización: {current} → {new_val}",
                target_file = mdp.name,
                patch_key   = "nsteps",
                patch_value = new_val,
                patch_old   = str(current),
                rationale   = "Más pasos dan más oportunidad de convergencia.",
                confidence  = 0.88,
            ))

        # Bajar emtol si es > 100
        m_emtol = re.search(r"^\s*emtol\s*=\s*([\d.]+)", content, re.MULTILINE)
        if m_emtol:
            current_tol = float(m_emtol.group(1))
            if current_tol > 500:
                new_tol = str(current_tol / 2)
                actions.append(RemediationAction(
                    action_type = ActionType.PATCH_MDP,
                    description = f"Reducir emtol de {current_tol} a {new_tol} kJ/mol/nm",
                    target_file = mdp.name,
                    patch_key   = "emtol",
                    patch_value = new_tol,
                    patch_old   = str(current_tol),
                    rationale   = "emtol más bajo fuerza convergencia más estricta.",
                    confidence  = 0.80,
                ))

        # Cambiar integrador a l-bfgs si ya es steep
        m_integ = re.search(r"^\s*integrator\s*=\s*(\w[\w-]*)", content, re.MULTILINE)
        if m_integ and m_integ.group(1).strip() == "steep":
            actions.append(RemediationAction(
                action_type = ActionType.PATCH_MDP,
                description = "Cambiar integrador: steep → l-bfgs (más robusto)",
                target_file = mdp.name,
                patch_key   = "integrator",
                patch_value = "l-bfgs",
                patch_old   = "steep",
                rationale   = "l-bfgs es más eficiente para sistemas atascados "
                              "en mínimos locales.",
                confidence  = 0.78,
            ))

    actions.append(RemediationAction(
        action_type = ActionType.RESET_STEP,
        description = "Limpiar outputs para retry",
        confidence  = 1.0,
        rationale   = "Partir desde la entrada original.",
    ))

    return RemediationPlan(
        step_id          = diag.step_id,
        diagnosis        = diag,
        actions          = actions,
        is_applicable    = True,
        requires_human   = False,
        max_retries      = 3,
        strategy         = "Más pasos + emtol menor + integrador más robusto",
        expected_outcome = "La minimización debería converger con más pasos y l-bfgs.",
    )


def _plan_poor_equilibration(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """Temperatura/presión inestable → tau_t más conservador + posres."""
    actions: list[RemediationAction] = []

    mdp_files = list(step_dir.glob("*.mdp"))

    if mdp_files:
        for mdp in mdp_files:
            content = mdp.read_text()

            m_tau = re.search(r"^\s*tau_t\s*=\s*([\d.\s]+)", content, re.MULTILINE)
            if m_tau:
                # tau_t puede ser múltiples valores (por grupo de termostato)
                vals = m_tau.group(1).strip().split()
                new_vals = []
                for v in vals:
                    try:
                        new_vals.append(str(round(float(v) * 2, 2)))
                    except ValueError:
                        new_vals.append(v)
                actions.append(RemediationAction(
                    action_type = ActionType.PATCH_MDP,
                    description = f"Doblar tau_t: {m_tau.group(1).strip()} → {' '.join(new_vals)}",
                    target_file = mdp.name,
                    patch_key   = "tau_t",
                    patch_value = " ".join(new_vals),
                    patch_old   = m_tau.group(1).strip(),
                    rationale   = "tau_t más grande amortigua las fluctuaciones de temperatura.",
                    confidence  = 0.82,
                ))
            break

    actions.append(RemediationAction(
        action_type = ActionType.RESET_STEP,
        description = "Limpiar outputs para retry",
        confidence  = 1.0,
        rationale   = "Partir desde checkpoint inicial.",
    ))

    return RemediationPlan(
        step_id          = diag.step_id,
        diagnosis        = diag,
        actions          = actions,
        is_applicable    = True,
        requires_human   = False,
        max_retries      = 2,
        strategy         = "tau_t más conservador + retry",
        expected_outcome = "La temperatura debería estabilizarse con tau_t más grande.",
    )


def _plan_artifact_validation(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """exit=0 but expected outputs missing — non-retryable, needs topology audit."""
    return RemediationPlan(
        step_id        = diag.step_id,
        diagnosis      = diag,
        actions        = [RemediationAction(
            action_type  = ActionType.LOG_ONLY,
            description  = "Fatal: pdb2gmx exited 0 but required output files are missing",
            rationale    = (
                "GROMACS succeeded but the expected artifact is absent. "
                "For multichain proteins, pdb2gmx generates chain-specific files "
                "(topol_Protein_chain_A.itp / posre_Protein_chain_A.itp) instead of "
                "a literal posre.itp. SimForge must validate the include graph "
                "dynamically, not assume fixed filenames."
            ),
            confidence   = 1.0,
            is_reversible = True,
        )],
        is_applicable  = False,
        requires_human = True,
        max_retries    = 0,
        strategy       = "No retry — inspect expected_outputs in metadata.json",
        expected_outcome = (
            "Remove literal posre.itp from expected_outputs; use "
            "protein_topology_manifest.json to discover chain topology files."
        ),
    )


def _plan_topology_position_restraint_index(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """[ position_restraints ] atom index exceeds local moleculetype atom count — non-retryable."""
    return RemediationPlan(
        step_id        = diag.step_id,
        diagnosis      = diag,
        actions        = [RemediationAction(
            action_type  = ActionType.LOG_ONLY,
            description  = "Fatal: [ position_restraints ] contains atom index out of bounds for moleculetype",
            rationale    = (
                "GROMACS found an atom index in a [ position_restraints ] section that exceeds "
                "the local atom count of the enclosing moleculetype. "
                "This means global system atom indices were written instead of per-chain local indices. "
                "Fix: split strong_posre.itp into per-chain files (strong_posre_Protein_chain_A.itp, …) "
                "and write only local 1-based indices from each chain's [ atoms ] section. "
                "Retrying will not help — the restraint generator must produce local indices."
            ),
            confidence   = 1.0,
            is_reversible = True,
        )],
        is_applicable  = False,
        requires_human = True,
        max_retries    = 0,
        strategy       = "No retry — restraint generator must use local [ atoms ] indices per chain",
        expected_outcome = (
            "Per-chain strong_posre_Protein_chain_X.itp files with local indices will pass grompp "
            "without 'Atom index out of bounds'. "
            "The pre-grompp validator TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR will catch violations early."
        ),
    )


def _plan_topology_undefined_moleculetype(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """[ molecules ] entry has no [ moleculetype ] definition — non-retryable assembler bug."""
    return RemediationPlan(
        step_id        = diag.step_id,
        diagnosis      = diag,
        actions        = [RemediationAction(
            action_type  = ActionType.LOG_ONLY,
            description  = "Fatal: [ molecules ] entry has no matching [ moleculetype ] definition",
            rationale    = (
                "GROMACS cannot find a moleculetype for one or more entries in [ molecules ]. "
                "Common causes: (1) multichain protein — pdb2gmx generates 'Protein_chain_A 1' etc. "
                "but SimForge wrote 'Protein 1'; (2) lipid — pdb2gmx generated 'DPP 1' (combined "
                "topology covering all lipids) but SimForge counted GRO residues and wrote 'DPP 478'. "
                "Retrying will not help — the topology assembler must use pdb2gmx-authoritative entries."
            ),
            confidence   = 1.0,
            is_reversible = True,
        )],
        is_applicable  = False,
        requires_human = True,
        max_retries    = 0,
        strategy       = "No retry — topology assembler must use molecule_entries from manifest and lipid [ molecules ] from pdb2gmx",
        expected_outcome = (
            "Use molecule_entries from protein_topology_manifest.json (not hardcoded 'Protein 1'). "
            "Use [ molecules ] from lipid_topol.top (not GRO residue count). "
            "The _validate_moleculetypes check will catch mismatches before grompp is called."
        ),
    )


def _plan_topology_duplicate_defaults(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """Duplicate [ defaults ] in expanded topology — non-retryable FF double-include bug."""
    return RemediationPlan(
        step_id        = diag.step_id,
        diagnosis      = diag,
        actions        = [RemediationAction(
            action_type  = ActionType.LOG_ONLY,
            description  = "Fatal: forcefield.itp included more than once in assembled topology",
            rationale    = (
                "GROMACS found a second [ defaults ] section, meaning forcefield.itp "
                "was expanded twice. Root cause: pdb2gmx writes the FF include as an "
                "absolute path (e.g. #include \"/abs/.../forcefield.itp\"). "
                "_walk_includes resolves this absolute path and places it in molecule_itps "
                "alongside chain .itp files. The assembler then adds its own canonical FF "
                "include at the top, creating a duplicate. "
                "The fix requires updating _walk_includes to classify .ff/ paths as FF "
                "includes and guarding the molecule_itps loop in assemble_system_topology. "
                "Retrying will not fix this."
            ),
            confidence   = 1.0,
            is_reversible = True,
        )],
        is_applicable  = False,
        requires_human = True,
        max_retries    = 0,
        strategy       = "No retry — topology assembler code must be fixed",
        expected_outcome = (
            "The assembled topol.top must have exactly one #include for forcefield.itp. "
            "Ensure _walk_includes skips .ff/ paths and the molecule_itps loop filters "
            "out any entry with '.ff/' in its path."
        ),
    )


def _plan_topology_include_resolution(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """Broken #include path in assembled topology — non-retryable path normalization bug."""
    return RemediationPlan(
        step_id        = diag.step_id,
        diagnosis      = diag,
        actions        = [RemediationAction(
            action_type  = ActionType.LOG_ONLY,
            description  = "Fatal: topology include path cannot be resolved",
            rationale    = (
                "The assembled topol.top has a #include path that GROMACS cannot "
                "find. Typical cause: the protein_topology_manifest.json stores "
                "absolute paths, which are then naively concatenated with a "
                "relative prefix like '../01_generate_protein_topology/'. "
                "Retrying will not fix this — the assembler code must normalize "
                "every include via Path.resolve() + os.path.relpath()."
            ),
            confidence   = 1.0,
            is_reversible = True,
        )],
        is_applicable  = False,
        requires_human = True,
        max_retries    = 0,
        strategy       = "No retry — topology assembler path normalization must be fixed in code",
        expected_outcome = (
            "All includes in topol.top must be valid relative paths from the "
            "topology file's directory. Never concatenate a relative prefix with "
            "an absolute path from the manifest."
        ),
    )


def _plan_topology_preprocessing(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """OXT/terminus mismatch — deterministic preprocessing error, non-retryable."""
    return RemediationPlan(
        step_id        = diag.step_id,
        diagnosis      = diag,
        actions        = [RemediationAction(
            action_type  = ActionType.LOG_ONLY,
            description  = "Fatal: OXT/terminus mismatch — pdb2gmx run on mixed/embedded GRO",
            rationale    = (
                "This error is deterministic and cannot be fixed by retrying. "
                "pdb2gmx must run only on the original protein PDB "
                "(inputs/protein_1.pdb via generate_protein_topology), not on "
                "embed_in_bilayer/system.gro. Regenerate the workspace."
            ),
            confidence   = 1.0,
            is_reversible = True,
        )],
        is_applicable  = False,
        requires_human = True,
        max_retries    = 0,
        strategy       = "No retry — Phase 11 topology architecture required",
        expected_outcome = (
            "Use generate_protein_topology (pdb2gmx on protein PDB) + "
            "assemble_system_topology to build the correct combined topology."
        ),
    )


def _plan_missing_parameter(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """Parámetro FF faltante → fatal, no se puede remediar automáticamente."""
    return RemediationPlan(
        step_id        = diag.step_id,
        diagnosis      = diag,
        actions        = [RemediationAction(
            action_type = ActionType.LOG_ONLY,
            description = "Error fatal: parámetro de forcefield faltante",
            rationale   = "Este error requiere reparametrización manual del ligando. "
                          "SimForge no puede generar parámetros FF automáticamente.",
            confidence  = 1.0,
            is_reversible = True,
        )],
        is_applicable  = False,
        requires_human = True,
        max_retries    = 0,
        strategy       = "Intervención manual requerida",
        expected_outcome = "Revisar la parametrización del ligando en ParamChem/CGenFF "
                           "y re-ejecutar el step de parametrización.",
    )


def _plan_checkpoint(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """Checkpoint corrupto → eliminar .cpt y reiniciar desde input."""
    actions: list[RemediationAction] = []

    cpt_files = list(step_dir.glob("*.cpt"))
    for cpt in cpt_files:
        actions.append(RemediationAction(
            action_type = ActionType.DELETE_FILE,
            description = f"Eliminar checkpoint corrupto: {cpt.name}",
            target_file = cpt.name,
            rationale   = "Checkpoint corrupto → reiniciar desde el input original.",
            confidence  = 0.95,
        ))

    actions.append(RemediationAction(
        action_type = ActionType.RESET_STEP,
        description = "Limpiar estado para retry sin checkpoint",
        confidence  = 1.0,
        rationale   = "Sin checkpoint el executor usará el input original.",
    ))

    return RemediationPlan(
        step_id          = diag.step_id,
        diagnosis        = diag,
        actions          = actions,
        is_applicable    = True,
        requires_human   = False,
        max_retries      = 2,
        strategy         = "Eliminar checkpoint corrupto + retry desde input",
        expected_outcome = "El step debería correr desde el principio sin el .cpt corrupto.",
    )


def _plan_interactive_block(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """Proceso bloqueado esperando TTY → parchear script para eliminar la flag interactiva."""
    actions: list[RemediationAction] = []

    for script in sorted(step_dir.glob("*.sh")):
        try:
            content = script.read_text()
        except Exception:
            continue
        if "-princ" in content:
            actions.append(RemediationAction(
                action_type = ActionType.PATCH_SCRIPT,
                description = f"Eliminar -princ de {script.name} (bloquea esperando TTY)",
                target_file = script.name,
                patch_key   = "-princ",
                patch_value = "",
                patch_old   = "-princ",
                rationale   = "gmx editconf -princ solicita selección interactiva de grupo. "
                              "Sin terminal disponible el proceso queda bloqueado indefinidamente.",
                confidence  = 0.97,
            ))

    actions.append(RemediationAction(
        action_type = ActionType.RESET_STEP,
        description = "Limpiar outputs parciales antes de retry",
        confidence  = 1.0,
        rationale   = "El step fue terminado a la fuerza; limpiar estado antes de retry.",
    ))

    return RemediationPlan(
        step_id          = diag.step_id,
        diagnosis        = diag,
        actions          = actions,
        is_applicable    = bool(actions),
        requires_human   = not bool(actions),
        max_retries      = 1,
        strategy         = "Eliminar flag interactiva del script + retry",
        expected_outcome = "Sin -princ el editconf corre sin prompt. "
                           "La caja puede ser ligeramente mayor para péptidos elongados.",
    )


def _plan_generic_retry(diag: DiagnosisResult, step_dir: Path) -> RemediationPlan:
    """Fallback: reset + retry sin modificar nada."""
    return RemediationPlan(
        step_id        = diag.step_id,
        diagnosis      = diag,
        actions        = [RemediationAction(
            action_type = ActionType.RESET_STEP,
            description = "Reset limpio y retry sin modificaciones",
            rationale   = "Error no clasificado — intentar retry limpio antes "
                          "de escalar a intervención humana.",
            confidence  = 0.60,
        )],
        is_applicable  = True,
        requires_human = False,
        max_retries    = 1,
        strategy       = "Retry limpio",
        expected_outcome = "Si el error fue transitorio (I/O, red) debería resolverse. "
                           "Si persiste, escalar a revisión manual.",
    )


# Tabla de dispatch: categoría → función de planificación
_PLANNERS: dict[ErrorCategory, Callable] = {
    ErrorCategory.ARTIFACT_VALIDATION_ERROR:               _plan_artifact_validation,
    ErrorCategory.TOPOLOGY_POSITION_RESTRAINT_INDEX_ERROR: _plan_topology_position_restraint_index,
    ErrorCategory.TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR:   _plan_topology_undefined_moleculetype,
    ErrorCategory.TOPOLOGY_DUPLICATE_DEFAULTS_ERROR:       _plan_topology_duplicate_defaults,
    ErrorCategory.TOPOLOGY_INCLUDE_RESOLUTION_ERROR:  _plan_topology_include_resolution,
    ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR:       _plan_topology_preprocessing,
    ErrorCategory.LINCS_WARNING:       _plan_lincs,
    ErrorCategory.LINCS_FATAL:         _plan_lincs,
    ErrorCategory.NAN_ENERGY:          _plan_nan_energy,
    ErrorCategory.EXPLODING_SYSTEM:    _plan_nan_energy,
    ErrorCategory.FMAX_NOT_CONVERGED:  _plan_fmax,
    ErrorCategory.POOR_EQUILIBRATION:  _plan_poor_equilibration,
    ErrorCategory.MISSING_PARAMETER:   _plan_missing_parameter,
    ErrorCategory.ATOM_TYPE_MISMATCH:  _plan_missing_parameter,
    ErrorCategory.CORRUPT_CHECKPOINT:  _plan_checkpoint,
    ErrorCategory.INTERACTIVE_BLOCK:   _plan_interactive_block,
    ErrorCategory.CHARGE_IMBALANCE:    _plan_generic_retry,
    ErrorCategory.MISSING_INPUT_FILE:  _plan_generic_retry,
    ErrorCategory.TIMEOUT:             _plan_generic_retry,
    ErrorCategory.NONZERO_EXIT:        _plan_generic_retry,
    ErrorCategory.MISSING_OUTPUT:      _plan_generic_retry,
    ErrorCategory.UNKNOWN:             _plan_generic_retry,
}


# ═══════════════════════════════════════════════════════════════════════════════
# AdaptiveReasoner — API pública
# ═══════════════════════════════════════════════════════════════════════════════

class AdaptiveReasoner:
    """
    Motor de diagnóstico y planificación de remediación.

    Uso:
        reasoner = AdaptiveReasoner()

        # Después de que un step falla:
        diag = reasoner.diagnose(record, step_dir)
        plan = reasoner.plan_remediation(diag, step_dir)

        # Aplicar el plan:
        executor.apply(plan)
    """

    def diagnose(
        self,
        record:   StepExecutionRecord,
        step_dir: Path,
    ) -> DiagnosisResult:

        # ── Fuente primaria: diagnóstico rico ──────────────────────────
        rich_diag = getattr(record, "gromacs_diagnostic", None)

        if rich_diag is not None:
            return self.diagnose_from_gromacs(
                rich_diag,
                record,
            )

        # ── Texto a analizar ───────────────────────────────────────────
        combined = (record.stdout or "") + "\n" + (record.stderr or "")

        # ── Metadata del step ──────────────────────────────────────────
        meta_file = step_dir / "metadata.json"
        stage  = ""
        engine = ""

        if meta_file.exists():
            import json
            try:
                meta   = json.loads(meta_file.read_text())
                stage  = meta.get("stage", "")
                engine = meta.get("engine", "")
            except Exception:
                pass

        # ── Evaluar señales ───────────────────────────────────────────────────
        matched_signal: _Signal | None = None
        primary_line   = ""

        for signal in _SIGNAL_PATTERNS:
            hit, line = signal.match(combined)
            if hit:
                matched_signal = signal
                primary_line   = line
                break

        # ── Construir diagnóstico ─────────────────────────────────────────────
        if matched_signal:
            evidence = _extract_evidence_lines(
                combined,
                matched_signal.patterns,
                n = 5,
            )
            return DiagnosisResult(
                step_id        = record.step_id,
                step_dir       = str(step_dir),
                category       = matched_signal.category,
                severity       = matched_signal.severity,
                confidence     = matched_signal.confidence,
                primary_signal = primary_line,
                evidence_lines = evidence[:10],   # máximo 10 líneas
                exit_code      = record.exit_code,
                stage          = stage,
                engine         = engine,
                explanation    = matched_signal.explanation,
                reasoning      = (
                    f"Señal '{matched_signal.category.value}' detectada por "
                    f"patrón regex en stdout/stderr. "
                    f"Confianza: {matched_signal.confidence:.0%}."
                ),
            )

        # ── Fallback: outputs faltantes ────────────────────────────────────────
        if record.outputs_missing:
            return DiagnosisResult(
                step_id        = record.step_id,
                step_dir       = str(step_dir),
                category       = ErrorCategory.MISSING_OUTPUT,
                severity       = ErrorSeverity.RECOVERABLE,
                confidence     = 0.80,
                primary_signal = f"Outputs faltantes: {record.outputs_missing}",
                evidence_lines = [],
                exit_code      = record.exit_code,
                stage          = stage,
                engine         = engine,
                explanation    = "El step terminó sin error explícito pero "
                                 "no generó los archivos esperados.",
                reasoning      = "Diagnóstico por ausencia de outputs esperados.",
            )

        # ── Fallback: exit code ≠ 0 sin señal ─────────────────────────────────
        return DiagnosisResult(
            step_id        = record.step_id,
            step_dir       = str(step_dir),
            category       = ErrorCategory.UNKNOWN,
            severity       = ErrorSeverity.RECOVERABLE,
            confidence     = 0.50,
            primary_signal = f"Exit code: {record.exit_code}",
            evidence_lines = combined.splitlines()[-10:],
            exit_code      = record.exit_code,
            stage          = stage,
            engine         = engine,
            explanation    = "Error no identificado. Ver stdout/stderr completo.",
            reasoning      = "Ninguna señal conocida hizo match. Fallback a UNKNOWN.",
        )

    def plan_remediation(
        self,
        diagnosis: DiagnosisResult,
        step_dir:  Path,
    ) -> RemediationPlan:
        """
        Genera un RemediationPlan desde un DiagnosisResult.

        Dispatch por categoría a la función de planificación correspondiente.
        Si la categoría no tiene planner registrado → fallback genérico.
        """
        planner = _PLANNERS.get(diagnosis.category, _plan_generic_retry)
        return planner(diagnosis, step_dir)