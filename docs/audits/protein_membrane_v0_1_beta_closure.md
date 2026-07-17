# Cierre de fase — Protein–Membrane Builder v0.1-beta

**Fecha:** 2026-07-17
**Decisión:** cerrar la fase de validación geométrica estática como beta cualitativa.

La generación automática fue revisada en un receptor tipo GLP-1R y en un canal grande. En ambos casos produce sistemas iniciales visualmente razonables: el receptor ya no presenta el gran cráter anular artificial de iteraciones anteriores y el canal conserva tanto su envolvente lipídica externa como su lumen funcional abierto.

El resultado no se declara listo para producción desatendida. Es adecuado para validación dinámica corta y controlada; antes de MD larga siguen siendo obligatorios el QC visual y las métricas de membrana, proteína, solvente y consistencia topológica.

La fase corrigió la paridad de cutoffs de InflateGRO (1.4 nm inicial y 0.0 nm durante shrink), la sincronización coordenadas–topología, la exclusión pre-shrink consciente de TM/slab, la clasificación de lípidos atrapados y la integración conservadora de reparación anular. También estableció que los poros centrales legítimos de canales no deben clasificarse como vacíos anulares.

La siguiente fase evaluará estabilidad dinámica mediante minimización, NVT/NPT restringidos y equilibración corta, sin continuar parcheando geometría estática salvo defecto bloqueante.
