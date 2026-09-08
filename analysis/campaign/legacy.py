"""Read-only legacy annotations on campaign SystemRecords, never executors.

Residue templates are comparison references, never generated selections.
"""
from pathlib import Path
import csv
import math
import re

from analysis.campaign.structure.index_groups import parse_index_groups, resolve_index_group

REFERENCES = {
    'ActiveSite_AA': [58, 59, 62, 98, 101, 165, 195, 197, 233, 299],
    'Catalytic_AA': [197, 233, 300],
    'ActiveSite_HMG': [559, 655, 656, 690, 691, 692, 751, 755, 866],
    'Catalytic_HMG': [559, 691, 692, 866],
    'ActiveSite_LIP': [75, 77, 78, 79, 80, 85, 252, 255, 256, 258, 259, 260, 263, 264],
    'Catalytic_LP': [152, 176, 263],
    'ActiveSite_AG': [227, 228, 229, 232, 235, 236, 550, 551, 552, 557, 559],
    'Catalytic_AG': [228, 235],
}
TARGET_GROUPS = {'AA': ('ActiveSite_AA', 'Catalytic_AA'),
                 'AG': ('ActiveSite_AG', 'Catalytic_AG'),
                 'HMG-R': ('ActiveSite_HMG', 'Catalytic_HMG'),
                 'LP': ('ActiveSite_LIP', 'Catalytic_LP')}
THESIS_CORE = ['protein_rmsd', 'ligand_rmsd', 'active_site_mindist', 'active_site_contacts', 'catalytic_com_distance']
OPTIONAL_SHORT = ['protein_rmsf', 'sasa', 'hbonds', 'interaction_energy', 'mmpbsa', 'mmpbsa_decomposition']
SHORT = THESIS_CORE
MECHANISTIC = THESIS_CORE + OPTIONAL_SHORT + ['radius_of_gyration', 'pca', 'fel', 'clustering', 'dccm']
PATTERNS = {
    'protein_rmsd': [('rmsd_protein.xvg', 'rmsd')],
    'ligand_rmsd': [('rmsd_ligand.xvg', 'rmsd')],
    'protein_rmsf': [('rmsf*.xvg', 'rmsf|rms fluctuation')], 'sasa': [('sasa*.xvg', 'area|sasa')],
    'hbonds': [('hbonds*.xvg', 'hydrogen|h.?bond')],
    'active_site_mindist': [('mindist_lig_active.xvg', 'distance')],
    'active_site_contacts': [('contacts_lig_active.xvg', 'contact')],
    'catalytic_com_distance': [('dist_lig_catalytic.xvg', 'distance')],
    'radius_of_gyration': [('gyrate*.xvg', 'gyration'), ('rg*.xvg', 'gyration')],
    'interaction_energy': [('interaction_energy*.xvg', r'(coul|lj).*?(protein.*lig|lig.*protein)'),
                           ('energy_interaction*.xvg', r'(coul|lj).*?(protein.*lig|lig.*protein)')],
    'clustering': [('cluster*.xvg', r'cluster.*(size|population|number|distribution)'),
                   ('clusters*.xpm', 'cluster')],
    'dccm': [('dccm*.xpm', 'cross.?correlation|dccm')],
    'mmpbsa': [('binding_energy.xvg', 'binding.*energy|mmpbsa')],
    'mm_component': [('energy_MM.xvg', 'energy')],
    'polar_component': [('polar.xvg', 'energy|solvation')],
    'apolar_component': [('apolar.xvg', 'energy|solvation')],
    'pca_eigenvalues': [('eigenval*.xvg', 'eigenvalue')],
    'pca_projection': [('proj*.xvg', 'projection'), ('pca_proj*.xvg', 'projection')],
    'fel': [('fel*.xpm', 'free energy|gibbs'), ('gibbs*.xpm', 'free energy|gibbs')],
}


def _valid_xvg(path, header_pattern):
    header, rows, width = [], 0, None
    try:
        with path.open() as stream:
            for line in stream:
                line = line.strip()
                if not line:
                    continue
                # Command comments do not establish the meaning of the data.
                if line.startswith('#'):
                    continue
                if line.startswith('@'):
                    header.append(line)
                    continue
                values = [float(x) for x in line.split()]
                if len(values) < 2 or not all(math.isfinite(x) for x in values):
                    return False
                if width is not None and width != len(values):
                    return False
                width = len(values)
                rows += 1
        return rows >= 2 and bool(re.search(header_pattern, '\n'.join(header), re.I))
    except (OSError, ValueError):
        return False


def _valid_xpm(path, header_pattern):
    """Require a complete XPM matrix, not only a title and dimensions."""
    try:
        content = path.read_text(errors='replace')
        title = re.search(r'/\*\s*title:\s*"([^"]*)"', content, re.I)
        if 'XPM' not in content or not re.search(r'}\s*;', content) or not title or not re.search(header_pattern, title[1], re.I):
            return False
        lines = re.findall(r'^\s*"([^"\n]*)"', content, re.M)
        dimensions = next((i for i, row in enumerate(lines)
                           if re.fullmatch(r'\s*\d+\s+\d+\s+\d+\s+\d+\s*', row)), None)
        if dimensions is None:
            return False
        width, height, colors, cpp = map(int, lines[dimensions].split())
        if min(width, height, colors, cpp) <= 0:
            return False
        payload = lines[dimensions + 1:]
        if len(payload) != colors + height:
            return False
        keys = {row[:cpp] for row in payload[:colors]
                if len(row) >= cpp and re.search(r'\s+c\s+', row[cpp:])}
        return len(keys) == colors and all(len(row) == width * cpp and
            all(row[i:i + cpp] in keys for i in range(0, len(row), cpp))
            for row in payload[colors:])
    except (OSError, ValueError):
        return False


def existing_results(directory):
    """Conservative local result evidence; component energies alone are not PBSA."""
    directory = Path(directory)
    found, rejected = {}, []
    for family, patterns in PATTERNS.items():
        for pattern, header in patterns:
            for path in sorted(directory.glob(pattern)):
                if not path.is_file():
                    continue
                if path.suffix == '.xpm':
                    valid = _valid_xpm(path, header)
                else:
                    valid = _valid_xvg(path, header)
                    if valid:
                        metadata = '\n'.join(line for line in path.read_text().splitlines()
                                             if line.lstrip().startswith('@'))
                        # Related distributions/correlations are not the requested series.
                        if family in {'protein_rmsd', 'ligand_rmsd', 'protein_rmsf',
                                      'sasa', 'hbonds', 'radius_of_gyration'}:
                            valid = not re.search(r'autocorrelation|distribution|histogram', metadata, re.I)
                        if family == 'protein_rmsf':
                            valid &= not re.search(r'ligand|\bLIG\b', metadata, re.I)
                    if valid and family == 'interaction_energy':
                        # A single Coulomb/LJ component is only partial evidence.
                        metadata = '\n'.join(line for line in path.read_text().splitlines()
                                             if line.lstrip().startswith('@'))
                        pair = r'(protein.*lig|lig.*protein)'
                        valid = all(re.search(term + r'[^\n]*' + pair, metadata, re.I)
                                    for term in (r'coul', r'lj'))
                if valid:
                    found.setdefault(family, []).append(str(path))
                else:
                    rejected.append(str(path))
    for filename, family in [('energy_summary.csv', 'mmpbsa'),
                             ('residues_energy_summary.csv', 'mmpbsa_decomposition')]:
        path = directory / filename
        if not path.is_file():
            continue
        try:
            with path.open() as stream:
                rows = list(csv.reader(stream))
            header = [x.strip().lower() for x in rows[0]]
            data = [r for r in rows[1:] if r and any(x.strip() for x in r)]
            if family == 'mmpbsa':
                valid = {'energy', 'average'} <= set(header)
                label_column = header.index('energy') if 'energy' in header else 0
                valid &= {'total', 'polar-solvation', 'non-polar-solvation'} <= {
                    r[label_column].strip().lower() for r in data}
            else:
                valid = {'residue', 'total', 'polar', 'apolar'} <= set(header)
                label_column = header.index('residue') if 'residue' in header else 0
            labels = [r[label_column].strip().lower() for r in data]
            valid &= len(labels) == len(set(labels))
            valid &= len([h for h in header if h]) == len(set(h for h in header if h))
            numeric_columns = [i for i, name in enumerate(header) if i != label_column and name]
            valid &= bool(data) and bool(numeric_columns) and all(
                len(r) == len(header) and bool(r[label_column].strip()) and all(math.isfinite(float(r[i]))
                for i in numeric_columns) for r in data)
        except (OSError, ValueError, IndexError):
            valid = False
        if valid:
            found.setdefault(family, []).append(str(path))
        else:
            rejected.append(str(path))
    if 'pca_eigenvalues' in found and 'pca_projection' in found:
        eigenvalues = sorted(set(found['pca_eigenvalues']))
        projections = sorted(set(found['pca_projection']))
        # Do not combine alternative runs, or explicitly conflicting analysis windows.
        if len(eigenvalues) == len(projections) == 1:
            windows = []
            for filename in eigenvalues + projections:
                commands = '\n'.join(line for line in Path(filename).read_text().splitlines()
                                     if line.lstrip().startswith('#'))
                windows.append({flag: float(value) for flag, value in
                                re.findall(r'(?<!\S)-(b|e)\s+(-?\d+(?:\.\d+)?)', commands)})
            if all(windows[0][flag] == windows[1][flag]
                   for flag in windows[0].keys() & windows[1].keys()):
                found['pca'] = eigenvalues + projections
    return {family: sorted(set(paths)) for family, paths in sorted(found.items())}, sorted(set(rejected))


def annotate_legacy(rec, directory):
    """Return deterministic evidence without modifying files or invoking GROMACS."""
    directory = Path(directory)
    evidence, issues = [], []
    def fact(value, confidence, source):
        return {'value': value, 'confidence': confidence, 'provenance': source}

    candidates = sorted(directory.glob('*.ndx'))
    # Prefer a local index scope, never a sibling's index or an arbitrary first hit.
    if not candidates and rec.index_path:
        candidates = sorted(Path(rec.index_path).parent.glob('*.ndx'))
    index = candidates[0] if len(candidates) == 1 else None
    if len(candidates) != 1:
        issues.append('multiple index candidates' if candidates else 'missing index')
    groups = []
    if index:
        try:
            groups = parse_index_groups(index)
        except (OSError, ValueError) as exc:
            issues.append(f'invalid index: {exc}')
    names = {g.name for g in groups}
    atoms = []
    if rec.structure_path and Path(rec.structure_path).suffix == '.gro':
        from utils.gro_parser import parse_gro
        try:
            atoms = parse_gro(rec.structure_path).atoms
        except (OSError, ValueError, IndexError) as exc:
            issues.append(f'structure unavailable: {exc}')
    resnames = {a.residue_name.upper() for a in atoms}
    structural_tokens = set(resnames) | names
    compound_hits = set()
    topology = directory / 'topol.top'
    topology_content = topology.read_text(errors='replace') if topology.is_file() else ''
    includes = {Path(name).name for name in re.findall(r'#include\s+"([^"]+)"',
                '\n'.join(line.split(';', 1)[0] for line in topology_content.splitlines()))}
    # An unused parameter template does not prove a molecule is in the system.
    for path in sorted(directory.glob('*.itp')):
        if topology.is_file() and path.name not in includes:
            continue
        if not topology.is_file() and 'LIG' not in names:
            continue
        if path.stem in {f'A{i}' for i in range(1, 7)}:
            compound_hits.add(path.stem)
            evidence.append(f'{path.name}: compound identifier; ' +
                            ('included by topol.top' if topology.is_file() else 'corroborated by LIG index group'))
    section = ''
    for line in topology_content.splitlines():
        line = line.split(';', 1)[0].strip()
        if line.startswith('['):
            section = line.strip('[] ').lower()
        elif section == 'molecules' and line and not line.startswith('#'):
            fields = line.split()
            try:
                if len(fields) >= 2 and int(fields[1]) > 0:
                    structural_tokens.add(fields[0])
            except ValueError:
                issues.append('invalid topology molecule count')
    report = directory / 'assembly_report.yaml'
    if report.is_file():
        import yaml
        try:
            data = yaml.safe_load(report.read_text()) or {}
            ligand = data.get('ligand_mol_name')
            if ligand:
                structural_tokens.add(str(ligand))
            ligand_path = data.get('outputs', {}).get('ligand_itp', '')
            if Path(ligand_path).stem in {f'A{i}' for i in range(1, 7)}:
                compound_hits.add(Path(ligand_path).stem)
                evidence.append(f'{report.name}: outputs.ligand_itp={ligand_path}')
        except (OSError, ValueError, TypeError, AttributeError, yaml.YAMLError):
            issues.append('invalid assembly report')
    structural_compounds = structural_tokens & {f'A{i}' for i in range(1, 7)}
    compound_hits |= structural_compounds
    if structural_compounds:
        evidence.append(f'{directory}: index/residue/topology molecular names {sorted(structural_compounds)}')
    hinted = set(re.findall(r'(?<![A-Za-z0-9])A[1-6](?![A-Za-z0-9])', directory.name))
    compound = fact(next(iter(compound_hits)), .9, evidence.copy()) if len(compound_hits) == 1 else fact(None, 0, [])
    if not compound_hits and len(hinted) == 1:
        compound = fact(next(iter(hinted)), .45, [f'folder hint: {directory.name}'])
    if len(compound_hits) > 1 or (compound_hits and hinted and compound_hits != hinted):
        issues.append('conflicting compound identity')
    targets = {t for t, gs in TARGET_GROUPS.items() if names.intersection(gs)}
    target_hints = {t for t, pattern in {'AA': r'\bAA\b', 'AG': r'\bAG\b',
        'LP': r'\bLP\b|lipase', 'HMG-R': r'HMG'}.items()
        if re.search(pattern, directory.name.replace('_', '-'), re.I)}
    target = fact(next(iter(targets)), .85, [f'{index}: target-specific headers {sorted(names.intersection(TARGET_GROUPS[next(iter(targets))]))}']) if len(targets) == 1 else fact(None, 0, [])
    if not targets and len(target_hints) == 1:
        target = fact(next(iter(target_hints)), .45, [f'folder hint: {directory.name}'])
    if len(targets) > 1 or (targets and target_hints and targets != target_hints):
        issues.append('conflicting target identity')
    cofactor_names = sorted({x for x in structural_tokens if x.upper() in {'COA', 'HMG', 'NAD', 'NADP', 'NAP', 'FAD'}})
    cofactor = fact(cofactor_names, .95 if cofactor_names else 0,
                    [f'{directory}: index/residue/topology molecular names {cofactor_names}'] if cofactor_names else [])
    from analysis.campaign.legacy_duration import duration_evidence
    duration, duration_issues = duration_evidence(rec, directory)
    issues.extend(duration_issues)
    roles = {'protein': ('Protein', ()), 'backbone': ('Backbone', ()),
             'ligand': ('LIG', tuple(sorted(structural_tokens & {f'A{i}' for i in range(1, 7)}))),
             'cofactor': ('COA', tuple(cofactor_names))}
    sites = TARGET_GROUPS.get(target['value'])
    roles['active_site'] = (sites[0] if sites else 'ActiveSite', ('ActiveSite',))
    roles['catalytic_site'] = (sites[1] if sites else 'Catalytic', ('Catalytic',))
    mapping = {}
    selection_issues = set()
    for role, (name, aliases) in roles.items():
        try:
            group = resolve_index_group(index, name, aliases) if index and groups else None
        except ValueError as exc:
            issues.append(str(exc))
            selection_issues.add(str(exc))
            group = None
        entry = fact(group.name if group else None, .95 if group else 0,
                     [f'{index}: header position {group.group_id}; semantic name match'] if group else [])
        entry['group_id'] = group.group_id if group else None
        if group:
            entry['atom_count'] = len(group.atom_ids)
            if not group.atom_ids or (atoms and max(group.atom_ids) > len(atoms)):
                issues.append(f'{role}: empty or out-of-range atoms')
                selection_issues.add(f'{role}: empty or out-of-range atoms')
                entry['confidence'] = 0
            elif atoms:
                residues = sorted({atoms[i - 1].residue_number for i in group.atom_ids})
                entry['residue_numbers'] = residues
                entry['provenance'].append(f'{rec.structure_path}: atom positions (not wrapped atom serials)')
                if group.name in REFERENCES:
                    entry['legacy_reference_match'] = residues == REFERENCES[group.name]
                    entry['reference_note'] = 'Comparison only; GRO numbering may differ from historical numbering.'
        mapping[role] = entry
    required = ['protein', 'backbone', 'ligand', 'active_site', 'catalytic_site']
    missing = [r for r in required if not mapping[r]['confidence']]
    profile = None
    if cofactor_names:
        profile = 'mechanistic_cofactor'
    elif compound['value'] and mapping['protein']['value'] and mapping['ligand']['value']:
        ns = duration['value']
        if ns is not None and ns >= 100:
            profile = 'mechanistic'
        elif ns is not None and 0 < ns <= 50 and target['value']:
            profile = 'xanthone_short'
    if profile is None:
        issues.append('insufficient evidence for study profile')
    if compound['value'] is None or target['value'] is None:
        issues.append('unresolved compound or target')
    elif min(compound['confidence'], target['confidence']) < .5:
        issues.append('compound or target identity relies on folder hints')
    if duration['confidence'] < .5 and profile != 'mechanistic_cofactor':
        issues.append('production duration requires verification')
    if missing:
        message = 'missing required groups: ' + ', '.join(missing)
        issues.append(message)
        selection_issues.add(message)
    issues += [w.message for w in rec.warnings if w.severity in ('review_required', 'error')]
    existing, rejected = existing_results(directory)
    # Corpus consensus: every complete 25ns-DM core XVG spans 0..25000 ps.
    # This is a protocol prior for A6, not proof that this A6 directory was run
    # with that interval; HMG remains review-gated until an A6-specific result exists.
    thesis_window = {'value_ps': [0, 25000], 'confidence': .85,
                     'provenance': ['Resultados-Tesis-maestría/25ns-DM: validated A1-A5 core XVG axes'],
                     'status': 'consensus_reference'}
    if target['value'] == 'HMG-R' and compound['value'] == 'A6':
        thesis_window['status'] = 'review_required: no HMG A6 thesis result in corpus'
    planned = []
    global_issues = set(issues) - selection_issues
    catalog = MECHANISTIC if profile == 'mechanistic' else SHORT + OPTIONAL_SHORT if profile == 'xanthone_short' else []
    for analysis in catalog:
        requirements = ['protein', 'ligand'] if analysis not in ('protein_rmsd', 'protein_rmsf', 'radius_of_gyration', 'sasa', 'pca', 'fel', 'clustering', 'dccm') else ['protein']
        if analysis.startswith('active_site'):
            requirements += ['active_site']
        if analysis == 'catalytic_com_distance':
            requirements += ['catalytic_site']
        blocked = [r for r in requirements if r in missing]
        planned.append({'analysis': analysis,
                        'status': 'existing' if analysis in existing else 'review_required' if blocked or global_issues else 'planned',
                        'required': analysis in (MECHANISTIC if profile == 'mechanistic' else SHORT),
                        'required_roles': requirements, 'missing_roles': blocked,
                        'existing_files': existing.get(analysis, []),
                        'execution_implemented': False,
                        'availability': ('requires compatible energy terms; md.edr present' if (directory / 'md.edr').is_file() else 'requires compatible energy terms; md.edr absent') if analysis == 'interaction_energy' else 'optional when available' if analysis == 'mmpbsa_decomposition' else 'planned capability'})
    return {'system': directory.name, 'directory': str(directory),
            'compound': compound, 'target': target,
            'production_duration': duration,
            'available_trajectory_duration': duration if (rec.trajectory_inspection and rec.trajectory_inspection.total_duration_ps) else {'value': None, 'confidence': 0, 'provenance': []},
            'analysis_window': thesis_window,
            'analysis_window_duration_ns': {'value': 25.0, 'confidence': thesis_window['confidence'], 'provenance': thesis_window['provenance']},
            # Backwards-compatible alias; callers should prefer the explicit fields above.
            'duration_ns': duration,
            'profile': fact(profile, min(compound['confidence'], target['confidence'], duration['confidence']) if profile != 'mechanistic_cofactor' else cofactor['confidence'],
                            evidence + target['provenance'] + duration['provenance'] + cofactor['provenance']),
            'cofactor': cofactor, 'status': 'review_required' if issues or profile == 'mechanistic_cofactor' else 'resolved',
            'issues': sorted(set(issues)), 'index_candidates': [str(p) for p in candidates],
            'index': str(index) if index else None,
            'groups': [{'name': g.name, 'group_id': g.group_id, 'atom_count': len(g.atom_ids)} for g in groups],
            'selections': mapping, 'required_groups_found': not missing,
            'existing_results': existing, 'unvalidated_results': rejected, 'analyses': planned,
            'required_analyses': THESIS_CORE if profile == 'xanthone_short' else catalog,
            'optional_analyses': OPTIONAL_SHORT if profile == 'xanthone_short' else [],
            'note': 'requires dedicated analysis mapping' if profile == 'mechanistic_cofactor' else 'Dry-run only; no analysis executed.'}
