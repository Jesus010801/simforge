"""Legacy plans are metadata only, with per-file GROMACS group positions."""
import hashlib
import json
from pathlib import Path

import pytest

from analysis.campaign.legacy import annotate_legacy, existing_results, REFERENCES
from analysis.campaign.models import SystemRecord, TrajectoryArtifact
from analysis.campaign.structure.index_groups import resolve_index_group
from analysis.campaign.manifest import build_manifest
from analysis.campaign.orchestration.study_analyzer import run_inspect


def system(tmp_path, name='AA-A6', ligand_id=13, site_id=21, duration=25000):
    d = tmp_path / name
    d.mkdir()
    names = ['System', 'Protein', 'Protein-H', 'C-alpha', 'Backbone']
    names += [f'Group{i}' for i in range(5, 25)]
    names[ligand_id] = 'LIG'
    suffix = 'HMG' if 'HMG' in name else 'AA'
    names[site_id] = f'ActiveSite_{suffix}'
    names[24] = f'Catalytic_{suffix}'
    (d / 'index.ndx').write_text(''.join(f'[ {n} ]\n1 2\n' for n in names))
    (d / 'A6.itp').write_text('[ moleculetype ]\nLIG 3\n')
    (d / 'md.xtc').write_bytes(b'fake trajectory')
    (d / 'topol.top').write_text('#include "A6.itp"\n[ molecules ]\nProtein 1\nLIG 1\n')
    (d / 'md.mdp').write_text(f'integrator = md\ndt = 0.002\nnsteps = {duration / .002:.0f}\n')
    r = SystemRecord('test', 'test', 'rep01', index_path=str(d/'index.ndx'),
                     production_trajectory_paths=[str(d/'md.xtc')],
                     trajectory_artifacts=[TrajectoryArtifact(path=str(d/'md.xtc'), stage='production', end_time_ps=duration)])
    return d, r


def snapshot(root):
    return {str(p.relative_to(root)): (p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).hexdigest())
            for p in root.rglob('*') if p.is_file()}


def test_independent_group_positions(tmp_path):
    a, ra = system(tmp_path)
    b, rb = system(tmp_path, 'HMG-R-25ns-A6', ligand_id=14, site_id=23)
    assert resolve_index_group(a/'index.ndx', 'LIG').group_id == 13
    assert resolve_index_group(b/'index.ndx', 'LIG').group_id == 14
    assert resolve_index_group(b/'index.ndx', 'ActiveSite_HMG').group_id == 23
    assert annotate_legacy(rb, b)['selections']['active_site']['group_id'] == 23


def test_missing_catalytic(tmp_path):
    d, r = system(tmp_path)
    p = d/'index.ndx'
    p.write_text(p.read_text().replace('Catalytic_AA', 'Unknown'))
    result = annotate_legacy(r, d)
    assert not result['required_groups_found']
    item = next(a for a in result['analyses'] if a['analysis'] == 'catalytic_com_distance')
    assert item['status'] == 'review_required'
    assert item['missing_roles'] == ['catalytic_site']


def test_multiple_indices_and_backups(tmp_path):
    d, r = system(tmp_path)
    (d/'#index.ndx.1#').write_text('[ LIG ]\n1\n')
    assert annotate_legacy(r, d)['index'] == str(d/'index.ndx')
    (d/'other.ndx').write_text('[ LIG ]\n1\n')
    result = annotate_legacy(r, d)
    assert result['index'] is None
    assert len(result['index_candidates']) == 2
    assert result['status'] == 'review_required'
    (d/'index.ndx').unlink()
    (d/'other.ndx').unlink()
    assert annotate_legacy(r, d)['index'] is None


def test_profiles_and_name_conflict(tmp_path):
    d, r = system(tmp_path)
    assert annotate_legacy(r, d)['profile']['value'] == 'xanthone_short'
    d, r = system(tmp_path, 'HMG-R-200ns-A6', duration=200000)
    result = annotate_legacy(r, d)
    assert result['profile']['value'] == 'mechanistic'
    assert {'pca', 'fel', 'clustering', 'dccm'} <= {a['analysis'] for a in result['analyses']}
    d, r = system(tmp_path, 'HMG-R-25ns-A6', duration=200000)
    assert 'folder duration conflicts with production metadata' in annotate_legacy(r, d)['issues']


def test_cofactor_and_atom_names(tmp_path):
    d, r = system(tmp_path, 'system_A6_COA')
    with (d/'index.ndx').open('a') as f:
        f.write('[ COA ]\n1\n')
    (d/'protein.itp').write_text('[ atoms ]\n1 C 1 ALA A1 1 0 12\n')
    result = annotate_legacy(r, d)
    assert result['profile']['value'] == 'mechanistic_cofactor'
    assert result['compound']['value'] == 'A6'
    assert result['cofactor']['value'] == ['COA']
    assert result['analyses'] == []
    assert result['status'] == 'review_required'


def test_existing_results_require_content(tmp_path):
    (tmp_path/'binding_energy.xvg').write_text('0 1\n1 2\n')
    assert 'mmpbsa' not in existing_results(tmp_path)[0]
    (tmp_path/'energy_MM.xvg').write_text('@ title "Energy"\n0 -5\n1 -6\n')
    assert 'mmpbsa' not in existing_results(tmp_path)[0]
    (tmp_path/'energy_summary.csv').write_text('Energy,Average,Standard-Deviation,\nvDW,-100,2,\nPolar-solvation,25,1,\nNon-polar-solvation,-2,1,\nTotal,-77,2,\n')
    (tmp_path/'residues_energy_summary.csv').write_text('Residue,polar,apolar,total,\nALA-1,1,2,3,\n')
    found, rejected = existing_results(tmp_path)
    assert {'mmpbsa', 'mmpbsa_decomposition'} <= set(found)
    assert str(tmp_path/'binding_energy.xvg') in rejected
    (tmp_path/'binding_energy.xvg').write_text('@ title "Binding energy"\n0 -5\n1 -6\n')
    assert str(tmp_path/'binding_energy.xvg') in existing_results(tmp_path)[0]['mmpbsa']


def test_invalid_numeric_results(tmp_path):
    for content in ['@ title "RMSD"\n0 nan\n1 2\n', '@ title "RMSD"\n0 1\n', '@ title "RMSD"\n0 1\n1 2 3\n']:
        (tmp_path/'rmsd_protein.xvg').write_text(content)
        assert 'protein_rmsd' not in existing_results(tmp_path)[0]


def test_deterministic_read_only_manifest_and_inspect(tmp_path):
    d, r = system(tmp_path)
    before = snapshot(tmp_path)
    a = build_manifest(tmp_path, inspect_trajectories=False)
    b = build_manifest(tmp_path, inspect_trajectories=False)
    assert json.dumps(a.systems[0].legacy_study, sort_keys=True) == json.dumps(b.systems[0].legacy_study, sort_keys=True)
    assert SystemRecord.from_dict(a.systems[0].to_dict()).legacy_study == a.systems[0].legacy_study
    result = run_inspect(tmp_path, inspect_trajectories=False)
    assert result.output_files == []
    assert before == snapshot(tmp_path)
    assert not (tmp_path/'simforge_analysis').exists()


def test_exact_alias_and_duplicate_resolution(tmp_path):
    p = tmp_path/'index.ndx'
    p.write_text('[ LIG ]\n1\n[ A6 ]\n2\n')
    assert resolve_index_group(p, 'LIG', ('A6',)).group_id == 0
    assert resolve_index_group(p, 'ligand', ('A6',)).group_id == 1
    assert resolve_index_group(p, 'lig') is None
    with p.open('a') as f:
        f.write('[ LIG ]\n3\n')
    with pytest.raises(ValueError, match='ambiguous'):
        resolve_index_group(p, 'LIG')


def test_residue_recovery_uses_atom_position(tmp_path):
    d, r = system(tmp_path)
    p = d/'md.gro'
    residues = REFERENCES['Catalytic_AA']
    p.write_text('test\n3\n'+''.join(f'{res:5d}{"ALA":<5}{"CA":>5}{99999-i:5d}{0:8.3f}{0:8.3f}{0:8.3f}\n' for i,res in enumerate(residues))+'1 1 1\n')
    r.structure_path = str(p)
    index = d/'index.ndx'
    index.write_text(index.read_text().replace('[ Catalytic_AA ]\n1 2', '[ Catalytic_AA ]\n1 2 3'))
    selection = annotate_legacy(r, d)['selections']['catalytic_site']
    assert selection['residue_numbers'] == residues
    assert selection['legacy_reference_match'] is True


def test_existing_results_are_not_planned(tmp_path):
    d, r = system(tmp_path)
    (d/'binding_energy.xvg').write_text('@ title "Binding energy"\n0 -5\n1 -6\n')
    result = annotate_legacy(r, d)
    assert next(a for a in result['analyses'] if a['analysis'] == 'mmpbsa')['status'] == 'existing'


def test_empty_and_out_of_range_groups(tmp_path):
    d, r = system(tmp_path)
    p = d/'index.ndx'
    p.write_text(p.read_text().replace('[ LIG ]\n1 2', '[ LIG ]'))
    assert not annotate_legacy(r, d)['required_groups_found']


def test_export_does_not_overwrite(tmp_path):
    system(tmp_path)
    output = tmp_path/'export'
    run_inspect(tmp_path, inspect_trajectories=False, output_dir=output)
    before = snapshot(output)
    with pytest.raises(FileExistsError):
        run_inspect(tmp_path, inspect_trajectories=False, output_dir=output)
    assert before == snapshot(output)


def test_shared_discovery_selects_xtc_and_excludes_mdfit(tmp_path):
    d, r = system(tmp_path)
    (d/'md.trr').write_bytes(b'full precision')
    (d/'mdfit.xtc').write_bytes(b'derived')
    rec = build_manifest(tmp_path, inspect_trajectories=False).systems[0]
    assert rec.production_trajectory_paths == [str(d/'md.xtc')]
    assert rec.artifact(str(d/'mdfit.xtc')).stage == 'derived'


def test_selections_and_plan_cli_are_read_only(tmp_path):
    from typer.testing import CliRunner
    import typer
    from analysis.campaign.cli import study_plan_fn, study_selections_fn
    system(tmp_path)
    before = snapshot(tmp_path)
    app = typer.Typer()
    app.command('plan')(study_plan_fn)
    app.command('selections')(study_selections_fn)
    runner = CliRunner()
    outputs = []
    for command in ['plan', 'selections', 'plan']:
        result = runner.invoke(app, [command, str(tmp_path), '--json'])
        assert result.exit_code == 0, result.output
        outputs.append(json.loads(result.output))
    assert outputs[0] == outputs[1] == outputs[2]
    assert before == snapshot(tmp_path)



def test_unused_parameter_templates_do_not_change_identity(tmp_path):
    d, r = system(tmp_path)
    (d / 'COA.itp').write_text('[ moleculetype ]\nCOA 3\n')
    (d / 'A3.itp').write_text('[ moleculetype ]\nA3 3\n')
    result = annotate_legacy(r, d)
    assert result['cofactor']['value'] == []
    assert result['compound']['value'] == 'A6'
    assert result['profile']['value'] == 'xanthone_short'
    assert result['compound']['provenance']


def test_zero_count_cofactor_not_present(tmp_path):
    d, r = system(tmp_path)
    with (d / 'topol.top').open('a') as stream:
        stream.write('COA 0\n')
    assert annotate_legacy(r, d)['cofactor']['value'] == []
    with (d / 'topol.top').open('a') as stream:
        stream.write('COA 1\n')
    assert annotate_legacy(r, d)['cofactor']['value'] == ['COA']


def test_standalone_ligand_template_requires_index_corroboration(tmp_path):
    d, r = system(tmp_path)
    (d / 'topol.top').unlink()
    assert annotate_legacy(r, d)['compound']['confidence'] == .9
    p = d / 'index.ndx'
    p.write_text(p.read_text().replace('[ LIG ]', '[ Unknown ]'))
    assert annotate_legacy(r, d)['compound']['confidence'] < .5


def test_missing_site_does_not_block_unrelated_analysis(tmp_path):
    d, r = system(tmp_path)
    p = d / 'index.ndx'
    p.write_text(p.read_text().replace('Catalytic_AA', 'Unknown'))
    result = annotate_legacy(r, d)
    analyses = {a['analysis']: a for a in result['analyses']}
    assert result['status'] == 'review_required'
    assert analyses['catalytic_com_distance']['status'] == 'review_required'
    assert analyses['protein_rmsd']['status'] == 'planned'
    assert analyses['active_site_contacts']['status'] == 'planned'
