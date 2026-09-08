"""Conservative format validation for legacy result discovery."""
import pytest
from analysis.campaign.legacy import existing_results


@pytest.mark.parametrize('filename,title,family', [
    ('rmsd_protein.xvg', 'RMSD', 'protein_rmsd'),
    ('rmsd_ligand.xvg', 'RMSD', 'ligand_rmsd'),
    ('rmsf_protein.xvg', 'RMSF', 'protein_rmsf'),
    ('sasa_protein.xvg', 'Solvent accessible area', 'sasa'),
    ('rg_protein.xvg', 'Radius of gyration', 'radius_of_gyration'),
    ('hbonds_lig.xvg', 'Hydrogen bonds', 'hbonds'),
    ('mindist_lig_active.xvg', 'Minimum distance', 'active_site_mindist'),
    ('contacts_lig_active.xvg', 'Number of contacts', 'active_site_contacts'),
    ('dist_lig_catalytic.xvg', 'COM distance', 'catalytic_com_distance'),
    ('interaction_energy.xvg', 'Coul-SR:Protein-LIG LJ-SR:Protein-LIG', 'interaction_energy'),
    ('clusters_size.xvg', 'Cluster size distribution', 'clustering'),
])
def test_valid_numeric_families(tmp_path, filename, title, family):
    path = tmp_path / filename
    path.write_text(f'@ title "{title}"\n0 1\n1 2\n')
    assert existing_results(tmp_path)[0][family] == [str(path)]
    path.write_text(f'# generated command {title}\n@ title "Unknown"\n0 1\n1 2\n')
    assert family not in existing_results(tmp_path)[0]


@pytest.mark.parametrize('body', ['0 1\n', '0 1\n1 nan\n', '0 1\n1 2 3\n', ''])
def test_partial_xvg(tmp_path, body):
    (tmp_path / 'rmsd_protein.xvg').write_text('@ title "RMSD"\n' + body)
    assert 'protein_rmsd' not in existing_results(tmp_path)[0]


def test_pca_needs_both_components_and_deduplicates(tmp_path):
    (tmp_path / 'eigenval.xvg').write_text('@ title "Eigenvalues"\n1 2\n2 1\n')
    assert 'pca' not in existing_results(tmp_path)[0]
    (tmp_path / 'proj.xvg').write_text('@ title "Projection"\n0 1\n1 2\n')
    assert len(existing_results(tmp_path)[0]['pca']) == 2


@pytest.mark.parametrize('filename,title,family', [
    ('fel.xpm', 'Gibbs free energy', 'fel'),
    ('dccm.xpm', 'Dynamic cross-correlation', 'dccm'),
    ('clusters.xpm', 'Clusters', 'clustering'),
])
def test_matrix_requires_complete_payload(tmp_path, filename, title, family):
    path = tmp_path / filename
    header = f'/* XPM */\n/* title: "{title}" */\nstatic char *matrix[] = {{\n"2 2 2 1",\n'
    path.write_text(header)
    assert family not in existing_results(tmp_path)[0]
    path.write_text(header + '"a c #ffffff",\n"b c #000000",\n"ab",\n"ba"\n};\n')
    assert family in existing_results(tmp_path)[0]
    path.write_text(path.read_text().replace('"ba"', '"bz"'))
    assert family not in existing_results(tmp_path)[0]


def test_csv_column_order_and_missing_values(tmp_path):
    path = tmp_path / 'energy_summary.csv'
    path.write_text('Average,Energy\n3,Polar-solvation\n-1,Non-polar-solvation\n2,Total\n')
    assert 'mmpbsa' in existing_results(tmp_path)[0]
    path.write_text('Average,Energy\n3,Polar-solvation\n-1,Non-polar-solvation\n,Total\n')
    assert 'mmpbsa' not in existing_results(tmp_path)[0]


def test_different_cofactor_pairs_not_combined(tmp_path):
    pair = tmp_path / 'COA-ligand-MMPBSA'
    pair.mkdir()
    (pair / 'binding_energy.xvg').write_text('@ title "Binding energy"\n0 1\n1 2\n')
    assert 'mmpbsa' not in existing_results(tmp_path)[0]
    assert 'mmpbsa' in existing_results(pair)[0]


def test_deterministic_and_read_only_with_spaces(tmp_path):
    directory = tmp_path / 'system with spaces'
    directory.mkdir()
    path = directory / 'binding_energy.xvg'
    path.write_text('@ title "Binding energy"\n0 1\n1 2\n')
    before = {p.name: p.read_bytes() for p in directory.iterdir()}
    assert existing_results(directory) == existing_results(directory)
    assert {p.name: p.read_bytes() for p in directory.iterdir()} == before


def test_interaction_energy_requires_both_components(tmp_path):
    path = tmp_path / 'interaction_energy.xvg'
    path.write_text('@ title "Coul-SR:Protein-LIG"\n0 1\n1 2\n')
    assert 'interaction_energy' not in existing_results(tmp_path)[0]


@pytest.mark.parametrize('filename,title', [
    ('hbonds_lig.xvg', 'Hydrogen bond autocorrelation'),
    ('rmsd_protein.xvg', 'RMSD distribution'),
    ('rmsf_ligand.xvg', 'RMS fluctuation of ligand'),
])
def test_related_scientific_family_not_completed(tmp_path, filename, title):
    (tmp_path / filename).write_text(f'@ title "{title}"\n0 1\n1 2\n')
    assert not existing_results(tmp_path)[0]


def test_realistic_gromacs_rmsf_title(tmp_path):
    (tmp_path / 'rmsf.xvg').write_text('@ title "RMS fluctuation"\n1 0.1\n2 0.2\n')
    assert 'protein_rmsf' in existing_results(tmp_path)[0]


def test_pca_conflicting_explicit_windows_and_alternatives(tmp_path):
    (tmp_path / 'eigenval.xvg').write_text('# gmx covar -b 0 -e 100\n@ title "Eigenvalues"\n1 2\n2 1\n')
    projection = tmp_path / 'proj.xvg'
    projection.write_text('# gmx anaeig -b 50 -e 100\n@ title "Projection"\n0 1\n1 2\n')
    assert 'pca' not in existing_results(tmp_path)[0]
    projection.write_text(projection.read_text().replace('-b 50', '-b 0'))
    assert 'pca' in existing_results(tmp_path)[0]
    (tmp_path / 'proj_other.xvg').write_text(projection.read_text())
    assert 'pca' not in existing_results(tmp_path)[0]


def test_duplicate_csv_summary_rows_are_ambiguous(tmp_path):
    (tmp_path / 'energy_summary.csv').write_text('Energy,Average\nPolar-solvation,1\nNon-polar-solvation,2\nTotal,3\nTotal,4\n')
    assert 'mmpbsa' not in existing_results(tmp_path)[0]


def test_xpm_missing_closing_array_is_partial(tmp_path):
    (tmp_path / 'fel.xpm').write_text('/* XPM */\n/* title: "Free energy" */\n"1 1 1 1",\n"a c #ffffff",\n"a"\n')
    assert 'fel' not in existing_results(tmp_path)[0]
