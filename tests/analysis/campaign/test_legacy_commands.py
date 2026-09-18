import json
from pathlib import Path
from analysis.campaign.legacy_commands import build_xanthone_commands
from analysis.campaign.manifest import build_manifest

def test_dynamic_ids_and_no_execution(tmp_path):
    d=tmp_path / "A A-6"; d.mkdir()
    (d/'index.ndx').write_text("[ System ]\n1\n[ Protein ]\n1\n[ Backbone ]\n1\n[ LIG ]\n1\n[ ActiveSite_AA ]\n1\n[ Catalytic_AA ]\n1\n")
    (d/'md.xtc').write_bytes(b'x'); (d/'md.tpr').write_bytes(b'x')
    legacy={'system':'AA-A6','directory':str(d),'profile':{'value':'xanthone_short'},'status':'resolved','index':str(d/'index.ndx'),'selections':{r:{'group_id':i} for i,r in enumerate(['protein','backbone','ligand','active_site','catalytic_site'])},'existing_results':{}}
    plan=build_xanthone_commands(None,legacy,gmx='gmx-test')
    assert plan['execution']=='disabled; dry-run only'
    assert plan['commands'][1]['argv'][0]=='gmx-test'
    assert {r['analysis'] for r in plan['commands']} >= {'interaction_energy','mmpbsa'}
    assert plan['commands'][1]['group_ids']['protein']==0

def test_review_hmg_never_planned(tmp_path):
    d=tmp_path/'HMG'; d.mkdir()
    legacy={'system':'HMG-R-25ns-A6','directory':str(d),'profile':{'value':'mechanistic'},'status':'review_required','selections':{},'existing_results':{}}
    assert {x['status'] for x in build_xanthone_commands(None,legacy)['commands']} == {'REVIEW_REQUIRED'}

def test_analysis_rows_require_mdfit_and_never_fall_back_to_md_xtc(tmp_path):
    d = tmp_path / "AA-A6"; d.mkdir()
    (d / 'index.ndx').write_text("[ System ]\n1\n[ Protein ]\n1\n[ LIG ]\n1\n[ ActiveSite_AA ]\n1\n[ Catalytic_AA ]\n1\n")
    (d / 'md.xtc').write_bytes(b'x'); (d / 'md.tpr').write_bytes(b'x')
    legacy = {'system': 'AA-A6', 'directory': str(d), 'profile': {'value': 'xanthone_short'},
              'status': 'resolved', 'index': str(d / 'index.ndx'),
              'selections': {r: {'group_id': i} for i, r in enumerate(
                  ['protein', 'ligand', 'active_site', 'catalytic_site'])},
              'existing_results': {}}

    # no mdfit.xtc -> every analysis row is REVIEW_REQUIRED and names mdfit.xtc
    plan = build_xanthone_commands(None, legacy)
    analysis_rows = [r for r in plan['commands']
                     if r['analysis'] in ('protein_rmsd', 'ligand_rmsd',
                                          'active_site_mindist', 'catalytic_com_distance')]
    assert {r['status'] for r in analysis_rows} == {'REVIEW_REQUIRED'}
    assert all('mdfit.xtc' in r['reason'] for r in analysis_rows)
    for r in plan['commands']:
        f = r['argv'][r['argv'].index('-f') + 1] if '-f' in r['argv'] else ''
        assert not f.endswith('/md.xtc')

    # add mdfit.xtc -> the same rows now PLAN, bound to mdfit.xtc
    (d / 'mdfit.xtc').write_bytes(b'fit')
    plan = build_xanthone_commands(None, legacy)
    for r in plan['commands']:
        if r['analysis'] in ('protein_rmsd', 'ligand_rmsd', 'active_site_mindist'):
            assert r['status'] == 'PLANNED'
            assert r['argv'][r['argv'].index('-f') + 1].endswith('/mdfit.xtc')
    # the PBC/fit preprocessing still consumes the raw production trajectory
    center = plan['preprocessing'][0]['argv']
    assert center[center.index('-f') + 1].endswith('md.xtc')


def test_preprocessing_and_window_are_explicit(tmp_path):
    d=tmp_path/'AA-A6'; d.mkdir(); (d/'md.xtc').write_bytes(b'x'); (d/'md.tpr').write_bytes(b'x')
    legacy={'system':'AA-A6','directory':str(d),'profile':{'value':'xanthone_short'},'status':'resolved','index':str(d/'index.ndx'),'selections':{'protein':{'group_id':1},'ligand':{'group_id':13},'active_site':{'group_id':21},'catalytic_site':{'group_id':22}},'existing_results':{},'analysis_window':{'value_ps':[0,25000],'status':'consensus_reference'}}
    plan=build_xanthone_commands(None,legacy)
    assert [x['step'] for x in plan['preprocessing']] == ['center','fit']
    assert plan['preprocessing'][0]['argv'][5:12] == ['-pbc','res','-ur','compact','-center','-o','/tmp/never'] or '-pbc' in plan['preprocessing'][0]['argv']
    assert plan['analysis_window']['value_ps'] == [0,25000]
