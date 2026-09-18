import numpy as np
from validators.membrane_water_v2.identity import load_system,PROTEIN,LIPID,OTHER
from validators.membrane_water_v2.cell import PeriodicCell,Lattice
from validators.membrane_water_v2.clearance import ClearanceIndex
from validators.membrane_water_v2.models import SystemGeometry,ExclusionModel


def _make(path,waters=True,offset=0,atomshift=0,model="SOL"):
    atoms=[(100+offset,"ALA","CA",1.,1.,2.),(101+offset,"DPP","P",2.,2.,4.),(101+offset,"DPP","C50",2.,2.,3.5)]
    if waters:
        atoms += [(7+offset,model,"OW",3.,3.,3.),(7+offset,model,"HW1",3.01,3.,3.),
                  (7+offset,model,"HW2",3.,3.01,3.)]
        if model in {"TIP4","TIP5"}: atoms += [(7+offset,model,"MW1",3.,3.,3.)]
        if model=="TIP5": atoms += [(7+offset,model,"MW2",3.,3.,3.)]
    lines=["identity",str(len(atoms))]
    for i,(r,rn,an,x,y,z) in enumerate(atoms):
        lines.append(f"{r%100000:5d}{rn:<5}{an:>5}{(i+1+atomshift)%100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    lines.append("6.0 6.0 8.0");path.write_text("\n".join(lines)+"\n");return path


def test_geometry_identical_with_and_without_water(tmp_path):
    a=load_system(_make(tmp_path/"wet.gro",True)).geometry
    b=load_system(_make(tmp_path/"dry.gro",False)).geometry
    np.testing.assert_array_equal(a.coordinates,b.coordinates)
    np.testing.assert_array_equal(a.species,b.species)
    np.testing.assert_array_equal(a.radii,b.radii)


def test_residue_and_atom_numbering_do_not_define_identity(tmp_path):
    a=load_system(_make(tmp_path/"a.gro",True,offset=0,atomshift=0))
    b=load_system(_make(tmp_path/"b.gro",True,offset=99999,atomshift=99998))
    np.testing.assert_array_equal(a.geometry.coordinates,b.geometry.coordinates)
    assert [w.uid for w in a.waters]==[w.uid for w in b.waters]==[0]
    assert a.waters[0].residue_number!=b.waters[0].residue_number


def test_water_models_are_post_geometry_species(tmp_path):
    signatures=[]
    for i,model in enumerate(("SOL","TIP4","TIP5")):
        s=load_system(_make(tmp_path/f"w{i}.gro",True,model=model))
        signatures.append((s.geometry.coordinates.copy(),s.geometry.species.copy(),s.geometry.radii.copy()))
        assert s.waters[0].valid
    for a,b in zip(signatures[0],signatures[1]): np.testing.assert_array_equal(a,b)
    for a,b in zip(signatures[0],signatures[2]): np.testing.assert_array_equal(a,b)


def test_rotated_cell_clearance_is_equivariant():
    base=np.diag([6.,7.,8.]);theta=.6
    rot=np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta),0,np.cos(theta)]])
    p=np.array([[1.,2.,3.]])
    def geom(h,pt):
        c=PeriodicCell(h)
        return SystemGeometry(c,pt,np.array([.17]),np.array([1],np.uint8),("C",),np.array([0]),np.array([True]))
    model=ExclusionModel()
    ca,ha,_=ClearanceIndex(geom(base,p)).query(p,model)
    pr=p@rot.T;cr,hr,_=ClearanceIndex(geom(rot@base,pr)).query(pr,model)
    np.testing.assert_allclose(ca,cr,atol=1e-10)
    np.testing.assert_allclose(ha,hr,atol=1e-10)


def test_cell_period_matches_exact_box_under_spacing_and_phase():
    cell=PeriodicCell(np.diag([12.44332,12.44140,19.2]))
    for phase in ((.5,.5,.5),(.13,.71,.29)):
        grid=Lattice.create(cell,.17,phase)
        np.testing.assert_allclose(grid.steps*np.asarray(grid.shape),cell.lengths,rtol=0,atol=1e-12)
        p=np.array([[.001,4.,10.]])
        np.testing.assert_array_equal(grid.indices(p),grid.indices(p+cell.matrix[:,0]))

