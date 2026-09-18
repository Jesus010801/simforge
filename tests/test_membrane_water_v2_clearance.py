import numpy as np
from validators.membrane_water_v2.models import SystemGeometry, ExclusionModel
from validators.membrane_water_v2.cell import PeriodicCell
from validators.membrane_water_v2.clearance import ClearanceIndex, build_field
from validators.membrane_water_v2.void_components import components


def test_exact_periodic_unequal_spheres():
    cell = PeriodicCell(np.diag([4.,4.,4.]))
    g = SystemGeometry(cell,np.array([[.05,2,2],[1,2,2]]),np.array([.17,.30]),np.array([1,2]),("C","X"),np.array([0,1]),np.ones(2,bool))
    points = np.array([[3.98,2,2],[.6,2,2]])
    c, hard, nearest = ClearanceIndex(g).query(points,ExclusionModel())
    np.testing.assert_allclose(c,[-.10,.10],atol=1e-10)
    assert hard[0] < 0
    assert nearest.tolist()==[0,1]
    f=build_field(g,ExclusionModel(),.23)
    assert np.all(f.certified_free <= f.nominal_free)
    assert np.all(f.nominal_free <= f.possible_free)


def test_periodic_components():
    a=np.zeros((7,7,7),bool)
    a[0,3,3]=a[-1,3,3]=True
    a[3,3,3]=True
    labels,n=components(a)
    assert n==2
    assert labels[0,3,3]==labels[-1,3,3]!=labels[3,3,3]
