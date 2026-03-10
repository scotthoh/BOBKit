import bobkit.util as util
import bobkit.clipper as clipper
import pickle
import pytest
import pathlib
import faulthandler
faulthandler.enable()


#def find_unpicklable_object(obj):
#    try:
#        pickle.dumps(obj)
#        return None  # Object is picklable
#    except Exception as e:
#        # If it's a simple error, report the object immediately
#        if "bad_cast" in str(e) or not hasattr(obj, '__tuple__'):
#            print("bad cast")
#            return obj 
#        
#        # Otherwise, check the attributes recursively
#        for key, value in obj.__tuple__.items():
#            culprit = find_unpicklable_object(value)
#            if culprit:
#                print(f"Culprit found in attribute '{key}': {type(culprit)}")
#                return culprit
#        return obj # The top-level object is still the issue
@pytest.fixture
def cell_descr_instance():
    return clipper.Cell_descr(105.0, 105.0, 105.0, 90.0, 90.0, 90.0)

@pytest.fixture
def cell_instance(cell_descr_instance):
    return clipper.Cell(cell_descr_instance)

@pytest.fixture
def read_structure_instance(request, cell_instance):
    testdir = pathlib.Path(request.module.__file__).parent.parent / "test_data"
    mmol = clipper.MiniMol(clipper.Spacegroup.p1(), cell_instance)
    flag = util.read_structure(mmol, str(testdir / "pdb5ni1_cryst1.pdb"), False)  # noqa 501
    return flag, mmol

class TestSerial:
    def test_serial(self, read_structure_instance):
        flag, mmol = read_structure_instance
        # no merge_chain_parts
        assert flag
        assert len(mmol) == 12
        a = mmol.atom_list()
        d = pickle.dumps(a)
        l = pickle.loads(d)
        assert l[0].pos.x == a[0].pos.x
        assert l[0].pos.y == a[0].pos.y
        assert l[0].pos.z == a[0].pos.z

        cell = mmol.cell
        m = pickle.dumps(cell)
        l = pickle.loads(m)
        assert l.equals(cell)

        g = clipper.Grid_sampling(10, 11, 12)
        gm = pickle.dumps(g)
        gl = pickle.loads(gm)
        #assert gl.equals(g)
#res = s.model[0][0]
#m = pickle.dumps(res)
#l = pickle.loads(m)
#print(l.id)
#print(l.pos)
#m = pickle.dumps(s)
#l = pickle.loads(m)
