import bobkit.util as util
import bobkit.clipper as clipper
import pickle
import faulthandler
faulthandler.enable()


def find_unpicklable_object(obj):
    try:
        pickle.dumps(obj)
        return None  # Object is picklable
    except Exception as e:
        # If it's a simple error, report the object immediately
        if "bad_cast" in str(e) or not hasattr(obj, '__tuple__'):
            print("bad cast")
            return obj 
        
        # Otherwise, check the attributes recursively
        for key, value in obj.__tuple__.items():
            culprit = find_unpicklable_object(value)
            if culprit:
                print(f"Culprit found in attribute '{key}': {type(culprit)}")
                return culprit
        return obj # The top-level object is still the issue

s = util.read_structure("/Users/swh514/Work/data/emd_3488/pdb5ni1.pdb")
a = s.atom_list()
d = pickle.dumps(a)
l = pickle.loads(d)

cell = s.cell
m = pickle.dumps(cell)
l = pickle.loads(m)

g = clipper.Grid_sampling(10, 11, 12)
gm = pickle.dumps(g)
print(gm)
gl = pickle.loads(gm)
print(gl)
#res = s.model[0][0]
#m = pickle.dumps(res)
#l = pickle.loads(m)
#print(l.id)
#print(l.pos)
#m = pickle.dumps(s)
#l = pickle.loads(m)
