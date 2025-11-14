import bobkit.util as util
import bobkit.clipper as clipper
import pickle
import faulthandler
faulthandler.enable()

s = util.read_structure("/Users/swh514/Work/data/emd_3488/pdb5ni1.pdb")
a = s.atom_list()
d = pickle.dumps(a)
l = pickle.loads(d)

