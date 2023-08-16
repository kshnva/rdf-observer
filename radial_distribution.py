import numpy as np
from ase.visualize import view
from ase.geometry.analysis import Analysis
from ase.lattice.cubic import FaceCenteredCubic
from ase.build import bulk
from functools import lru_cache
import matplotlib.pyplot as plt
from ase.build import bulk
from ase.ga.data import PrepareDB
from ase.db import connect
from itertools import combinations as com
from ase.ga.utilities import get_rdf
import random
from scipy.interpolate import make_interp_spline
size=(10,10,10)
tot=np.prod(size)
atoms = connect('aucu.db').get(id=1).toatoms()*size
random_numbers = [random.randint(0, tot-1) for _ in range(int(tot/3))]
for ele in random_numbers:
	atoms[ele].symbol='Cu'
ana=Analysis(atoms)
print(atoms.get_chemical_formula(empirical=False))
rad=ana.get_rdf(rmax=10.00,nbins=50,elements=('Au','Cu'), return_dists=True)
x=rad[0][1]
y=rad[0][0]
plt.plot(x,y, color='blue', linewidth=1)
plt.xlabel('r',size=15)
plt.ylabel('g(r)',size=15)
plt.title('Radial Distribution Function between Au and Cu ')
plt.show()