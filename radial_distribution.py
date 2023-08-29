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
from ase.ga.utilities import get_atoms_distribution,get_mic_distance
from clease.settings import Concentration,CEBulk
from clease.calculator import attach_calculator
c=0
def alternate(atoms):
    cell=atoms
    i=1
    Cu_no=0
    while(i <len(cell)):
        cell.symbols[i]='Cu'
        i+=2
        Cu_no+=1
    return cell,Cu_no
def half(atoms):
    cell=atoms
    i=0
    while(i<(len(cell)/2)):
        cell.symbols[i]='Cu'
        i+=1    
    return cell
size=(3,3,3)
tot=np.prod(size)
atoms = connect('aucu.db').get(id=1).toatoms()*size
alt,Cu_no=alternate(atoms)
max_distance=15
number_of_bins=100
elements=[]
bins=np.linspace(0,max_distance,number_of_bins)
bins = np.round(bins,2)
frequency=[0 for i in range(number_of_bins)]
atom_positions=[]
view(alt)
def get_atom_positions(atom_positions,atomic_number,atoms):
    for i in range(len(atoms)):
        if (atoms.numbers[i] in atomic_number):
            atom_positions.append(i)
    return atom_positions
atom_positions=get_atom_positions(atom_positions=atom_positions,atomic_number=[29,79],atoms=alt)

def get_data(atom_positions,bins,atoms):
    dist=atoms.get_distances(0,atom_positions,vector=False)
    for d in dist:
        for i,bin in enumerate(bins):
            if(d<=bin):
                frequency[i]+=1
                break
    return(frequency)
frequency=get_data(atom_positions=atom_positions,bins=bins,atoms=alt)
def plot_graph(frequency,bins):
    plt.plot(bins,frequency, color='blue', linewidth=1)
    plt.xlabel('Distance (in Å)',size=15)
    plt.ylabel('No.of Atoms ',size=15)
    plt.title("Radial Distribution Function")
    plt.show()
plot_graph(frequency=frequency,bins=bins)





'''
#view(alt)
#view(half(atoms))
print(alt.get_chemical_formula())
max_distance=15
number_of_bins=100
no_count_types=[29]
atoms_distribution=get_atoms_distribution(alt,center=(0,0,0),number_of_bins=number_of_bins,max_distance=max_distance,no_count_types=[29])
#print(atoms_distribution)
atom_positions=[]

for i in range(1,len(alt)):
    if(atoms.symbols[i]=='Au'):
        Au_positions.append(i)
distances=[]
for i, a in enumerate(atoms_distribution):
        dist=i*(max_distance/number_of_bins)
        distances.append(dist)
        if(a !=0):
            print(dist,"-->",a)
print(sum(atoms_distribution))
#for i in range (len(Au_positions)):
#    print("Distance from Atom 0 to Atom Number {} is {}".format(Au_positions[i],alt.get_distances(0,Au_positions[i],vector=False)))
bins=[0]*number_of_bins
distances=np.linspace(0,max_distance,number_of_bins)
d=[]

for i in range(1,len(atoms)):
    distance=alt.get_distances(0,i,vector=False)
    if(distance>max_distance):
        continue
    d.append(distance)
for dist in d:
    for i in range(len(bins)):
        if(dist<=distances[i]):
            bins[i]+=1
            break
print(bins,"-->",distances)
import matplotlib.pyplot as plt
plt.plot(distances,bins, color='blue', linewidth=1)
plt.xlabel('r',size=15)
plt.ylabel('g(r)',size=15)
plt.title('Radial Distribution Function between Au and Cu ')
plt.show()
'''
            