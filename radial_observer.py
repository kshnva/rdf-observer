from clease.datastructures import MCStep, SystemChanges
from clease.montecarlo.observers import MCObserver
from ase.geometry import get_layers
from ase.db import connect
from clease.tools import update_db
from clease.calculator import attach_calculator
from clease.montecarlo import Montecarlo
from clease.settings import Concentration,CEBulk
import ase
from clease.datastructures import MCStep
import json 
import numpy as np
import matplotlib.pyplot as plt 
from ase.geometry.analysis import Analysis
def get_radial(atoms):
    ana=Analysis(atoms)
    rad=ana.get_rdf(rmax=6.50,nbins=100,elements=None, return_dists=True)
    x=rad[0][1]
    y=rad[0][0]
    print(y)    
    plt.plot(x,y, color='blue', linewidth=1)
    plt.xlabel('r / Å',size=15)
    plt.ylabel('g(r)',size=15)
    plt.savefig(fname="Radial Distribution Function")
class LayerMonitor(MCObserver):


    def __init__(self, atoms: ase.Atoms):
        super().__init__()
        self.atoms = atoms
        self.reset()

    def reset(self) -> None:
        self.emin_atoms: ase.Atoms = self.atoms.copy()
        self.lowest_energy = np.inf


        
    def observe_step(self, mc_step: MCStep) -> None:
        """
        Check if the current state has lower energy and store the current
        state if it has a lower energy than the previous state.

        mc_step: MCStep
             Instance of MCStep with information on the latest step.
        """
        if mc_step.energy < self.lowest_energy:

            self.lowest_energy = mc_step.energy
            self.emin_atoms.numbers = self.atoms.numbers
            get_radial(self.atoms)
                
    

def load_json_file(file_path):#Load the ECI values from the JSON file
    with open(file_path, 'r') as file:
        json_data = json.load(file)
    return json_data 
file_path = 'rdf-observer/eci_l1.json'  
eci = load_json_file(file_path)
conc=Concentration(basis_elements=[['Au','Cu']])
conc.set_conc_formula_unit(formulas=["Au<x>Cu<1-x>"], variable_range={"x": (0, 1)})
settings=CEBulk(crystalstructure='fcc',
               a=3.8,
               supercell_factor=27,
               concentration=conc,
               db_name="aucu.db",
               max_cluster_dia=[6.0,5,5])

size=(6,6,6)
tot=np.product(size)
T=500
atoms = connect('aucu.db').get(id=1).toatoms()*size
atoms = attach_calculator(settings, atoms=atoms, eci=eci)
for j in range (0,int(tot/3)):
	atoms[j].symbol='Cu'
monitor = LayerMonitor(atoms)
mc = Montecarlo(atoms, T)
mc.attach(monitor, interval=50)
mc.run(steps=500)