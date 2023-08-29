from clease.datastructures import MCStep
from clease.montecarlo.observers import MCObserver
import ase
from clease.datastructures import MCStep
import numpy as np
import matplotlib.pyplot as plt 

def get_atom_positions(atom_positions,atomic_number,atoms):

    '''
        Gets the positions of the atoms whose RDF is to be calculated
    '''
    
    for i in range(len(atoms)):
        if (atoms.numbers[i] in atomic_number):#Check if the element is desired in the RDF to be produced
            atom_positions.append(i)
    
    return atom_positions


def get_data(atom_positions,frequency,bins,atoms):

    '''
        Gets the No. of atoms between a particular range of distance from the reference atom
    '''

    dist=atoms.get_distances(0,atom_positions,vector=False)

    for d in dist:
        for i,bin in enumerate(bins):
            if(d<=bin):
                frequency[i]+=1
                break

    return(frequency)


def save_plot_graph(frequency,bins):

    '''
        Plots RDF from the reference atom
    '''
   
    plt.plot(bins,frequency, color='blue', linewidth=1)
    plt.xlabel('Distance (in Å)',size=15)
    plt.ylabel('No.of Atoms ',size=15)
    plt.title("Radial Distribution Function")
    plt.savefig("Radial Distribution Function")


def get_radial(atoms,number_of_bins=100,max_distance=6.50):

    elements=[29,79]#Atomic Number of Elements to be included in the RDF
    bins=np.linspace(0,max_distance,number_of_bins)
    bins=np.round(bins,2)
    frequency=[0]*number_of_bins
    atom_positions=[]
    atom_positions=get_atom_positions(atom_positions=atom_positions,atomic_number=elements,atoms=atoms)
    frequency=get_data(atom_positions=atom_positions,frequency=frequency,bins=bins,atoms=atoms)
    save_plot_graph(frequency=frequency,bins=bins)
    
           

class RDFObserver(MCObserver):


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
        state if it has a lower energy than the previous state then get the 
        RDF for the Atoms Object .

        """
        if mc_step.energy < self.lowest_energy:

            self.lowest_energy = mc_step.energy
            self.emin_atoms.numbers = self.atoms.numbers
            get_radial(self.atoms)