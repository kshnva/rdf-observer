import pandas as pd 
import numpy as np 
from clease.structgen import NewStructures
from clease.settings import Concentration,CEBulk
from ase.db import connect 
import os
from clease import Evaluate
import clease.plot_post_process as pp
import matplotlib.pyplot as plt
from ase.calculators.emt import EMT
from ase.db import connect
from clease.tools import update_db
from clease.calculator import attach_calculator
from ase.visualize import view
from clease.montecarlo import Montecarlo
from clease.montecarlo import SGCMonteCarlo
from clease.montecarlo.observers import EnergyEvolution
from clease.montecarlo.observers import Snapshot
import json 
import random
from clease.montecarlo import MCEvaluator
from ase.calculators.emt import EMT
def load_json_file(file_path):#Load the ECI values from the JSON file
    with open(file_path, 'r') as file:
        json_data = json.load(file)
    return json_data 
def run_mc(cell,T):#Run the Monte Carlo Simulation for the cell 
	internal_energy=[]
	mc=Montecarlo(cell,T)
	Au_conc.append(mc.count_atoms()['Au']/tot)
	for mc_step in mc.irun(500):
		if mc_step.step%10==0:
			internal_energy.append(mc_step.energy)
	mean_energy.append(np.mean(internal_energy/tot))
	formula.append(cell.get_chemical_formula(empirical=False))
	temperature.append(T)
def find_pure_energy(atoms, string,T):#Finding the pure energy of the lattices
	for i in range(len(atoms)):
			atoms[i].symbol=string
	atoms = attach_calculator(settings, atoms=atoms, eci=eci)
	formula.append(atoms.get_chemical_formula(empirical=False))
	mean_energy.append(atoms.get_total_energy()/tot)
	temperature.append(T)
	if(string=='Au'):
		Au_conc.append(1)
		return (atoms.get_total_energy()/tot)
	else:
		Au_conc.append(0)
		return(atoms.get_total_energy()/tot)
def get_free_energy(U,T):
	Kb=1.380649e-23
	beta=1/(Kb*T)
	N=1000
#	Z=math.exp(-U*beta*N)
	F=(Kb*T*U*beta*N)
	free_energy.append(F)
file_path = 'eci_l1.json'  
json_dict = load_json_file(file_path)
eci=json_dict
chem_pot={}
chem_pot['c1_0']=eci['c1_0']#chemical potentional for SGC Monte Carlo 
conc=Concentration(basis_elements=[['Au','Cu']])
conc.set_conc_formula_unit(formulas=["Au<x>Cu<1-x>"], variable_range={"x": (0, 1)})
settings=CEBulk(crystalstructure='fcc',
               a=3.8,
               supercell_factor=27,
               concentration=conc,
               db_name="aucu.db",
               max_cluster_dia=[6.0,5,5])
formula=[]
mean_energy=[]
free_energy=[]
size=(10,10,10)
tot=np.product(size)
atoms = connect('aucu.db').get(id=1).toatoms()*size
temp=[500]
temperature=[]
Au_conc=[]
Au_pure,Cu_pure=0,0
c=0
formation=[]
for T in temp:#for different temperature ranges
	i=1
	Au_pure=find_pure_energy(atoms,'Au',T)	
	while(i>=1 and i<=len(atoms)-1):
		cell = connect('aucu.db').get(id=1).toatoms()*size
		for j in range (0,i):
			cell[j].symbol='Cu'
		cell = attach_calculator(settings, atoms=cell, eci=eci)
		run_mc(cell,T)
		i=i+10
	Cu_pure=find_pure_energy(atoms,'Cu',T)
for i in range(len(mean_energy)):
	u=((mean_energy[i]-(Au_conc[i]*Au_pure)-((1-Au_conc[i])*Cu_pure))*96.484934)
	formation.append(u)
	get_free_energy(mean_energy[i],temperature[i])
for i in range (len(free_energy)):
	free_energy[i]=(free_energy[i]-(Au_conc[i]*free_energy[0])-((1-Au_conc[i])*free_energy[-1]))*96.484934
plt.plot(Au_conc,formation)
plt.xlabel('Au Conc')
plt.ylabel('Formation Energy (kJ/mole)')
plt.axhline(y = 0.0, color = 'r', linestyle = '-')
plt.legend()
plt.savefig("Formation1")
