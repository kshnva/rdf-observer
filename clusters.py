import pandas as pd 
import numpy as np 
from sklearn.model_selection import train_test_split 
from clease.structgen import NewStructures
from clease.settings import Concentration,CEBulk
from ase.db import connect 
from clease.structgen import NewStructures
import os
from clease import Evaluate
import clease.plot_post_process as pp
import matplotlib.pyplot as plt
from ase.calculators.emt import EMT
from ase.db import connect
from clease.tools import update_db
import json 
def generate_probe():
	ns = NewStructures(settings, generation_number=1, struct_per_gen=10)
	ns.generate_probe_structure()
def eval():
	eva = Evaluate(settings=settings, scoring_scheme='k-fold', nsplits=10)
	eva.set_fitting_scheme(fitting_scheme='l1')
	alpha = eva.plot_CV(alpha_min=1E-7, alpha_max=1.0, num_alpha=50)
	eva.set_fitting_scheme(fitting_scheme='l1', alpha=alpha)
	eva.fit()
	fig = pp.plot_fit(eva)
	plt.savefig("L1Regularization.png")
	plt.show()
	fig = pp.plot_eci(eva)
	plt.show()
	eva.save_eci(fname='eci_l1')
conc=Concentration(basis_elements=[['Au','Cu']])
conc.set_conc_formula_unit(formulas=["Au<x>Cu<1-x>"], variable_range={"x": (0, 1)})#Setting the concentration ranges
settings=CEBulk(crystalstructure='fcc',
               a=3.8,
               supercell_factor=27,
               concentration=conc,
               db_name="aucu.db",
               max_cluster_dia=[6.0,5,5])#Setting for the Cluster Expansion
from clease.structgen import NewStructures
ns = NewStructures(settings, generation_number=0, struct_per_gen=10)
ns.generate_initial_pool()
generate_probe()
calc = EMT()
db_name = "aucu.db"
db = connect(db_name)
for row in db.select(converged=False):
  atoms = row.toatoms()
  atoms.calc = calc
  atoms.get_potential_energy()
  update_db(uid_initial=row.id, final_struct=atoms, db_name=db_name)
eval()
