# -*- coding: utf-8 -*-

# Libraries
import os
import sys
from cobra.io import read_sbml_model, write_sbml_model

# Determine the current file's directory and the project root.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")

# Add the project root to sys.path to access top-level folders like 'functions' and 'models'
if project_root not in sys.path:
    sys.path.append(project_root)

from functions.functions_network_consistency import *
from functions.functions_merge_metabolic_networks import *



# Input parameters
GEM1 = os.path.join(project_root, "models", "THG-beta2.xml")
GEM2 = os.path.join(project_root, "models", "Human Database.xml")
Output = os.path.join(project_root, "models", "THG-2023-02-25.xml")


# Loading the input data
network_1 = read_sbml_model(GEM1) # target network (i.e. network from GEM)
network_2 = read_sbml_model(GEM2) # source network (i.e. netwokr from DBs)


# Step 1
network_3 = network_metabolites_merge_3 (network_1, network_2)
network_3_metabolites = network_3[0]
equivalent_metabolites = network_3[1]


# Step 2
network_3_genes = network_genes_merge_2 (network_3_metabolites, network_2)
network_3_genes_model = network_3_genes[0]
equivalent_genes = network_3_genes[1]


# Step 3
network_3_reactions = network_reactions_merge_6 (network_3_genes_model, network_2, equivalent_metabolites, equivalent_genes)
network_3_reactions_model = network_3_reactions[0]
updated_reactions_from_network_2 = network_3_reactions[1]
new_reactions_from_network_2 = network_3_reactions[2]


# Step 4 : Sanitize model
new_network_4_metabolites = delete_isolated_metabolites(network_3_reactions[0])
new_network_4_reactions = delete_not_used_reactions(new_network_4_metabolites[0])
new_network_4_genes = delete_not_used_genes(new_network_4_reactions)

write_sbml_model(new_network_4_reactions, Output)
