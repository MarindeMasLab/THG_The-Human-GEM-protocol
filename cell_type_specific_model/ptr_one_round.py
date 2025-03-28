# Import libs
from __future__ import print_function
import cobra
import cobra.io
import cobra.core
from cobra.flux_analysis import gapfill
import sys
import re
import copy
import pandas as pd 
import numpy as np
import time
from multiprocessing import freeze_support
import concurrent.futures
import pickle
import memote.support.consistency as consistency
from memote.utils import annotate, wrapper
import pdb
import os
import multiprocessing

# Determine the current file's directory and the project root.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")
# Add the project root to sys.path to access top-level folders like 'functions' and 'models'
if project_root not in sys.path:
    sys.path.append(project_root)

# Define functions that will be used by "add_transport_and_sink_reactions" function
def ident_final_dead_end_metabolites(model):
    """ 
    Identify dead-end metabolites that are at the edges of an isolated path, also incluying those in an isolated futile cycle

    Inputs
    ----------
        - model: Model name (SBML format)
    Output
    ----------
        - DeadEnd_M: dead-end metabolites
    """
    # Create a copy of the model and work with the copy
    model_copy = model.copy()
    # Identify Blocked reactions
    blockedReactions = cobra.flux_analysis.find_blocked_reactions(model_copy) 

    # Identify dead-end metabolites
    DeadEnd_S = []
    DeadEnd_P = []
    DeadEnd_M = [] #will only contain metabolites in isolated futile cycles       
    for x in blockedReactions:
        IthBR = getattr(model_copy.reactions, x)
        if abs(IthBR.upper_bound)+abs(IthBR.lower_bound) != 0: # if the reaction has an upper and lower bound different from 0
            for y in range(0,len(IthBR.metabolites)): 
                YthMet = list(IthBR.metabolites.items())[y][0]
                if len(YthMet.reactions) == 1:
                    YthSth = list(IthBR.metabolites.items())[y][1]
                    if YthSth < 0: DeadEnd_S.append(YthMet.id)
                    else:  DeadEnd_P.append(YthMet.id)                
    print("Original DeadEnd_S len: " ,len(DeadEnd_S))
    print("Original DeadEnd_P len: " ,len(DeadEnd_P))


    #To detect isolated futile cycles
    allR = [] #list with all the reactions in the model
    for x in model_copy.reactions: allR.append(x.id)

    blocked_reactions_set = set(blockedReactions)

    # Find the difference between the two sets to obtain the unblocked Rs
    unblocked_reactions = set(allR) - blocked_reactions_set

    print("Length of all reactions: ", len(allR))
    print("Length of blocked reactions: ", len(blocked_reactions_set))
    print("Length of unblocked reactions: ", len(unblocked_reactions))

    unblocked_reactions_list = list(unblocked_reactions)         
    lowerbounds=[] #list with all the lower bounds of unblocked Reactions
    upperbounds=[] #list with all the upper bounds of unblocked Reactions   

    for IthR in unblocked_reactions: #not analyzing blocked reactions

        IthR = getattr(model_copy.reactions, IthR)
        if IthR in model_copy.exchanges or IthR in model_copy.demands or IthR in model_copy.sinks:
            IthR.upper_bound=100.0
            if IthR.reversibility:
                IthR.lower_bound=-100.0 #Reversible reaction
            else:
                IthR.lower_bound=0.0 #Irreversible reaction
        else:
            IthR.upper_bound=1000.0
            if IthR.reversibility:
                IthR.lower_bound=-1000.0 #Reversible reaction
            else:
                IthR.lower_bound=0.0 #Irreversible reaction
        lowerbounds.append(IthR.lower_bound)
        upperbounds.append(IthR.upper_bound)
                
    # Perform FVA and get the minimum and maximum fluxes
    fva_result = cobra.flux_analysis.flux_variability_analysis(model_copy, reaction_list=unblocked_reactions_list) #excluding BR
    
    # Convert bounds to pandas DataFrame with reaction IDs as index
    lower_bounds_df = pd.DataFrame(lowerbounds, index=fva_result.index, columns=['lower_bound'])
    upper_bounds_df = pd.DataFrame(upperbounds, index=fva_result.index, columns=['upper_bound'])
    # Find reactions where min flux equals lower bound and max flux equals upper bound, i.e., reactions that operate at their limits

    # Comparison conditions using np.isclose()
    condition_min = np.isclose(fva_result['minimum'], lower_bounds_df['lower_bound'])
    condition_max = np.isclose(fva_result['maximum'], upper_bounds_df['upper_bound'])

    # Apply the conditions to filter the DataFrame and get the index as a list

    selected_reactions = fva_result[condition_min & condition_max].index.tolist()

    # Extract metabolites from selected reactions: adding all the metabolites involved in selected reactions
    selected_metabolites = set(metabolite for reaction_id in selected_reactions for metabolite in model_copy.reactions.get_by_id(reaction_id).metabolites)
    selected_metabolites = list(selected_metabolites)

    # Initialize two dictionaries to store the metabolites that act as reactants and products
    reactants = {}
    products = {}

    # Loop over the selected reactions
    for reaction_id in selected_reactions:
        reaction = model_copy.reactions.get_by_id(reaction_id)
        for metabolite, coefficient in reaction.metabolites.items():
            # If the coefficient is negative, the metabolite is a reactant
            if coefficient < 0:
                reactants[metabolite] = reactants.get(metabolite, 0) + 1
            # If the coefficient is positive, the metabolite is a product
            elif coefficient > 0:
                products[metabolite] = products.get(metabolite, 0) + 1


    selected_metabolites_2R_ids = []

    for metabolite in selected_metabolites:
        if reactants.get(metabolite, 0) == 1 and products.get(metabolite, 0) == 1 and len(metabolite.reactions)== 2:
            list_R=list(metabolite.reactions)
            list_R1_S_id = sorted([reactant.id for reactant in list_R[0].reactants])
            list_R1_P_id = sorted([product.id for product in list_R[0].products])
            list_R2_S_id = sorted([reactant.id for reactant in list_R[1].reactants])
            list_R2_P_id = sorted([product.id for product in list_R[1].products])


            if list_R1_S_id==list_R2_P_id and list_R1_P_id==list_R2_S_id:
                #Check that all the metabolites in the reactions are in only 2 reactions
                IFC=True

                for metabolite in list_R1_S_id + list_R1_P_id:
                    if len(model_copy.metabolites.get_by_id(metabolite).reactions)!=2:
                        IFC=False

                if IFC:
                    for metabolite in list_R1_S_id + list_R1_P_id:
                        selected_metabolites_2R_ids.append(metabolite)

    DeadEnd_M.extend(selected_metabolites_2R_ids)
    # Remove duplicates
    DeadEnd_M = list(set(DeadEnd_M))

    print("DeadEnd_M len: " ,len(DeadEnd_M))
    del(model_copy)
              
    return DeadEnd_M, DeadEnd_S, DeadEnd_P

def modify_reaction_bounds(model):
    """
    Modify the bounds of non-biomass reactions in a model if their bounds are 0,0.
    
    This function iterates through all the reactions in the given model and checks if a reaction is not biomass and its bounds are 0,0.
    If these conditions are met, the function modifies the bounds of the reaction based on its reversibility.
    
    Inputs
    ----------
        model: A Cobra model representing a metabolic network.
        
    Output
    -------
        model: The modified Cobra model object with updated reaction bounds.
    """
    for reaction in model.reactions:
        # Check if the reaction is not biomass and its bounds are 0,0
        if 'biomass' not in reaction.name.lower() and all(substrate.name.lower() != 'biomass' for substrate in reaction.reactants) and reaction.lower_bound == 0 and reaction.upper_bound == 0:
            # If the reaction is reversible, modify upper and lower boundaries
            if reaction.reversibility:
                reaction.lower_bound = -1000
                reaction.upper_bound = 1000
            # If the reaction is irreversible, modify upper boundary
            else:
                reaction.upper_bound = 1000
    return model

def search_metabolite_locations(metid, metabolites_in_model, model):
    """ 
    creates a dictionary with the metabolite and the compartments where the metabolite can be found 

    Inputs
    ----------
        - metID: Metabolite ID
        - model: cobra model
        - metabolites_in_model: list of metabolites in the model

    Output
    ----------
        - metabolite_in_compartments: list of metabolite cobra objects describing the same metabolite in different compartments
    """    
    metabolite_in_compartments = []
    metabolite_in_model = getattr(model.metabolites, metid) # retrieves the metabolite object for the given metid from the model.
    metabolite_compartment = metabolite_in_model.compartment
    metabolite_id = re.sub(metabolite_compartment, '', metabolite_in_model.id) # replaces the metabolite_compartment with an empty string in the metabolite_in_model.id. This action isolates the core ID of the metabolite, allowing comparisons across different compartments.
    for x in metabolites_in_model:
        x_metabolite = getattr(model.metabolites, x)
        x_compartment = x_metabolite.compartment
        x_id = re.sub(x_compartment, '', x_metabolite.id)
        if x_id == metabolite_id: metabolite_in_compartments = metabolite_in_compartments + [x_metabolite] #Check if the modified ID matches the ID of the target metabolite (metabolite_id). If it does, this indicates that this metabolite (despite being in a different compartment) is the same as the one we're searching for.
    return metabolite_in_compartments #A list of COBRA metabolite objects. Each object represents an occurrence of the specified metabolite (metid) in a different compartment of the model.


def identify_connected_species(species, allowed_connections, model):
    """ 
    Identify those compartments that can be connected with a transport reaction and eliminate the not allowed connections

    Inputs
    ----------
        - species: dictionary containing lists of the same metabolite (metabolite cobra object) in different compartments 
        - allowed_connected_compartment: list with pairs of compartment ids allowed to be connected
        - model: cobra model
    Output
    ----------
        - allowed_connected_species: pairs of metabolites allowed to be connected by a transport reaction since they are adjacent
    """
    allowed_connected_species = []
    metabolites_in_model = [x.id for x in model.metabolites] # Metabolite IDs
    for x in species:
        t = search_metabolite_locations(x, metabolites_in_model, model)
        for p in range(0,len(t)):
            for q in range(p+1,len(t)): # iterates over the remaining metabolite objects in t, starting from the next element after p to avoid duplicate comparisons and self-pairing
                if (t[q].id in x or t[p].id in x) and ((t[p].compartment,t[q].compartment) in allowed_connections or (t[q].compartment,t[p].compartment) in allowed_connections): #checks if any of the metabolites in the pair is in the soecies list and if the compartments of the two metabolites are in the allowed_connections list
                    allowed_connected_species = allowed_connected_species + [[t[p], t[q]]]      
    return allowed_connected_species

def product_transport_metabolite_to(metabolite_pair, ref_rxn_id, DeadEnd_P):

    """
    Given a metabolite pair, create the transport reaction from the metabolite in the Dead_End_P list to the other metabolite.
    
    Return the transport reaction and the transported metabolite's ID. The transported metabolite is also added to
    metabolites_list and, if provided, metabolites_by_size.
    
    Inputs
    ----------
    metabolite_pair : list of metabolites in the reaction (reaction objects)
    ref_rxn_id : reference number to generate automatic reaction IDs
    DeadEnd_P : list of dead-end products

    Output
    -------
    transport_reaction : new transport reaction (cobra reaction object)
    """

    transport_reaction = cobra.Reaction("MAR"+str(ref_rxn_id))
    transport_reaction.name = "Transport "+metabolite_pair[0].name+" from "+[y.compartment for y in metabolite_pair][0]+" to "+[y.compartment for y in metabolite_pair][1]
    transport_reaction.subsystem = 'Putative Transport Reaction'
    transport_reaction.lower_bound = 0
    transport_reaction.upper_bound = 1000.0
    transport_reaction.notes = {'Confidence Level': '4', 'AUTHORS': 'Putative transport reaction automatically generated'}

    #Check if any of the metabolites in the pair is in the DeadEnd_P list
    # if metabolite_pair[0].id in DeadEnd_P and metabolite_pair[1].id in DeadEnd_P:
    #     print("Both mets are in the DeadEnd_P list, no TR will be created because no feasible solution is possible")
    if metabolite_pair[0].id in DeadEnd_P: #then, this is the met that has the sink met associated
        transport_reaction.add_metabolites({metabolite_pair[0]:-1.0, metabolite_pair[1]:1.0})  #one consumed, one produced. 
        transport_reaction.gene_reaction_rule = ''
    elif metabolite_pair[1].id in DeadEnd_P: #then, this is the met that has the sink met associated
        transport_reaction.add_metabolites({metabolite_pair[1]:-1.0, metabolite_pair[0]:1.0})  #one consumed, one produced. 
        transport_reaction.gene_reaction_rule = ''

    return transport_reaction

def substrate_transport_metabolite_to(metabolite_pair, ref_rxn_id, DeadEnd_S):

    """
    Given a metabolite pair, create the transport reaction from the other metabolite to the metabolite in the Dead_End_S list.
    
    Return the transport reaction and the transported metabolite's ID. The transported metabolite is also added to
    metabolites_list and, if provided, metabolites_by_size.
    
    Inputs
    ----------
    metabolite_pair : list of metabolites in the reaction (reaction objects)
    ref_rxn_id : reference number to generate automatic reaction IDs
    DeadEnd_S : list of dead-end substrates

    Output
    -------
    transport_reaction : new transport reaction (cobra reaction object)
    """

    transport_reaction = cobra.Reaction("MAR"+str(ref_rxn_id))
    transport_reaction.name = "Transport "+metabolite_pair[0].name+" from "+[y.compartment for y in metabolite_pair][1]+" to "+[y.compartment for y in metabolite_pair][0]
    transport_reaction.subsystem = 'Putative Transport Reaction'
    transport_reaction.lower_bound = 0
    transport_reaction.upper_bound = 1000.0
    transport_reaction.notes = {'Confidence Level': '4', 'AUTHORS': 'Putative transport reaction automatically generated'}

    #Check if any of the metabolites in the pair is in the DeadEnd_S list
    # if metabolite_pair[0].id in DeadEnd_S and metabolite_pair[1].id in DeadEnd_S:
    #     print("Both mets are in the DeadEnd_S list, no TR will be created because no feasible solution is possible")
    if metabolite_pair[0].id in DeadEnd_S: #then, this is the met that has the sink met associated
        transport_reaction.add_metabolites({metabolite_pair[1]:-1.0, metabolite_pair[0]:1.0})  #one consumed, one produced. 
        transport_reaction.gene_reaction_rule = ''
    elif metabolite_pair[1].id in DeadEnd_S: #then, this is the met that has the sink met associated
        transport_reaction.add_metabolites({metabolite_pair[0]:-1.0, metabolite_pair[1]:1.0})  #one consumed, one produced. 
        transport_reaction.gene_reaction_rule = ''

    return transport_reaction

def wrapped_gapfill (model_target, model_source, lower_bound, penalties, exchange_reactions, demand_reactions):
    ''''gapfill with decorator'''    
    return gapfill(model_target, model_source, lower_bound, penalties, exchange_reactions, demand_reactions)

# Already defined functions
def jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

# Identify dead-end metabolites and add putative transport reactions, only to dead end S and P
def add_transport_and_sink_reactions (model_original, compartments, sheetname):
    """
    For a given GEM:
        - Identify the dead-end metabolites that are associated to only one reaction as either substrate or product
        - Identify the dead-metabolites in other compartments
        - Add a transport reaction between the dead-end metabolites and their equivalent in other locations (if not exists already in the model)
        - Add sink reactions to the dead-end metabolites (if not exists already in the model)
    
    Inputs
    ----------
    model_original : metabolic Cobra model
    compartments : List of compartments stored in a excel file and their interactions. 
                   Compartment interaction is defined by a square binary matrix of size nxn where n is the number of compartments
                    - 1: the compartment in the ith row and jth column are adjacent so can be connected by a transport reaction
                    - 0: the compartment in the ith row and jth column are not adjacent so cannot be connected by a transport reaction 
                   Additional row on top of the matrix and addition column on the left containing the compartment names 
    sheetname : name of the sheet where the above mentioned table is                 
    
    Output
    -------
    model_sink : cobra model with only sink reactions added
    model_transport : cobra model with only transport reactions added
    new_transport_reactions : list of newly added putative transport reactions 
    DeadEnd_M : list of dead-end metabolites
    DeadEnd_S : list of dead-end substrates
    DeadEnd_P : list of dead-end products
    """
    
    ## Models
    
    print('Loading model')
    model = model_original.copy() # Contains only sink reactions
    
    ## Compartments
    
    compartment_connection= pd.read_excel(compartments, sheet_name=sheetname ,skiprows=[0,1],header=None) # Adjacent compartments where transport reactions can be added
    compartments = compartment_connection[0]
    allowed_connections = [] # Compartments that can be connected
    for y in range(1, np.shape(compartment_connection)[1]):
        for x in range(y,len(compartment_connection[0])):
            if compartment_connection[y][x] == 1: allowed_connections.append((compartment_connection[0][x],compartment_connection[0][y-1]))
    

    ## Metabolites
    
    print('Identify dead-end metabolites')
    DeadEnd_M, DeadEnd_S, DeadEnd_P = ident_final_dead_end_metabolites(model)


    allowed_dead_end_substrate_connection = identify_connected_species(DeadEnd_S, allowed_connections, model) # Identify Dead-End substrates and metabolites in other compartments that can be potentially connected
    allowed_dead_end_product_connection = identify_connected_species(DeadEnd_P, allowed_connections, model) # Identify Dead-End products and metabolites in other compartments that can be potentially connected

    print("Connected species number for DeadEnd_S: ", len(allowed_dead_end_substrate_connection))
    print("Connected species number for DeadEnd_P: ", len(allowed_dead_end_product_connection))

    
    ## Reference reaction for annotation    
    
    ref_id = max([int(re.findall("^[A-Z,a-z,_,-]*([0-9]*)", x.id)[0]) for x in model.reactions])
    if not ref_id:
        ref_id = max([int(re.findall("^([0-9]*)[A-Z,a-z,_,-]", x.id)[0]) for x in model.reactions])
    r = re.sub(re.findall("^[A-Z,a-z,_,-]*([0-9]*)", model.reactions[0].id)[0],str(ref_id),model.reactions[0].id)
    pattern_replace = ref_id #determining a reference ID based on the numeric parts of reaction IDs in the model. It considers two patterns: one where the numeric part is at the end of the ID and another where it is at the beginning. The maximum of these numeric parts is chosen as the reference ID, and it is used for pattern replacement in the subsequent part of the code. The purpose is likely to ensure consistent and unique annotations within the function
    
    ## Add reactions
    
    # Add sink reactions substrates
    
    print('Adding sink reactions')
    new_reverse_reactions = []
    sink_reactions = []
    associated_sink_metabolites = []
    reaction_associated_to_dead_end_metabolites = [str(model.metabolites.get_by_id(x).reactions).split(" ")[1] for x in DeadEnd_S+DeadEnd_P] # list of reactions to be modified
    model.compartments.update({'sink': 'Sink'}) # model compartments are updated to include a new compartment named 'sink' with the label 'Sink'
    l = len(reaction_associated_to_dead_end_metabolites)
    count=1
    print('Number of reactions to be modified:', len(set(reaction_associated_to_dead_end_metabolites)))
    print('Without the set:', len(reaction_associated_to_dead_end_metabolites))
    
    for reaction_to_dead_end in set(reaction_associated_to_dead_end_metabolites):
        current_reaction = model.reactions.get_by_id(reaction_to_dead_end)
        if current_reaction.boundary == True: #if the reaction is a boundary reaction, create a new reaction and work on it            
            ref_id = ref_id +1
            ref_rxn = re.sub(str(pattern_replace), str(ref_id), r) 
            new_reaction = cobra.Reaction(ref_rxn)
            new_reaction.id = ref_rxn     
            new_reaction.name = current_reaction.name
            new_reaction.subsystem = 'Non-boundary reaction'
            new_reaction.lower_bound = -1000.0
            new_reaction.upper_bound = 1000.0
            if current_reaction.reactants:
                metabolite_in_reaction = current_reaction.reactants[0]
            else:
                metabolite_in_reaction = current_reaction.products[0]
            new_reaction.add_metabolites({getattr(model.metabolites, metabolite_in_reaction.id):-1.0})           
            new_reaction.notes = {'Confidence Level': '4', 'AUTHORS': 'Sink transport reaction automatically generated'}
            new_reaction.gene_reaction_rule = ''
            model.add_reactions([new_reaction.copy()])
            reaction_to_dead_end = ref_rxn
        reactions_in_loop = [reaction_to_dead_end]
        count = count + 1

        if current_reaction.lower_bound < 0:
            model.reactions.get_by_id(reaction_to_dead_end).lower_bound = 0 # make the reaction irreversible
            current_reaction.lower_bound=0
            backwards_reaction = current_reaction.reaction.split(' --> ')[1]+' --> '+current_reaction.reaction.split(' --> ')[0]            
            reverse_reaction = [x for x in model.reactions if x.reaction == backwards_reaction] # Check if the reverse reaction is already in the model
            if not reverse_reaction: # create a new reaction in the reverse direction
                reverse_reaction=""                
                reverse_reaction = current_reaction.copy()
                reverse_reaction.id = reverse_reaction.id+'_reverse'
                reverse_reaction.name = reverse_reaction.name+' Reverse'
                model.add_reactions([reverse_reaction.copy()]) # add reverse reaction to the model
                model.reactions.get_by_id(reverse_reaction.id).reaction = backwards_reaction
                new_reverse_reactions = new_reverse_reactions+[reverse_reaction.id] #to concatenate str to list, we need to make the str a list
                reactions_in_loop = reactions_in_loop+[reverse_reaction.id]
            else:
                reactions_in_loop = reactions_in_loop+[reverse_reaction[0].id] #here, RR would be a list so we need to extract it this way
           
        associated_metabolites = current_reaction.reactants+current_reaction.products

        for x in associated_metabolites:
            ref_id = ref_id +1
            ref_rxn = re.sub(str(pattern_replace), str(ref_id), r)
            metabolites_in_model = [x.id for x in model.metabolites] 
            if not x.id+'_sink' in metabolites_in_model and not re.findall('_sink', x.id) and not x.id+'_sink' in associated_sink_metabolites: # if the sink version of the metabolite in not in the model and the metabolite is not a sink metabolite
                if x.id in DeadEnd_S:
                    sink_reaction = cobra.Reaction(ref_rxn)
                    sink_reaction.id = ref_rxn     
                    sink_reaction.name = "Sink "+x.name
                    sink_reaction.subsystem = 'Putative Sink Reaction'
                    sink_reaction.lower_bound = 0
                    sink_reaction.upper_bound = 0.0
                    #If metabolite is a dead end substrate, it is a product in the sink reaction
                    sink_reaction.add_metabolites({getattr(model.metabolites, x.id):1.0})
                    sink_reaction.notes = {'Confidence Level': '4', 'AUTHORS': 'Sink transport reaction automatically generated'}
                    sink_reaction.gene_reaction_rule = ''
                    model.add_reactions([sink_reaction.copy()])
                    sink_reactions = sink_reactions + [sink_reaction.id]

                elif x.id in DeadEnd_P:
                    sink_reaction = cobra.Reaction(ref_rxn)
                    sink_reaction.id = ref_rxn     
                    sink_reaction.name = "Sink "+x.name
                    sink_reaction.subsystem = 'Putative Sink Reaction'
                    sink_reaction.lower_bound = 0
                    sink_reaction.upper_bound = 0.0
                    #If metabolite is a dead end product, it is a substrate in the sink reaction
                    sink_reaction.add_metabolites({getattr(model.metabolites, x.id):-1.0})           
                    sink_reaction.notes = {'Confidence Level': '4', 'AUTHORS': 'Sink transport reaction automatically generated'}
                    sink_reaction.gene_reaction_rule = ''
                    model.add_reactions([sink_reaction.copy()])
                    sink_reactions = sink_reactions + [sink_reaction.id]
                      
                    
    model_sink = model.copy() # model with only sink reactions added
    model_transport = cobra.Model("model_transport") # model to store transport reactions

    # Add transport reactions

    print('Adding transport reactions')
    new_transport_reactions = []
    metabolites_in_1_compound_reaction = list([x.metabolites for x in model.reactions if len(x.metabolites) == 2])

    #For DeadEnd_S:
    for metabolite_pair in allowed_dead_end_substrate_connection:

        # Create a reaction transporting the nth metabolite between compartments

        ref_id = ref_id +1
        ref_rxn = re.sub(str(pattern_replace), str(ref_id), r)
        new_transport_reaction = copy.deepcopy(substrate_transport_metabolite_to(metabolite_pair, ref_rxn, DeadEnd_S))

        # Check if the reaction is already in the model

        for y in range(0,len(metabolites_in_1_compound_reaction)):
            new_metabolites = [x.id for x in list(new_transport_reaction.metabolites)]
            yth_metabolites = [x.id for x in list(metabolites_in_1_compound_reaction[y])]
            overlapping_test = jaccard(new_metabolites, yth_metabolites)
            if overlapping_test == 1: # If a coincidence is found stop the loop
                break
        if overlapping_test != 1:   
            model.add_reactions([new_transport_reaction]) # if no coincidences have found then add the reaction to the model
            model_transport.add_reactions([new_transport_reaction]) # if no coincidences have found then add the reaction to the model           
            new_transport_reactions = new_transport_reactions + [new_transport_reaction]

    #For DeadEnd_P:
    for metabolite_pair in allowed_dead_end_product_connection:

        # Create a reaction transporting the nth metabolite between compartments

        ref_id = ref_id +1
        ref_rxn = re.sub(str(pattern_replace), str(ref_id), r)
        new_transport_reaction = copy.deepcopy(product_transport_metabolite_to(metabolite_pair, ref_rxn, DeadEnd_P))

        # Check if the reaction is already in the model

        for y in range(0,len(metabolites_in_1_compound_reaction)):
            new_metabolites = [x.id for x in list(new_transport_reaction.metabolites)]
            yth_metabolites = [x.id for x in list(metabolites_in_1_compound_reaction[y])]
            overlapping_test = jaccard(new_metabolites, yth_metabolites)
            if overlapping_test == 1: # If a coincidence is found stop the loop
                break
        if overlapping_test != 1:   
            model.add_reactions([new_transport_reaction]) # if no coincidences have found then add the reaction to the model
            model_transport.add_reactions([new_transport_reaction]) # if no coincidences have found then add the reaction to the model           
            new_transport_reactions = new_transport_reactions + [new_transport_reaction]

    new_transport_reactions = list(set(new_transport_reactions))
    del(model)
    return model_sink, model_transport, new_transport_reactions, DeadEnd_S, DeadEnd_P, DeadEnd_M

# Identify dead-end metabolites and add putative transport reactions, now only to the dead-end metabolites in futile cycles
def add_transport_and_sink_reactions_futile_cycles (model_original, compartments, sheetname):
    """
    For a given GEM:
        - Identify the dead-end metabolites that are associated to only one reaction as either substrate or product
        - Identify the dead-metabolites in other compartments
        - Add a transport reaction between the dead-end metabolites and their equivalent in other locations (if not exists already in the model)
        - Add sink reactions to the dead-end metabolites (if not exists already in the model)
    
    Inputs
    ----------
    model_original : metabolic Cobra model
    compartments : List of compartments stored in a excel file and their interactions. 
                   Compartment interaction is defined by a square binary matrix of size nxn where n is the number of compartments
                    - 1: the compartment in the ith row and jth column are adjacent so can be connected by a transport reaction
                    - 0: the compartment in the ith row and jth column are not adjacent so cannot be connected by a transport reaction 
                   Additional row on top of the matrix and addition column on the left containing the compartment names 
    sheetname : name of the sheet where the above mentioned table is                 
    
    Output
    -------
    model_sink : cobra model with only sink reactions added
    model_transport : cobra model with only transport reactions added
    new_transport_reactions : list of newly added putative transport reactions 
    DeadEnd_M : list of dead-end metabolites
    DeadEnd_S : list of dead-end substrates
    DeadEnd_P : list of dead-end products
    """
    
    ## Models
    
    print('Loading model')
    model = model_original.copy() # Contains only sink reactions
    
    ## Compartments
    
    compartment_connection= pd.read_excel(compartments, sheet_name=sheetname ,skiprows=[0,1],header=None) # Adjacent compartments where transport reactions can be added
    compartments = compartment_connection[0]
    allowed_connections = [] # Compartments that can be connected
    for y in range(1, np.shape(compartment_connection)[1]):
        for x in range(y,len(compartment_connection[0])):
            if compartment_connection[y][x] == 1: allowed_connections.append((compartment_connection[0][x],compartment_connection[0][y-1]))


    ## Metabolites
    
    print('Identify dead-end metabolites')
    DeadEnd_M, DeadEnd_S, DeadEnd_P = ident_final_dead_end_metabolites(model)


    allowed_futilecycles= identify_connected_species(DeadEnd_M, allowed_connections, model) # Identify Dead-End metabolites in futile cycles and metabolites in other compartments that can be potentially connected

    print("Connected species number for DeadEnd_M: ", len(allowed_futilecycles))
    
    ## Reference reaction for annotation    
    
    ref_id = max([int(re.findall("^[A-Z,a-z,_,-]*([0-9]*)", x.id)[0]) for x in model.reactions])
    if not ref_id:
        ref_id = max([int(re.findall("^([0-9]*)[A-Z,a-z,_,-]", x.id)[0]) for x in model.reactions])
    r = re.sub(re.findall("^[A-Z,a-z,_,-]*([0-9]*)", model.reactions[0].id)[0],str(ref_id),model.reactions[0].id)
    pattern_replace = ref_id #determining a reference ID based on the numeric parts of reaction IDs in the model. It considers two patterns: one where the numeric part is at the end of the ID and another where it is at the beginning. The maximum of these numeric parts is chosen as the reference ID, and it is used for pattern replacement in the subsequent part of the code. The purpose is likely to ensure consistent and unique annotations within the function
    
    ## Add reactions
    
    # Add sink reactions substrates
    
    print('Adding sink reactions')
    new_reverse_reactions = []
    sink_reactions = []
    associated_sink_metabolites = []
    reaction_associated_to_dead_end_metabolites = [str(model.metabolites.get_by_id(x).reactions).split(" ")[1] for x in DeadEnd_M] # list of reactions to be modified
    model.compartments.update({'sink': 'Sink'}) # model compartments are updated to include a new compartment named 'sink' with the label 'Sink'
    l = len(reaction_associated_to_dead_end_metabolites)
    count=1
    print('Number of reactions to be modified:', len(set(reaction_associated_to_dead_end_metabolites)))
    print('Without the set:', len(reaction_associated_to_dead_end_metabolites))
    
    for reaction_to_dead_end in set(reaction_associated_to_dead_end_metabolites):
        current_reaction = model.reactions.get_by_id(reaction_to_dead_end)
        if current_reaction.boundary == True: #if the reaction is a boundary reaction, create a new reaction and work on it            
            ref_id = ref_id +1
            ref_rxn = re.sub(str(pattern_replace), str(ref_id), r) 
            new_reaction = cobra.Reaction(ref_rxn)
            new_reaction.id = ref_rxn     
            new_reaction.name = current_reaction.name
            new_reaction.subsystem = 'Non-boundary reaction'
            new_reaction.lower_bound = -1000.0
            new_reaction.upper_bound = 1000.0
            if current_reaction.reactants:
                metabolite_in_reaction = current_reaction.reactants[0]
            else:
                metabolite_in_reaction = current_reaction.products[0]
            new_reaction.add_metabolites({getattr(model.metabolites, metabolite_in_reaction.id):-1.0})           
            new_reaction.notes = {'Confidence Level': '4', 'AUTHORS': 'Sink transport reaction automatically generated'}
            new_reaction.gene_reaction_rule = ''
            model.add_reactions([new_reaction.copy()])
            reaction_to_dead_end = ref_rxn
        reactions_in_loop = [reaction_to_dead_end]
        count = count + 1

        if current_reaction.lower_bound < 0:
            model.reactions.get_by_id(reaction_to_dead_end).lower_bound = 0 # make the reaction irreversible
            current_reaction.lower_bound=0
            backwards_reaction = current_reaction.reaction.split(' --> ')[1]+' --> '+current_reaction.reaction.split(' --> ')[0]            
            reverse_reaction = [x for x in model.reactions if x.reaction == backwards_reaction] # Check if the reverse reaction is already in the model
            if not reverse_reaction: # create a new reaction in the reverse direction
                reverse_reaction=""                
                reverse_reaction = current_reaction.copy()
                reverse_reaction.id = reverse_reaction.id+'_reverse'
                reverse_reaction.name = reverse_reaction.name+' Reverse'
                model.add_reactions([reverse_reaction.copy()]) # add reverse reaction to the model
                model.reactions.get_by_id(reverse_reaction.id).reaction = backwards_reaction
                new_reverse_reactions = new_reverse_reactions+[reverse_reaction.id] #to concatenate str to list, we need to make the str a list
                reactions_in_loop = reactions_in_loop+[reverse_reaction.id]
            else:
                reactions_in_loop = reactions_in_loop+[reverse_reaction[0].id] #here, RR would be a list so we need to extract it this way
           
        associated_metabolites = current_reaction.reactants+current_reaction.products

        for x in associated_metabolites:
            ref_id = ref_id +1
            ref_rxn = re.sub(str(pattern_replace), str(ref_id), r)
            metabolites_in_model = [x.id for x in model.metabolites] 
            if not x.id+'_sink' in metabolites_in_model and not re.findall('_sink', x.id) and not x.id+'_sink' in associated_sink_metabolites: # if the sink version of the metabolite in not in the model and the metabolite is not a sink metabolite
                if x.id in DeadEnd_S:
                    sink_reaction = cobra.Reaction(ref_rxn)
                    sink_reaction.id = ref_rxn     
                    sink_reaction.name = "Sink "+x.name
                    sink_reaction.subsystem = 'Putative Sink Reaction'
                    sink_reaction.lower_bound = 0
                    sink_reaction.upper_bound = 0.0
                    #If metabolite is a dead end substrate, it is a product in the sink reaction
                    sink_reaction.add_metabolites({getattr(model.metabolites, x.id):1.0})
                    sink_reaction.notes = {'Confidence Level': '4', 'AUTHORS': 'Sink transport reaction automatically generated'}
                    sink_reaction.gene_reaction_rule = ''
                    model.add_reactions([sink_reaction.copy()])
                    sink_reactions = sink_reactions + [sink_reaction.id]

                elif x.id in DeadEnd_P:
                    sink_reaction = cobra.Reaction(ref_rxn)
                    sink_reaction.id = ref_rxn     
                    sink_reaction.name = "Sink "+x.name
                    sink_reaction.subsystem = 'Putative Sink Reaction'
                    sink_reaction.lower_bound = 0
                    sink_reaction.upper_bound = 0.0
                    #If metabolite is a dead end product, it is a substrate in the sink reaction
                    sink_reaction.add_metabolites({getattr(model.metabolites, x.id):-1.0})           
                    sink_reaction.notes = {'Confidence Level': '4', 'AUTHORS': 'Sink transport reaction automatically generated'}
                    sink_reaction.gene_reaction_rule = ''
                    model.add_reactions([sink_reaction.copy()])
                    sink_reactions = sink_reactions + [sink_reaction.id]
                
                elif x.id in DeadEnd_M: #if the metabolite is in a futile cycle, treat it as a dead-end substrate
                    sink_reaction = cobra.Reaction(ref_rxn)
                    sink_reaction.id = ref_rxn     
                    sink_reaction.name = "Sink "+x.name
                    sink_reaction.subsystem = 'Putative Sink Reaction'
                    sink_reaction.lower_bound = 0
                    sink_reaction.upper_bound = 0.0
                    #If metabolite is a dead end substrate, it is a product in the sink reaction
                    sink_reaction.add_metabolites({getattr(model.metabolites, x.id):1.0})           
                    sink_reaction.notes = {'Confidence Level': '4', 'AUTHORS': 'Sink transport reaction automatically generated'}
                    sink_reaction.gene_reaction_rule = ''
                    model.add_reactions([sink_reaction.copy()])
                    sink_reactions = sink_reactions + [sink_reaction.id]                
                    
    model_sink = model.copy() # model with only sink reactions added
    model_transport = cobra.Model("model_transport") # model to store transport reactions

    # Add transport reactions

    print('Adding transport reactions')
    new_transport_reactions = []
    metabolites_in_1_compound_reaction = list([x.metabolites for x in model.reactions if len(x.metabolites) == 2])

    #For DeadEnd_M (futile cycles):
    for metabolite_pair in allowed_futilecycles:
        # Create a reaction transporting the nth metabolite between compartments
        ref_id = ref_id +1
        ref_rxn = re.sub(str(pattern_replace), str(ref_id), r)
        new_transport_reaction = copy.deepcopy(substrate_transport_metabolite_to(metabolite_pair, ref_rxn, DeadEnd_M))

        # Check if the reaction is already in the model

        for y in range(0,len(metabolites_in_1_compound_reaction)):
            new_metabolites = [x.id for x in list(new_transport_reaction.metabolites)]
            yth_metabolites = [x.id for x in list(metabolites_in_1_compound_reaction[y])]
            overlapping_test = jaccard(new_metabolites, yth_metabolites)
            if overlapping_test == 1: # If a coincidence is found stop the loop
                break
        if overlapping_test != 1:   
            model.add_reactions([new_transport_reaction]) # if no coincidences have found then add the reaction to the model
            model_transport.add_reactions([new_transport_reaction]) # if no coincidences have found then add the reaction to the model           
            new_transport_reactions = new_transport_reactions + [new_transport_reaction]
    new_transport_reactions = list(set(new_transport_reactions))

    del(model)
    return model_sink, model_transport, new_transport_reactions, DeadEnd_S, DeadEnd_P, DeadEnd_M

# Identify minimum number of transport reactions to reduce the number of dead-end metabolites as much as possible


def process_reaction(x, model_sink, model_transport, tolerance):
    """
    This function performs a gap-filling calculation on a metabolic model to identify and add reactions to fill gaps in the metabolic network.
    
    Inputs
    ----------
    x (str): Sink reaction id associated to transported deadend metabolites
    model_sink (Model): The metabolic model to be gap-filled.
    model_transport (Model): The model of transport reactions.
    tolerance (float): The tolerance for the gap-filling calculation.

    Output
    -------
    add_new_reactions (list): List of new reactions added to the model.
    """

    
    try:
            
        # Create a copy of the model
        model_sink_x = model_sink.copy()
        model_sink_x_metabolites = model_sink_x.metabolites
        
        # Initialize list to store new reactions
        add_new_reactions = []
        # Force the xth sink reaction to be active
        model_sink_x.reactions.get_by_id(x).upper_bound = 1000
        model_sink_x.reactions.get_by_id(x).lower_bound = tolerance * 10

        # Perform gap-filling calculation
        solution = gapfill(model_sink_x, model_transport, iterations=1, lower_bound=tolerance, penalties=None, exchange_reactions=False, demand_reactions=False)

        # If no solution is found, check if an extracellular metabolite exists in the model
        if not solution[0]:
            # Extract the extracellular metabolite name associated with the xth sink reaction
            xth_extracel_metabolite = re.sub("[a-z]_sink","e", str(model_sink_x.reactions.get_by_id(x).metabolites).split(' ')[1])
            
            # If extracellular metabolite exists, add exchange reaction
            if xth_extracel_metabolite in model_sink_x_metabolites:
                xth_metabolite = model_sink_x_metabolites.get_by_id(xth_extracel_metabolite)
                
                # Find the associated exchange reaction for the extracellular metabolite
                ex_rxn_in_xth_extracel_metabolite = [x for x in xth_metabolite.reactions if len(x.metabolites) == 1]
                
                # If the extracellular metabolite has an associated exchange reaction, change the boundaries
                if ex_rxn_in_xth_extracel_metabolite[0]:
                    ex_rxn_in_xth_extracel_metabolite[0].upper_bound = 1000
                    ex_rxn_in_xth_extracel_metabolite[0].lower_bound = -1000
                else:
                    # Add a transport reaction associated to the extracellular metabolite
                    ref_rxn = 'EX_'+xth_extracel_metabolite
                    ex_reaction = cobra.Reaction(ref_rxn)
                    ex_reaction.id = ref_rxn     
                    ex_reaction.name = "Exchange " + xth_metabolite.name
                    ex_reaction.subsystem = 'Putative Exchange Reaction'
                    ex_reaction.upper_bound = 1000.0
                    ex_reaction.lower_bound = -1000.0
                    ex_reaction.add_metabolites({getattr(model_sink_x_metabolites, xth_metabolite.id):-1.0})
                    ex_reaction.notes = {'Confidence Level': '4', 'AUTHORS': 'Exchange reaction automatically generated'}
                    ex_reaction.gene_reaction_rule = ''
                    model_sink_x.add_reactions([ex_reaction.copy()])
                    del(ex_reaction)
                
                del(ex_rxn_in_xth_extracel_metabolite)
                del(xth_metabolite)
                # Perform gap-filling calculation again
                solution = gapfill(model_sink_x, model_transport, iterations=1, lower_bound=tolerance, penalties=None, exchange_reactions=False, demand_reactions=False)
            del(xth_extracel_metabolite)

        # If no solution is found, change the direction of the reaction and perform gap-filling calculation again
        if not solution[0]:
            # Change the substrates by the products and viceversa
            print("Reaction: ", model_sink_x.reactions.get_by_id(x))
            print("Before change: ", model_sink_x.reactions.get_by_id(x).reactants, model_sink_x.reactions.get_by_id(x).products)

            for metabolite, coefficient in model_sink_x.reactions.get_by_id(x).metabolites.items():
                coefficient=-coefficient
                model_sink_x.reactions.get_by_id(x).add_metabolites({metabolite: coefficient}, combine=False)
            del(metabolite)
            del(coefficient)
            print("After change: ", model_sink_x.reactions.get_by_id(x).reactants, model_sink_x.reactions.get_by_id(x).products)

            # Force the xth sink reaction to be active
            model_sink_x.reactions.get_by_id(x).upper_bound = 1000
            model_sink_x.reactions.get_by_id(x).lower_bound = tolerance * 10
            # Perform gap-filling calculation
            solution = gapfill(model_sink_x, model_transport, iterations=1, lower_bound=tolerance, penalties=None, exchange_reactions=False, demand_reactions=False)

            if not solution[0]:
                # Extract the extracellular metabolite name associated with the xth sink reaction
                xth_extracel_metabolite = re.sub("[a-z]_sink","e", str(model_sink_x.reactions.get_by_id(x).metabolites).split(' ')[1])

                # If extracellular metabolite exists, add exchange reaction
                if xth_extracel_metabolite in model_sink_x_metabolites: # if extracellular specie exists add excahnge reaction
                    xth_metabolite = model_sink_x_metabolites.get_by_id(xth_extracel_metabolite)
                    ex_rxn_in_xth_extracel_metabolite = [x for x in xth_metabolite.reactions if len(x.metabolites) == 1]                
                    if ex_rxn_in_xth_extracel_metabolite: # if the extracellular specia has an associated Ex reaction change the boundaries
                        ex_rxn_in_xth_extracel_metabolite[0].upper_bound = 1000
                        ex_rxn_in_xth_extracel_metabolite[0].lower_bound = -1000
                    else: # add a transport reaction associated to the extracellular metabolite
                        ref_rxn = 'EX_'+xth_extracel_metabolite
                        ex_reaction = cobra.Reaction(ref_rxn)
                        ex_reaction.id = ref_rxn     
                        ex_reaction.name = "Exchange "+xth_metabolite.name
                        ex_reaction.subsystem = 'Putative Exchange Reaction'
                        ex_reaction.upper_bound = 1000.0
                        ex_reaction.lower_bound = -1000.0
                        ex_reaction.add_metabolites({getattr(model_sink_x_metabolites, xth_metabolite.id):-1.0})
                        ex_reaction.notes = {'Confidence Level': '4', 'AUTHORS': 'Exchange reaction automatically generated'}
                        ex_reaction.gene_reaction_rule = ''
                        model_sink_x.add_reactions([ex_reaction.copy()])
                        del(ex_reaction)
                    del(ex_rxn_in_xth_extracel_metabolite)
                    del(xth_metabolite)
                    solution = gapfill(model_sink_x, model_transport, iterations=1, lower_bound=tolerance, penalties=None, exchange_reactions=False, demand_reactions=False)

                del(xth_extracel_metabolite)
                

        # Save solution (if any)
        if solution[0]:
            print (x, solution, "Solution found: ", model_sink_x.reactions.get_by_id(x))
            print("Sink Reaction: ", model_sink_x.reactions.get_by_id(x))
            del(model_sink_x)
            if solution[0] not in add_new_reactions: #If TR is not already in the list, add it
                add_new_reactions.extend(solution[0])
            del(solution)
            del(x)

         # Delete the model copy to free RAM
    except Exception as e:
        print('\naborted', x,'[[]]')
        print('Exception:', str(e))
        print("Reaction: ", model_sink_x.reactions.get_by_id(x))
         # Delete the model copies to free RAM
        del(model_sink_x)
        del(model_transport)
        del(model_sink)
        del(model_sink_x_metabolites)
        del(x)

    return add_new_reactions


def identify_transport_reaction_to_dead_end_metabolite(model_sink, model_transport):
    """
    Runs a paralellization of the gap-filling calculation from the previous function (process_reaction)

    Inputs 
    ----------
        model_sink (object): The model containing sink reactions.
        model_transport (object): The model containing transport reactions.

    Output
    -------
        list: A list of new reactions associated with the transport of dead-end metabolites.
    """

    #Define list of metabolites transported in the model_transport and sink reactions associated with the transported dead-end metabolites (generates a list of IDs for reactions that have "Sink" in their name and involve at least one metabolite that is also in model_transport.metabolites)
    sink_reactions_associated_to_transported_deadend_metabolites = [x.id for x in model_sink.reactions if "Sink" in x.name and any(y.id in [z.id for z in model_transport.metabolites] for y in x.metabolites)]

    print("Number of sink reactions associated with the transported dead-end metabolites: ", len(sink_reactions_associated_to_transported_deadend_metabolites))
    print("Set of sink reactions associated with the transported dead-end metabolites: ", len(set(sink_reactions_associated_to_transported_deadend_metabolites)))

    # Set the tolerance value
    tolerance = 0.001

    # Disable traceback
    sys.tracebacklimit = 0

    # Initialize lists to store all new reactions and processed subsets
    all_new_reactions = []
    processed_subsets = []

    # Define file paths for saving progress
    reactions_file = os.path.join(project_root, 'files', "processed_reactions.pkl")
    subsets_file = os.path.join(project_root, 'files', "processed_subsets.pkl")

    # Check if the files already exist to resume progress
    if os.path.exists(reactions_file):
        # Load existing reactions if available
        with open(reactions_file, "rb") as f:
            all_new_reactions = pickle.load(f)
        print(f"Loaded {len(all_new_reactions)} reactions from {reactions_file}")

    if os.path.exists(subsets_file):
        # Load processed subsets if available
        with open(subsets_file, "rb") as f:
            processed_subsets = pickle.load(f)
        print(f"Loaded {len(processed_subsets)} processed subsets from {subsets_file}")


    #Create a subset of sink reactions to parallelize the process by parts
    
    # Convert the set to a list to enable indexing and slicing
    reactions_list = list(set(sink_reactions_associated_to_transported_deadend_metabolites))
    #sort the list
    reactions_list.sort()

    # Define the size of each subset
    subset_size = 200

    # Calculate the total number of subsets needed
    num_subsets = (len(reactions_list) + subset_size - 1) // subset_size  # Ceiling division

    # Initialize a list to store subsets
    subsets = []

    # Loop to create subsets
    for i in range(num_subsets):
        start_index = i * subset_size
        end_index = (i + 1) * subset_size
        subset = reactions_list[start_index:end_index]
        subsets.append(subset)

    # Iterate over each subset in subsets
    for subset_index, subset in enumerate(subsets):
        # Skip the subset if it has already been processed
        if subset_index in processed_subsets:
            print(f"Subset {subset_index + 1} already processed. Skipping...")
            continue

        print(f"Processing Subset {subset_index + 1}...")
        
        # Process each sink reaction in parallel using ProcessPoolExecutor
        with concurrent.futures.ProcessPoolExecutor(max_workers=6) as executor: #Check workers depending on the number of cores of the current computer
            futures = [executor.submit(process_reaction, x, model_sink, model_transport, tolerance) for x in subset]
            for future in concurrent.futures.as_completed(futures):
                all_new_reactions.extend(future.result())

        #Save the progress to pickle files
        
        # Save the updated all_new_reactions list to the reactions file
        with open(reactions_file, "wb") as f:
            pickle.dump(all_new_reactions, f)
        
        print("TR printed: ", all_new_reactions)

        # Mark the current subset as processed and save the updated list to the subsets file
        processed_subsets.append(subset_index)
        with open(subsets_file, "wb") as f:
            pickle.dump(processed_subsets, f)

        print("Subset printed: ", processed_subsets)

        print(f"Subset {subset_index + 1} processed and saved to {reactions_file} and {subsets_file}.\n")

    # Clean up the models
    del(model_transport)
    del(model_sink)

    # Return the list of new reactions
    return all_new_reactions


"""Stoichiometric consistency tests for an instance of ``cobra.Model``."""

@annotate(title="Stoichiometric Consistency", format_type="percent")
def test_stoichiometric_consistency(model):
    """
    Expect that the stoichiometry is consistent.

    Stoichiometric inconsistency violates universal constraints:
    1. Molecular masses are always positive, and
    2. On each side of a reaction the mass is conserved.
    A single incorrectly defined reaction can lead to stoichiometric
    inconsistency in the model, and consequently to unconserved metabolites.
    Similar to insufficient constraints, this may give rise to cycles which
    either produce mass from nothing or consume mass from the model.

    Implementation:
    This test uses an implementation of the algorithm presented in
    section 3.1 by Gevorgyan, A., M. G Poolman, and D. A Fell.
    "Detection of Stoichiometric Inconsistencies in Biomolecular Models."
    Bioinformatics 24, no. 19 (2008): 2245.
    doi: 10.1093/bioinformatics/btn425

    """
    ann = test_stoichiometric_consistency.annotation
    is_consistent = consistency.check_stoichiometric_consistency(model)
    ann["data"] = is_consistent
    ann["metric"] = 1.0 - float(is_consistent)
    ann["message"] = wrapper.fill(
        """This model's stoichiometry {}""".format(
            "consistent" if is_consistent else "inconsistent"
        )
    )
    assert is_consistent, ann["message"]


'''
Encapsulate the whole code into a function: Will call it putative_tr
'''
def putative_tr(gem, compartments, sheetname, bm, factor):

    '''
    Runs the code to reduce the number of dead-end metabolites in the model by adding transport reactions 

    Inputs 
    ----------
        gem (str): Name of the GEM SBML file to use as the model.
        compartments (str): Name of the XLSX file with the list of compartments.
        sheetname (str): Name of the sheet inside the XLSX file that defines the compartments for the model used.
        bm (str): The ID for the biomass in the current model.
        factor (float): The factor defining the minimum allowed for biomass (e.g. 0.5 will define that the minimum biomass allowed is 1/2 of the optimal).

    '''
    freeze_support()
    starttimedev = time.time()

    # # open processed_subsets.pkl and print the content (debug)
    # with open("processed_subsets.pkl", "rb") as f:
    #     processed_subsets = pickle.load(f)
    # print("Processed subsets: ", processed_subsets)

    
    model_original = cobra.io.read_sbml_model(gem)
    print("Model solved using: ", model_original.solver) #debug

    model_original = modify_reaction_bounds(model_original)
    #First round, consider DeadEnd_S and DeadEnd_P only
    model_sink, model_transport, new_transport_reactions, DeadEnd_S, DeadEnd_P, DeadEnd_M = add_transport_and_sink_reactions(model_original, compartments, sheetname)

    # Optimize BM as reference

    model_sink.objective = bm
    model_sink.reactions.get_by_id(bm).upper_bound = 1000
    optimal_bm_solution_sink = model_sink.optimize()
    optimal_bm_sink = optimal_bm_solution_sink.objective_value
    model_sink.reactions.get_by_id(bm).lower_bound = optimal_bm_sink*factor # minimum allowed for BM # Define minimum BM (1/2 of optimal)

    # Identify transport reactions to reduce the number of dead-end metabolites

    add_new_reactions = identify_transport_reaction_to_dead_end_metabolite (model_sink, model_transport)

    # Add new transport reactions to the original model

    add_new_reactions_unique = set(add_new_reactions)
    model_transport_final = cobra.Model("model_transport_final") #model to store transport reactions relevant to the original model

    print("Number of new reactions that may be tested for stoichiometry: ", len(add_new_reactions_unique))
    final_reactions_added=[] #list to store the reactions that will be added to the model after testing

    print("Reactions to be tested for stoichiometry: ", add_new_reactions_unique)

    #Uncomment this after the MEMOTE issue is corrected, then remove the default value of 'n' for check_stochiometry
    #check_stochiometry=input("Do you want to check the stoichiometry before adding the new reactions? (y/n): ")
    check_stochiometry='n'
    #ask the user for input until the input is y or n:

    while check_stochiometry not in ['y', 'n']:
        print("Invalid input. Please enter 'y' or 'n'.")
        check_stochiometry=input("Do you want to check the stoichiometry before adding the new reactions? (y/n): ")


    if check_stochiometry == 'y': #check stoichiometry before adding the new reactions
        for x in add_new_reactions_unique:

            model_original.add_reactions([x.copy()])
            is_consistent = consistency.check_stoichiometric_consistency(model_original)  
            if is_consistent:
                print("The stoichiometry of the reaction is consistent. The following reaction will be added: ", x.id)
                model_original.reactions.get_by_id(x.id).lower_bound = -1000
                model_original.reactions.get_by_id(x.id).upper_bound = 1000
                model_transport_final.add_reactions([x.copy()])
                final_reactions_added.append(x.id)
            else:
                print("The stoichiometry of the reaction is not consistent. The following reaction will not be added: ", x.id)
                model_original.remove_reactions([x.id])

    elif check_stochiometry == 'n': #add the new reactions without checking the stoichiometry
        for x in add_new_reactions_unique:

            model_original.add_reactions([x.copy()])

            #make these TR reversible:
            model_original.reactions.get_by_id(x.id).lower_bound = -1000
            model_original.reactions.get_by_id(x.id).upper_bound = 1000

            #create also model with only TR
            model_transport_final.add_reactions([x.copy()])

            print("The following reaction will be added: ", x.id)
            final_reactions_added.append(x.id)


    print("Number of new reactions added after testing: ", len(final_reactions_added))

    # Find number of dead-end metabolites after adding transport reactions
    DeadEnd_M_left, DeadEnd_S_left, DeadEnd_P_left=ident_final_dead_end_metabolites(model_original)
    combined_dead_end_metabolites = DeadEnd_S_left + DeadEnd_P_left + DeadEnd_M_left
    #Original dead-end metabolites
    print("Dead-end substrates initially: ", len(DeadEnd_S))
    print("Dead-end products initially: ", len(DeadEnd_P))
    print("Dead-end metabolites initially: ", len(DeadEnd_M))

    #Dead-end metabolites after adding transport reactions
    print("Number of dead-end metabolites left: ", len(combined_dead_end_metabolites))
    print("Dead-end substrates left: ", len(DeadEnd_S_left))
    print("Dead-end products left: ", len(DeadEnd_P_left))
    print("Dead-end metabolites left: ", len(DeadEnd_M_left))



    print("\n\nFirst round of TR completed")


    #Only if there are dead-end metabolites left, we will consider them in the second round. Otherwise, we will stop here.
    if len(DeadEnd_M_left)!=0:
        #Second round, consider DeadEnd_M only

        #We use a different function to add transport and sink reactions to the model considering only the dead-end metabolites
        model_sink_2, model_transport_2, new_transport_reactions_2, DeadEnd_S_2, DeadEnd_P_2, DeadEnd_M_2 = add_transport_and_sink_reactions_futile_cycles(model_original, compartments, sheetname)

        # Optimize BM as reference

        model_sink_2.objective = bm
        model_sink_2.reactions.get_by_id(bm).upper_bound = 1000
        optimal_bm_solution_sink = model_sink_2.optimize()
        optimal_bm_sink = optimal_bm_solution_sink.objective_value
        model_sink_2.reactions.get_by_id(bm).lower_bound = optimal_bm_sink*factor # minimum allowed for BM # Define minimum BM (1/2 of optimal)
        
        # Identify transport reactions to reduce the number of dead-end metabolites

        add_new_reactions_2 = identify_transport_reaction_to_dead_end_metabolite (model_sink_2, model_transport_2)

        # Add new transport reactions to the original model

        add_new_reactions_unique_2 = set(add_new_reactions_2)

        print("Number of new reactions that may be tested for stoichiometry: ", len(add_new_reactions_unique_2))
        final_reactions_added_2=[] #list to store the reactions that will be added to the model after testing

        #Uncomment this after the MEMOTE issue is corrected, then remove the default value of 'n' for check_stochiometry
        #check_stochiometry=input("Do you want to check the stoichiometry before adding the new reactions? (y/n): ")
        check_stochiometry='n'
        #ask the user for input until the input is y or n:

        while check_stochiometry not in ['y', 'n']:
            print("Invalid input. Please enter 'y' or 'n'.")
            check_stochiometry=input("Do you want to check the stoichiometry before adding the new reactions? (y/n): ")


        if check_stochiometry == 'y': #check stoichiometry before adding the new reactions
            for x in add_new_reactions_unique_2:

                model_original.add_reactions([x.copy()])
                is_consistent = consistency.check_stoichiometric_consistency(model_original)  
                if is_consistent:
                    print("The stoichiometry of the reaction is consistent. The following reaction will be added: ", x.id)
                    model_original.reactions.get_by_id(x.id).lower_bound = -1000
                    model_original.reactions.get_by_id(x.id).upper_bound = 1000
                    model_transport_final.add_reactions([x.copy()])
                    final_reactions_added_2.append(x.id)
                else:
                    print("The stoichiometry of the reaction is not consistent. The following reaction will not be added: ", x.id)
                    model_original.remove_reactions([x.id])

        elif check_stochiometry == 'n': #add the new reactions without checking the stoichiometry
            for x in add_new_reactions_unique_2:

                model_original.add_reactions([x.copy()])

                #make these TR reversible:
                model_original.reactions.get_by_id(x.id).lower_bound = -1000
                model_original.reactions.get_by_id(x.id).upper_bound = 1000

                #create also model with only TR
                model_transport_final.add_reactions([x.copy()])

                print("The following reaction will be added: ", x.id)
                final_reactions_added_2.append(x.id)


        print("Number of new reactions added after testing: ", len(final_reactions_added_2))

        # Find number of dead-end metabolites after adding transport reactions
        DeadEnd_M_left_2, DeadEnd_S_left_2, DeadEnd_P_left_2=ident_final_dead_end_metabolites(model_original)
        combined_dead_end_metabolites_2 = DeadEnd_S_left_2 + DeadEnd_P_left_2 + DeadEnd_M_left_2
        print("Number of dead-end metabolites left: ", len(combined_dead_end_metabolites_2))
        print("Dead-end substrates left: ", len(DeadEnd_S_left_2))
        print("Dead-end products left: ", len(DeadEnd_P_left_2))
        print("Dead-end metabolites left: ", len(DeadEnd_M_left_2))   
        
        #Total reactions added after testing for stoichiometry
        total_reactions_added = final_reactions_added + final_reactions_added_2
        #Total reactions added before testing for stoichiometry
        add_new_reactions_total = add_new_reactions + add_new_reactions_2

    # Save the model

        '''Save security copy in files'''

        cobra.io.write_sbml_model(model_original, os.path.join(project_root, 'models', 'model_full_THG.xml')) #save the model with all reactions

        cobra.io.write_sbml_model(model_transport_final, os.path.join(project_root, 'models', 'model_transport_final_THG.xml')) #save the model with only transport reactions

        with open(os.path.join(project_root, 'files', 'new_transport_reactions_untested_THG.pkl'), 'wb') as file:
            pickle.dump(add_new_reactions_total, file) #save the list of new transport reactions still to be tested for stoichiometry

        with open(os.path.join(project_root, 'files', 'final_reactions_added_THG.pkl'), 'wb') as file:
            pickle.dump(total_reactions_added, file) #save the list of reactions that will be added to the model after testing for stoichiometry

        with open(os.path.join(project_root, 'files', 'DeadEnd_S_THG.pkl'), 'wb') as file:
            pickle.dump(DeadEnd_S, file) #save the list of dead-end substrates found before adding the transport reactions

        with open(os.path.join(project_root, 'files', 'DeadEnd_P_THG.pkl'), 'wb') as file:
            pickle.dump(DeadEnd_P, file) #save the list of dead-end products found before adding the transport reactions

        with open(os.path.join(project_root, 'files', 'futilecycle_mets_THG.pkl'), 'wb') as file:
            pickle.dump(DeadEnd_M_2, file) #save the list of dead-end metabolites in isolated futile cycles found before adding the transport reactions

        with open(os.path.join(project_root, 'files', 'Total_DeadEnd_M_left.pkl'), 'wb') as file:
            pickle.dump(combined_dead_end_metabolites_2, file) #save the list of total dead-end metabolites left after adding the transport reactions

        #If code runs without errors, delete the files with the security copies of the subsets and reactions from the parallelization of the gap-filling calculation
        os.remove(os.path.join(project_root, 'files', "processed_reactions.pkl"))
        os.remove(os.path.join(project_root, 'files', "processed_subsets.pkl"))

    else: #if there are no dead-end metabolites left after the first round

        #Total reactions added after testing for stoichiometry
        total_reactions_added = final_reactions_added 
        #Total reactions added before testing for stoichiometry
        add_new_reactions_total = add_new_reactions 
    
        # Save the model

        '''Save security copy in files'''

        cobra.io.write_sbml_model(model_original, os.path.join(project_root, 'models', 'model_full_THG_try2.xml')) #save the model with all reactions

        cobra.io.write_sbml_model(model_transport_final, os.path.join(project_root, 'models', 'model_transport_final_THG.xml')) #save the model with only transport reactions

        with open(os.path.join(project_root, 'files','new_transport_reactions_untested_THG.pkl'), 'wb') as file:
            pickle.dump(add_new_reactions_total, file) #save the list of new transport reactions still to be tested for stoichiometry

        with open(os.path.join(project_root, 'files','final_reactions_added_THG.pkl'), 'wb') as file:
            pickle.dump(total_reactions_added, file) #save the list of reactions that will be added to the model after testing for stoichiometry

        with open(os.path.join(project_root, 'files','DeadEnd_S_THG.pkl'), 'wb') as file:
            pickle.dump(DeadEnd_S, file) #save the list of dead-end substrates found before adding the transport reactions

        with open(os.path.join(project_root, 'files','DeadEnd_P_THG.pkl'), 'wb') as file:
            pickle.dump(DeadEnd_P, file) #save the list of dead-end products found before adding the transport reactions

        with open(os.path.join(project_root, 'files','futilecycle_mets_THG.pkl'), 'wb') as file:
            pickle.dump(DeadEnd_M, file) #save the list of dead-end metabolites in isolated futile cycles found before adding the transport reactions

        with open(os.path.join(project_root, 'files','Total_DeadEnd_M_left.pkl'), 'wb') as file:
            pickle.dump(combined_dead_end_metabolites, file) #save the list of total dead-end metabolites left after adding the transport reactions

        #If code runs without errors, delete the files with the security copies of the subsets and reactions from the parallelization of the gap-filling calculation
        os.remove(os.path.join(project_root, 'files', "processed_reactions.pkl"))
        os.remove(os.path.join(project_root, 'files', "processed_subsets.pkl"))
    
    #To measure the total execution time of the program: 
    endtimedev = time.time() 
    print("Time for dev: ", endtimedev - starttimedev) #in seconds

    

  
'''
Now, invoke the function in the main code
'''

if __name__ == '__main__':

    freeze_support()

    gem = os.path.join(project_root, "models", "THG_EC_endoB_exchange.xml")  
    sheetname = "EndoB"

    compartments = os.path.join(project_root, 'files', 'ListOfCompartments_sept2024.xlsx')

    bm = 'MAR13082'
    factor = 0.5
    putative_tr(gem, compartments, sheetname, bm, factor)

