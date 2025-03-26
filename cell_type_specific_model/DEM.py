#load pkl file 'Total_DeadEnd_M_left_round2.pkl'

import pickle
import pdb
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
import os
import multiprocessing
import re
import numpy as np
import pandas as pd
import copy
from ptr_multi_round import jaccard, search_metabolite_locations, identify_connected_species, product_transport_metabolite_to, substrate_transport_metabolite_to, process_reaction, ident_final_dead_end_metabolites

# Determine the current file's directory and the project root.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")
# Add the project root to sys.path to access top-level folders like 'functions' and 'models'
if project_root not in sys.path:
    sys.path.append(project_root)

# # Use add_transport_and_sink_reactions from main code, but simplified because we already have the DEM
def add_transport_and_sink_reactions_simplified (DeadEnd_S, DeadEnd_P, model_original, compartments, sheetname):

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
                    #If metabolite is a dead end substrate, it is a substrate in the sink reaction
                    sink_reaction.add_metabolites({getattr(model.metabolites, x.id):-1.0}) 
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
                    #If metabolite is a dead end product, it is a product in the sink reaction
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

    for x in new_transport_reactions:
        print(x.id, x.reaction)
    del(model)
    return model_sink, model_transport, new_transport_reactions, DeadEnd_S, DeadEnd_P



def identify_transport_reaction_to_dead_end_metabolite(model_sink, model_transport, round_id):
    """
    Process the gap-filling in parallel and store progress in round-specific checkpoint files.
    """
    # Define round-specific file names.
    reactions_file = f"processed_reactions_glycocalix.pkl"
    reactions_file = os.path.join(project_root, 'files', reactions_file)
    subsets_file = f"processed_subsets_glycocalix.pkl"
    subsets_file = os.path.join(project_root, 'files', subsets_file)

    # Try to load previous progress.
    all_new_reactions = []
    processed_subsets = []
    if os.path.exists(reactions_file):
        with open(reactions_file, "rb") as f:
            all_new_reactions = pickle.load(f)
        print(f"Loaded {len(all_new_reactions)} reactions from {reactions_file}")
    if os.path.exists(subsets_file):
        with open(subsets_file, "rb") as f:
            processed_subsets = pickle.load(f)
        print(f"Loaded {len(processed_subsets)} processed subsets from {subsets_file}")

        #Print content of subsets
        print("Content of processed subsets: ", processed_subsets)

        #Print number of reactions in all_new_reactions
        print("Number of reactions in all_new_reactions: ", len(all_new_reactions))

    # Build the list of sink reactions to process.
    sink_reactions_associated = [
        x.id for x in model_sink.reactions 
        if "Sink" in x.name and any(y.id in [z.id for z in model_transport.metabolites] for y in x.metabolites)
    ]
    print("Number of sink reactions:", len(sink_reactions_associated))

    tolerance = 1e-7 #current model tolerance, need to set it as a variable (debug)
    sys.tracebacklimit = 0

    # Prepare a sorted list and divide into subsets.
    reactions_list = sorted(set(sink_reactions_associated))
    subset_size = 200
    num_subsets = (len(reactions_list) + subset_size - 1) // subset_size
    subsets = [reactions_list[i * subset_size:(i + 1) * subset_size] for i in range(num_subsets)]

    # Process each subset.
    for subset_index, subset in enumerate(subsets):
        if subset_index in processed_subsets:
            print(f"Subset {subset_index + 1} already processed. Skipping...")
            continue
        print(f"Processing Subset {subset_index + 1} in round {round_id}...")
        with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
            futures = [executor.submit(process_reaction, x, model_sink, model_transport, tolerance) for x in subset]
            for future in concurrent.futures.as_completed(futures):
                all_new_reactions.extend(future.result())

        # Save progress after processing each subset.
        with open(reactions_file, "wb") as f:
            pickle.dump(all_new_reactions, f)
        processed_subsets.append(subset_index)
        with open(subsets_file, "wb") as f:
            pickle.dump(processed_subsets, f)
        print(f"Subset {subset_index + 1} of round {round_id} processed.")

    # When the round is complete, delete the checkpoint files.
    if os.path.exists(reactions_file):
        os.remove(reactions_file)
    if os.path.exists(subsets_file):
        os.remove(subsets_file)
    
    del(model_transport)
    del(model_sink)
    return all_new_reactions



def main():

    freeze_support()

    with open(os.path.join(project_root, 'files','Total_DeadEnd_M_left_round2_endoB.pkl'), 'rb') as f:
        dead_end_mets = pickle.load(f)

    print("Number of dead-end metabolites left after round 2: ", len(dead_end_mets))

    model = cobra.io.read_sbml_model(os.path.join(project_root, 'models','model_full_THG_round2_endoB.xml'))

    # # See which of the dead-end metabolites are in glycocalix (compartment 'y')

    glyco_dem_mets=[]
    for met in model.metabolites:
        if met.id in dead_end_mets and met.compartment == 'y':
            glyco_dem_mets.append(met)

    # print("Number of dead-end metabolites in glycocalix: ", len(glyco_dem_mets))
    # print("Dead-end metabolites in glycocalix: ", [met.id for met in glyco_dem_mets])

    # #open final_reactions_added_THG_round2.pkl 

    # with open('final_reactions_added_THG_round2.pkl', 'rb') as f:
    #     final_reactions_added_THG_round2 = pickle.load(f)

    # # check reactions in model and the compartments they involve
    # for rxn_id in final_reactions_added_THG_round2:
    #     rxn = model.reactions.get_by_id(rxn_id)
    #     comps = {met.compartment for met in rxn.metabolites}
    #     print(rxn_id, comps)

    # ''' 
    # TR glycocalix, both connect glycocalix to extracellular space only
    # None of the added TR in PTR connects glycocalix to cell membrane
    # MARMAR67714
    # MARMAR90226
    # '''

    # # Find the transport reactions that involve glycocalix
    # glyco_transport_reactions = [rxn for rxn in model.reactions if rxn.id.startswith('MARMAR') and 'y' in {met.compartment for met in rxn.metabolites}]
    # print("Number of transport reactions that involve glycocalix: ", len(glyco_transport_reactions))
    # print("Transport reactions that involve glycocalix: ", [rxn.id for rxn in glyco_transport_reactions])

    # # Display all the reactions that involve glycocalix

    glyco_reactions = [rxn for rxn in model.reactions if 'y' in {met.compartment for met in rxn.metabolites}]
    # glyco_mets = [met for met in model.metabolites if met.compartment == 'y']
    # print("Number of reactions that involve glycocalix: ", len(glyco_reactions))
    # print("Number of metabolites in glycocalix: ", len(glyco_mets))

    # for rxn in glyco_reactions:
    #     print(rxn.id, ":", rxn.reaction)

    # #For dead-end metabolites in glycocalix, create a set of their ids minus the last character (to remove the compartment)
    # met_ids = {met.id[:-1] for met in glyco_dem_mets}


    # # Find metabolites in the model whose id starts by these met_ids
    # eq_mets= []
    # for id in met_ids:
    #     print("For id: ", id, ". Equivalent metabolites: ")
    #     for met in model.metabolites:
    #         if met.id.startswith(id) and (met.compartment == 'a' or met.compartment=='e'): #glycocalix only connects to a and e compartments
    #             print(met.id, met.compartment, sep= " , ")
    #             eq_mets.append(met)
    #     print()

    # # Iterate over the equivalent metabolites and find which ones are dead-ends 
    # eq_dead_end_mets = []
    # for met in eq_mets:
    #     if met.id in dead_end_mets:
    #         eq_dead_end_mets.append(met)

    # # And find which ones are not dead-ends
    # eq_met_not_dem = list(set(eq_mets) - set(eq_dead_end_mets))
    # print("Number of equivalent metabolites that are dead-ends: ", len(eq_dead_end_mets))
    # print("Number of equivalent metabolites that are not dead-ends: ", len(eq_met_not_dem))


    # Now, create transport reactions between the equivalent metabolites that are not dead-ends and glycocalix dead-end metabolites

    # Figure out how to create TR
    # Create a sink reaction as well
    # Once these are done, run gapfilling with these TR reactions as the universal model. See if the gapfilling is successful.

    # Divde dem into substrates and products
    DeadEnd_S = set()
    DeadEnd_P = set()

    for met in glyco_dem_mets:
        print("Metabolite: ", met.id)
        if len(met.reactions)==1:
            rxn = list(met.reactions)[0] 
            if met in rxn.reactants:
                DeadEnd_S.add(met.id)
            elif met in rxn.products:
                DeadEnd_P.add(met.id)
        elif len(met.reactions)>1:
            print("Metabolite is in more than one reaction. Check manually.")

    print("Number of dead-end substrates: ", len(DeadEnd_S))
    print("Number of dead-end products: ", len(DeadEnd_P))

    DeadEnd_P= list(DeadEnd_P)
    DeadEnd_S= list(DeadEnd_S)


    compartments = os.path.join(project_root, 'files', 'ListOfCompartments_sept2024.xlsx')
    sheetname = "EndoB"

    bm = 'MAR13082'
    factor = 0.5

    model_original = model.copy()
    #First round, consider DeadEnd_S and DeadEnd_P only
    model_sink, model_transport, new_transport_reactions, DeadEnd_S, DeadEnd_P = add_transport_and_sink_reactions_simplified(DeadEnd_S, DeadEnd_P, model_original, compartments, sheetname)

    # # Optimize BM as reference

    model_sink.objective = bm
    model_sink.reactions.get_by_id(bm).upper_bound = 1000
    optimal_bm_solution_sink = model_sink.optimize()
    optimal_bm_sink = optimal_bm_solution_sink.objective_value
    model_sink.reactions.get_by_id(bm).lower_bound = optimal_bm_sink*factor # minimum allowed for BM # Define minimum BM (1/2 of optimal)
    
    model_transport_final = cobra.Model("model_transport_final") #model to store transport reactions relevant to the original model
    final_reactions_added = []
    
    # # Identify transport reactions to reduce the number of dead-end metabolites

    # add_new_reactions = identify_transport_reaction_to_dead_end_metabolite (model_sink, model_transport, round_id=1)

    # # Add new transport reactions to the original model

    # add_new_reactions_unique = set(add_new_reactions)
    # print("Number of new reactions that may be tested for stoichiometry: ", len(add_new_reactions_unique))

    # for x in add_new_reactions_unique:
            
    #             print("Testing reaction: ", x.id)
    #             model_original.add_reactions([x.copy()])
    #             sol = model_original.optimize()
    #             print("Objective value: ", sol.objective_value, "Status: ", sol.status)

    #             #make these TR reversible:
    #             model_original.reactions.get_by_id(x.id).lower_bound = -1000
    #             model_original.reactions.get_by_id(x.id).upper_bound = 1000

    #             if sol.objective_value > 0: #if the reaction is feasible, add it to the model
    #                 print("The following reaction will be added: ", x.id)
    #                 final_reactions_added.append(x.id)
    #                 model_transport_final.add_reactions([x.copy()])
    #             else:
    #                 print("The following reaction will not be added: ", x.id)
    #                 model_original.remove_reactions([x.id])

    # print("Number of new reactions added after testing: ", len(final_reactions_added))


    # For endoB, the add_new_reactions_unique is empty, as gapfilling was not successful for any of them
    # I have commented out the code that would try to add these reactions to the model (as it would not work)
    # We will now try to add manually all the candidate transport reactions (stored in new_transport_reactions)

    for x in new_transport_reactions:
            
                print("Testing reaction: ", x.id)
                model_original.add_reactions([x.copy()])
                sol = model_original.optimize()
                print("Objective value: ", sol.objective_value, "Status: ", sol.status)

                if sol.objective_value > 0: #if the reaction is feasible, add it to the model
                    print("The following reaction will be added: ", x.id)
                    #make these TR reversible:
                    model_original.reactions.get_by_id(x.id).lower_bound = -1000
                    model_original.reactions.get_by_id(x.id).upper_bound = 1000
                    final_reactions_added.append(x.id)
                    model_transport_final.add_reactions([x.copy()])
                else:
                    print("The following reaction will not be added: ", x.id)
                    model_original.remove_reactions([x.id])

    print("Number of new reactions added after testing: ", len(final_reactions_added))

    # Find number of dead-end metabolites after adding transport reactions
    _, DeadEnd_S_left, DeadEnd_P_left=ident_final_dead_end_metabolites(model_original)
    combined_dead_end_metabolites = DeadEnd_S_left + DeadEnd_P_left
    #Original dead-end metabolites
    print("Dead-end substrates initially: ", len(DeadEnd_S))
    print("Dead-end products initially: ", len(DeadEnd_P))

    #Dead-end metabolites after adding transport reactions
    print("Number of dead-end metabolites left: ", len(combined_dead_end_metabolites))
    print("Dead-end substrates left: ", len(DeadEnd_S_left))
    print("Dead-end products left: ", len(DeadEnd_P_left))

    print("\n\nFirst round of TR completed")

    # Then save the model and other files:
    cobra.io.write_sbml_model(model_original, os.path.join(project_root, 'models','model_after_round1_glycocalix.xml'))  # Save the model with all reactions.
    cobra.io.write_sbml_model(model_transport_final, os.path.join(project_root, 'models','model_transport_after_round1_glycocalix.xml')) # Save only the transport reactions.
    with open(os.path.join(project_root, 'files','final_reactions_added_round1_glycocalix.pkl'), 'wb') as file:
        pickle.dump(final_reactions_added, file)
    with open(os.path.join(project_root, 'files','DeadEnd_S_glycocalix.pkl'), 'wb') as file: # Save the dead-end metabolites found before the first round of TR.
        pickle.dump(DeadEnd_S, file)
    with open(os.path.join(project_root, 'files','DeadEnd_P_glycocalix.pkl'), 'wb') as file:
        pickle.dump(DeadEnd_P, file)
    with open(os.path.join(project_root, 'files','Total_DeadEnd_M_left_glycocalix.pkl'), 'wb') as file: # Save the dead-end metabolites found after the first round of TR.
        pickle.dump(combined_dead_end_metabolites, file)


    # # # run FVA on the reactions that involve glycocalix
    # # fva = cobra.flux_analysis.flux_variability_analysis(model, reaction_list=[rxn.id for rxn in glyco_reactions])

    # # # Find reactions that are blocked in the glycocalix
    # # blocked_glyco_reactions = {rxn_id for rxn_id, v in fva.items() if v['minimum'] == 0 and v['maximum'] == 0}
    # # print("Number of blocked reactions in glycocalix: ", len(blocked_glyco_reactions))


    # print(new_transport_reactions)

    # # save new_transport_reactions as a pickle file
    # with open('possible_glycocalix_transport_reactions.pkl', 'wb') as f: pickle.dump(new_transport_reactions, f)

    # # save model transport
    # cobra.io.write_sbml_model(model_transport, "model_transport_glycocalix.xml")


    # load model_transport
    # model_transport = cobra.io.read_sbml_model("model_transport_glycocalix.xml")

    # # # load possible_glycocalix_transport_reactions.pkl
    # # with open('possible_glycocalix_transport_reactions.pkl', 'rb') as f: 
    # #     possible_transport_reactions = pickle.load(f)

    # for rxn in model_transport.reactions:
    #     rxn.lower_bound = -1000
    #     rxn.upper_bound = 1000

    #     #add rxn to original model
    #     model.add_reactions([rxn])

    # model.objective = bm
    # optimal_bm_solution = model.optimize()
    # print("Optimal BM: ", optimal_bm_solution.objective_value)

    # #save model
    # cobra.io.write_sbml_model(model, "model_full_THG_glycocalix.xml")

    # #find number of blocked reactions in model
    # br= cobra.flux_analysis.find_blocked_reactions(model)

    # #find number of blocked reactions in glycocalix
    # br_glyco=[]

    # glyco_reactions_id = [rxn.id for rxn in glyco_reactions]
    
    # br_glyco = [rxn for rxn in br if rxn in glyco_reactions_id]

    # print("Number of blocked reactions in glycocalix: ", len(br_glyco))


    # # If feasible, see how many DeadEnd_S and DeadEnd_P are left in glycocalix

    # _, DeadEnd_S_left, DeadEnd_P_left = ident_final_dead_end_metabolites(model)

    # # save DeadEnd_S_left and DeadEnd_P_left as pickle files for further analysis
    # with open('DeadEnd_S_left_glycocalix_notonly.pkl', 'wb') as f: pickle.dump(DeadEnd_S_left, f)
    # with open('DeadEnd_P_left_glycocalix_notonly.pkl', 'wb') as f: pickle.dump(DeadEnd_P_left, f)

    # # find which of these are in glycocalix
    # DeadEnd_S_left_glycocalix = []
    # DeadEnd_P_left_glycocalix = []

    # for met in model.metabolites:
    #     if met.id in DeadEnd_S_left and met.compartment == 'y':
    #         DeadEnd_S_left_glycocalix.append(met)

    #     if met.id in DeadEnd_P_left and met.compartment == 'y':
    #         DeadEnd_P_left_glycocalix.append(met)


    # print("Before new TR (in glycocalix): ")
    # print("Number of DeadEnd_S: ", len(DeadEnd_S))
    # print("Number of DeadEnd_P: ", len(DeadEnd_P))
    # print("After new TR (in glycocalix): ")
    # print("Number of DeadEnd_S left in glycocalix: ", len(DeadEnd_S_left_glycocalix))
    # print("Number of DeadEnd_P left in glycocalix: ", len(DeadEnd_P_left_glycocalix))

    pdb.set_trace()



if __name__ == "__main__":
    freeze_support()

    main()