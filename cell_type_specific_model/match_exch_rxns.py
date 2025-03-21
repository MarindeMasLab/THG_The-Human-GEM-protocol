from cobra.io import read_sbml_model, write_sbml_model
import pandas as pd
import pickle
import logging
from copy import deepcopy
import pdb
import gurobipy
import os
import sys

# Determine the current file's directory and the project root.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")
# Add the project root to sys.path to access top-level folders like 'functions' and 'models'
if project_root not in sys.path:
    sys.path.append(project_root)

def match_exch_rxns(model_new, model_base, match_df, tol=1e-9):
    '''
    Match exchange reactions between two models.

    Parameters
    ----------
    model_new : cobra.Model
        The model to modify.
    model_base : cobra.Model
        The model to match to.
    match_df : pandas.DataFrame
        Dataframe with columns 'rxns' and 'rxnRecon3DID'.
    
    Returns
    -------
    cobra.Model
        The modified model.
    '''

    model_modified = model_new.copy()

    objective_reaction = model_modified.objective.expression
    logging.info(f"Objective reaction: {objective_reaction}")
    reaction_ids_new_exchange = [reaction.id for reaction in model_new.boundary]
    logging.info(f"Number of exchange reactions in new model: {len(reaction_ids_new_exchange)}")
    reactions_found = 0
    reactions_not_found = []
    common_rs_list = [] # list of reactions that are found in both models, debug


    # first modify all the reactions that have been found
    for index, row in match_df.iterrows():
        new_name = row['rxns']
        if new_name not in reaction_ids_new_exchange:
            continue
        base_name = str(row['rxnRecon3DID'])
        if base_name != "":
            base_name = base_name.replace('[', '(').replace(']', ')')
            
            try:
                non_modified_new_bounds = deepcopy(model_modified.reactions.get_by_id(new_name).bounds)
                
                # This is the target of the try except code block
                base_bounds = deepcopy(model_base.reactions.get_by_id(base_name).bounds)
                model_modified.reactions.get_by_id(new_name).lower_bound = base_bounds[0]
                model_modified.reactions.get_by_id(new_name).upper_bound = base_bounds[1]
                
                # All of this should not return errors if the code before works
                opt = model_modified.optimize()
                if opt.objective_value == None or opt.objective_value < tol:
                    logging.warning(f'Reaction not feasible:{new_name}')

                    # Relax bounds by increments of 0.1 until feasible
                    increment = 0.1
                    start_bounds = deepcopy(model_modified.reactions.get_by_id(new_name).bounds)
                    while opt.objective_value == None or opt.objective_value < tol:
                        
                        # make new increment bounds
                        increment_lb = round(start_bounds[0] - increment * abs(base_bounds[0]), 7)
                        increment_ub = round(start_bounds[1] + increment * abs(base_bounds[1]), 7)

                        # if the new bounds are too small, set to 0
                        model_modified.reactions.get_by_id(new_name).lower_bound = 0 if abs(increment_lb) < tol else increment_lb
                        model_modified.reactions.get_by_id(new_name).upper_bound = 0 if abs(increment_ub) < tol else increment_ub 
                        opt = model_modified.optimize()

                        # If increment is too large, break
                        if increment >= 10: # 10 is arbitrary, means bounds can increase 10 times
                            bounds_not_found = True
                            break 
                        increment += 0.1
                    else:
                        bounds_not_found = False
                        logging.warning(f'Reaction feasible after increment:{new_name}')
                        logging.info(f'Bounds in the base model:{base_bounds}')
                        logging.info(f'Bounds in the new model:{model_modified.reactions.get_by_id(new_name).bounds}')
                        logging.info(f'Objective value:{opt.objective_value}')
                    # If increment is too large, revert bounds to original
                    if bounds_not_found:
                        logging.warning(f'Reaction not feasible after increment:{new_name}')
                        model_modified.reactions.get_by_id(new_name).lower_bound = non_modified_new_bounds[0]
                        model_modified.reactions.get_by_id(new_name).upper_bound = non_modified_new_bounds[1]
                        opt_after = model_modified.optimize()
                        if opt_after.objective_value == None:
                            logging.error(f'Reaction not feasible after revert:{new_name}')
                else:
                    reactions_found += 1
                    #print(reactions_found)
                    common_rs_list.append(new_name) # debug

            except:
                #print('Reaction not found in old model:', new_name)
                reactions_not_found.append(new_name)
    # after modifying all the reactions that have been found, set the bounds of the remaining reactions to 0
    # and check if the model is still feasible
    # if not, revert the bounds to the original bounds
    for new_name in reactions_not_found:   
        non_modified_new_bounds = deepcopy(model_modified.reactions.get_by_id(new_name).bounds)
        if non_modified_new_bounds[0] == 0 or non_modified_new_bounds[1] == 0:
            logging.error(f'Bound already 0:{new_name}')
        #opt_debug_test = model_modified.optimize()
        #logging.debug(f'Objective value debug:{opt_debug_test.objective_value}')
        model_modified.reactions.get_by_id(new_name).lower_bound = 0
        model_modified.reactions.get_by_id(new_name).upper_bound = 0
        opt_0_bounds = model_modified.optimize()


        if opt_0_bounds.objective_value == None or opt_0_bounds.objective_value < tol:
            logging.warning(f'Model not feasible when set to 0 bounds:{new_name}')
            if non_modified_new_bounds[0] == 0 or non_modified_new_bounds[1] < tol:
                logging.error(f'Saved bounds 0:{new_name}')
            model_modified.reactions.get_by_id(new_name).lower_bound = non_modified_new_bounds[0]
            model_modified.reactions.get_by_id(new_name).upper_bound = non_modified_new_bounds[1]
            opt_0_after = model_modified.optimize()
            logging.debug(f'Objective value after revert:{opt_0_after.objective_value}')
            if opt_0_after.objective_value == None or opt_0_after.objective_value < tol:
                logging.error(f'Model not feasible after 0 revert:{new_name}')
                logging.error(f'Bounds in the new model:{model_modified.reactions.get_by_id(new_name).bounds}')
                logging.error(f'Saved bounds:{non_modified_new_bounds}')

        

    logging.info(f'Reactions found in old model: {reactions_found}')
    logging.info(f'Reactions not found in old model: {len(reactions_not_found)}')

    metab_df = get_metab_df(model_modified, reactions_not_found)

    return model_modified, metab_df , common_rs_list 

def get_metab_df(model, reactions):
    '''
    Get a dataframe of metabolites and the reactions they are involved in.

    Parameters
    ----------
    model : cobra.Model
        The model to use.
    reactions : list
        List of reaction ids.
    '''
    # Calculate the total number of metabolites
    total_metabolites = sum(len(model.reactions.get_by_id(reaction_id).metabolites) for reaction_id in reactions)
    
    # Preassign list lengths
    metab_ids = [None] * total_metabolites
    metab_names = [None] * total_metabolites
    metab_kegg = [None] * total_metabolites
    metab_seed = [None] * total_metabolites
    metab_bigg = [None] * total_metabolites
    metab_chebi = [None] * total_metabolites
    metab_pubchem = [None] * total_metabolites
    metab_vmh = [None] * total_metabolites
    metab_inchi = [None] * total_metabolites
    metab_inchikey = [None] * total_metabolites

    reaction_ids = [None] * total_metabolites
    reaction_names = [None] * total_metabolites
    reaction_ubs = [None] * total_metabolites
    reaction_lbs = [None] * total_metabolites
    reaction_bigg = [None] * total_metabolites


    index = 0
    for reaction_id in reactions:
        reaction = model.reactions.get_by_id(reaction_id)
        metab = list(reaction.metabolites.keys())[0]
        for metab in reaction.metabolites:
            metab_ids[index] = metab.id
            metab_names[index] = metab.name
            metab_kegg[index] = metab.annotation.get('kegg.compound')
            metab_seed[index] = metab.annotation.get('seed.compound')
            metab_bigg[index] = metab.annotation.get('bigg.metabolite')
            metab_chebi[index] = metab.annotation.get('chebi')
            metab_pubchem[index] = metab.annotation.get('pubchem.compound')
            metab_vmh[index] = metab.annotation.get('vmh.metabolite')
            metab_inchi[index] = metab.annotation.get('inchi')
            metab_inchikey[index] = metab.annotation.get('inchikey')

            reaction_ids[index] = reaction_id
            reaction_names[index] = reaction.name
            reaction_ubs[index] = reaction.upper_bound
            reaction_lbs[index] = reaction.lower_bound
            reaction_bigg[index] = reaction.annotation.get('bigg.reaction')
            index += 1

    metab_df = pd.DataFrame({
        'Metabolite IDs': metab_ids, 
        'Metabolite names': metab_names, 
        'Metabolite KeggIDs': metab_kegg, 
        'Metabolite SeedIDs': metab_seed, 
        'Metabolite BiGGIDs': metab_bigg, 
        'Metabolite ChebiIDs': metab_chebi,
        'Metabolite PubchemIDs': metab_pubchem,
        'Metabolite VMHIDs': metab_vmh,
        'Metabolite Inchi': metab_inchi,
        'Metabolite Inchikey': metab_inchikey,

        'Reaction IDs': reaction_ids,
        'Reaction names': reaction_names, 
        'Reaction UBs': reaction_ubs, 
        'Reaction LBs': reaction_lbs
    })
    return metab_df

if __name__ == '__main__':
     
    match_df = pd.read_csv(os.path.join(project_root, 'files','reactions.tsv'), sep='\t',header=0)
    model_base = read_sbml_model(os.path.join(project_root, 'models','EC_3006_1.xml'))
    model_new = read_sbml_model(os.path.join(project_root, 'models','THG-beta2.xml'))  
    
    model_new.solver = 'gurobi'
    opt = model_new.optimize()
    logging.basicConfig(level=logging.DEBUG)
    print('Initial objective value:', opt.objective_value)
    # model_new = pickle.load(open('pipeline/models/THG.pkl', 'rb'))
    # model_base = pickle.load(open('pipeline/models/iEC3006.pkl', 'rb'))
    return_model, metab_df, common_rs = match_exch_rxns(model_new, model_base, match_df) 
    sol = return_model.optimize()
    print('Final new model objective value:', sol.objective_value)
    sol_base = model_base.optimize()
    print('Base model objective value:', sol_base.objective_value)
    write_sbml_model(return_model, os.path.join(project_root, 'models','THG_EC_10_03.xml'))
    #metab_df.to_csv('pipeline/data/EC_THG_not_found_metabolites.csv', index=False)
    metab_df.to_excel(os.path.join(project_root, 'files','EC_THG_not_found_metabolites.xlsx'), index=False)

    #write common reactions to file
    with open(os.path.join(project_root, 'files','common_rs.txt'), 'w') as f:
        for item in common_rs:
            f.write("%s\n" % item)

    print('Done')
