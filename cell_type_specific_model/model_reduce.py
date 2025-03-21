import scipy.io
import numpy as np
import pdb
from cobra.io import read_sbml_model, write_sbml_model
from cobra.flux_analysis import gapfill
from cobra.flux_analysis import flux_variability_analysis
import os
import pickle
import matplotlib.pyplot as plt

'''
Analyzes metabolic reaction activity across samples, reduces the model by removing low-activity reactions, and restores functionality through gap-filling.

Inputs (set by default)
----------
    solutions_mat (str): Path to .mat file with flux solutions matrix, where rows are reactions and columns are samples with flux values.
    sbml_model (str): Path to GEM SBML model file.
    reaction_id_biomass (str): ID of the biomass reaction to check its presence across samples.
    tolerance (float): Tolerance level for gap-filling (e.g., 0.001), defining the lower bound for reaction activity.

Outputs
----------
    model_reduced.xml: SBML file for reduced model excluding zero-activity reactions.
    model_after_gapfilling.xml: SBML file for reduced model after gap-filling.
    compartments_info.txt: Text file with reaction presence and blocked reaction data by compartment.

Description
----------
    1. Loads flux data from iMAT, discretizes fluxes, and computes mean presence of each reaction across samples.
    2. Groups reactions by presence levels and identifies consistently absent or present reactions.
    3. Identifies blocked reactions in the model by compartment and organizes reactions by presence group within compartments.
    4. Creates a reduced model without absent reactions and saves it.
    5. Removes genes and metabolites not associated with any reaction in the reduced model.
    6. Performs gap-filling to ensure functionality in the reduced model, saving the gap-filled version for further analysis.
'''

def model_reduce(model, solutions_imat, output_model_path):

    # Define base path for transcriptomics files
    base_path = os.path.join('transcriptomics')
   
    # #After running iMAT on our model, we have a .mat file with all the reactions as rows, and samples as columns. The values are the fluxes of the reactions in the samples.
    
    #Check if data is in csv or mat format
    if solutions_imat.endswith('.mat'):
        sol_matrix = scipy.io.loadmat(solutions_imat)
        sol_matrix = sol_matrix['all_Solutions_matrix5']
    elif solutions_imat.endswith('.csv'):
        sol_matrix = np.loadtxt(solutions_imat, delimiter=',')
    else:
        print("Data format not recognized, please provide a .mat or .csv file")
        return
    

    #Discretize the matrix, i.e. if the flux is 0, set it to 0, if it is not 0, set it to 1
    bool_matrix = np.where(sol_matrix == 0, 0, 1)

    #Calculate the mean presence of each reaction in the samples
    mean_presence = np.mean(bool_matrix, axis=1)

    #Histogram of the presence of the reactions
    hist, bins = np.histogram(mean_presence, bins=10)
    print(bins)
    print(hist)
    plt.hist(mean_presence, bins=10)
    plt.xlabel('Mean presence of reactions')
    plt.ylabel('Frequency')
    plt.title('Mean presence of reactions in samples')

    mean_presence_plot_path = os.path.join(base_path, 'mean_presence.png')
    plt.savefig(mean_presence_plot_path)

    plt.close()

    #Find reactions with mean presence of 1 and 0
    mean_presence_1 = np.where(mean_presence == 1)
    print("reactions present in all samples: ",len(mean_presence_1[0])) #number of reactions with mean presence of 1, i.e. all samples
    mean_presence_0 = np.where(mean_presence == 0)
    print("reactions present in 0 samples: ", len(mean_presence_0[0])) #number of reactions with mean presence of 0


    #Create a dictionary with the reactions as keys and the mean presence as values
    reaction_presence = {model.reactions[i].id: mean_presence[i] for i in range(len(mean_presence))}

    print("Biomass presence: ", reaction_presence["MAR13082"]) #Check if the biomass reaction is present in all samples

    #Do a flux variability analysis on the model and see which Rs have min and max fluxes of 0
    #Check if a file with the reactions with min and max fluxes of 0 has already been created , called zero_flux_reactions.pk. If it has, load the file, if it hasn't, do the flux variability analysis and save the file
    zero_flux_reactions_path = os.path.join(base_path, 'zero_flux_reactions.pk')

    if os.path.exists(zero_flux_reactions_path):
        
        flux_pk = open(zero_flux_reactions_path, "rb")
        zero_flux_rs = pickle.load(flux_pk)
        flux_pk.close()

    else:
        fva = flux_variability_analysis(model)

        #Get the reactions with min and max fluxes of 0 (blocked reactions)

        zero_flux_rs = fva[(fva['minimum'] == 0) & (fva['maximum'] == 0)].index.tolist()

        #Save the reactions with min and max fluxes of 0 in a file
        flux_pk = open(zero_flux_reactions_path, "wb")
        pickle.dump(zero_flux_rs, flux_pk)
        flux_pk.close()
    
    print("Reactions with min and max fluxes of 0: ", len(zero_flux_rs)) 

    #See which compartment the blocked reactions are in

    #Print also the info in a text file called compartments_info.txt
    compartments_info_path = os.path.join(base_path, 'compartments_info.txt')

    compartments_info = open(compartments_info_path, "w")
    compartments_info.write("Information about blocked reactions in the model: \n\n")

    for comp in model.compartments:
        print("For compartment: ", model.compartments[comp])
        compartments_info.write("For compartment: "+ model.compartments[comp] + "\n")
        reactions = [r for r in model.reactions if r.compartments == {comp}]
        reactions = [r.id for r in reactions]
        zero_flux_rs_comp = [r for r in zero_flux_rs if r in reactions]
        print("Number of reactions with min and max fluxes of 0 in this compartment: ", len(zero_flux_rs_comp))
        compartments_info.write("Number of reactions with min and max fluxes of 0 in this compartment: "+ str(len(zero_flux_rs_comp)) + "\n")

        print("Number of reactions in this compartment: ", len(reactions))
        compartments_info.write("Number of reactions in this compartment: "+ str(len(reactions)) + "\n")

        if len(reactions) != 0:
            print("Percentage of reactions with min and max fluxes of 0 in this compartment: ", len(zero_flux_rs_comp)/len(reactions)*100)
            compartments_info.write("Percentage of reactions with min and max fluxes of 0 in this compartment: "+ str(len(zero_flux_rs_comp)/len(reactions)*100) + "\n")
        else:
            print("Percentage of reactions with min and max fluxes of 0 in this compartment: 0")
            compartments_info.write("Percentage of reactions with min and max fluxes of 0 in this compartment: 0" + "\n")

        print("")
        compartments_info.write("\n")

    #Define the 5 groups: group1: 0 mean presence, group2: >0-0.1 mean presence, group3: 0.1-0.9 mean presence, group4: 0.9-<1 mean presence, group5: 1 mean presence

    #Group all reactions in the model in these 5 groups, using reaction_presence dictionary

    group1 = [] #0 mean presence
    group2 = [] #>0-0.1 mean presence
    group3 = [] #0.1-0.9 mean presence
    group4 = [] #0.9-<1 mean presence
    group5 = [] #1 mean presence

    for k, v in reaction_presence.items():
        if v == 0:
            group1.append(k)
        elif v > 0 and v <= 0.1:
            group2.append(k)
        elif v > 0.1 and v < 0.9:
            group3.append(k)
        elif v >= 0.9 and v < 1:
            group4.append(k)
        elif v == 1:
            group5.append(k)

    #Print the number of reactions in each group

    print("Number of reactions in each group: \nGroup1: ", len(group1), "\nGroup2: ", len(group2), "\nGroup3: ", len(group3), "\nGroup4: ", len(group4), "\nGroup5: ", len(group5))
    compartments_info.write("Number of reactions in each group: \nGroup1: "+ str(len(group1)) + "\nGroup2: "+ str(len(group2)) + "\nGroup3: "+ str(len(group3)) + "\nGroup4: "+ str(len(group4)) + "\nGroup5: "+ str(len(group5)) + "\n")
    
    #For each compartment, calculate the number and % of reactions in each group

    compartments_info.write("\n\nInformation about the groups in each compartment: \n\n")
    for comp in model.compartments:
        print("For compartment: ", model.compartments[comp])
        compartments_info.write("For compartment: "+ model.compartments[comp] + "\n")
        reactions = [r for r in model.reactions if r.compartments == {comp}]
        reactions = [r.id for r in reactions]

        #Find intersection between groups and reactions using set intersection
        group1_comp = set(group1).intersection(reactions)
        group2_comp = set(group2).intersection(reactions)
        group3_comp = set(group3).intersection(reactions)
        group4_comp = set(group4).intersection(reactions)
        group5_comp = set(group5).intersection(reactions)


        print("Number of reactions in this compartment: ", len(reactions))
        print("Number of reactions in group1: ", len(group1_comp))
        print("Number of reactions in group2: ", len(group2_comp))
        print("Number of reactions in group3: ", len(group3_comp))
        print("Number of reactions in group4: ", len(group4_comp))
        print("Number of reactions in group5: ", len(group5_comp))

        compartments_info.write("Number of reactions in this compartment: "+ str(len(reactions)) + "\n")
        compartments_info.write("Number of reactions in group1: "+ str(len(group1_comp)) + "\n")
        compartments_info.write("Number of reactions in group2: "+ str(len(group2_comp)) + "\n")
        compartments_info.write("Number of reactions in group3: "+ str(len(group3_comp)) + "\n")
        compartments_info.write("Number of reactions in group4: "+ str(len(group4_comp)) + "\n")
        compartments_info.write("Number of reactions in group5: "+ str(len(group5_comp)) + "\n")


        if len(reactions) != 0:

            print("Percentage of reactions in group1: ", len(group1_comp)/len(reactions)*100)
            print("Percentage of reactions in group2: ", len(group2_comp)/len(reactions)*100)
            print("Percentage of reactions in group3: ", len(group3_comp)/len(reactions)*100)
            print("Percentage of reactions in group4: ", len(group4_comp)/len(reactions)*100)
            print("Percentage of reactions in group5: ", len(group5_comp)/len(reactions)*100)

            compartments_info.write("Percentage of reactions in group1: "+ str(len(group1_comp)/len(reactions)*100) + "\n")
            compartments_info.write("Percentage of reactions in group2: "+ str(len(group2_comp)/len(reactions)*100) + "\n")
            compartments_info.write("Percentage of reactions in group3: "+ str(len(group3_comp)/len(reactions)*100) + "\n")
            compartments_info.write("Percentage of reactions in group4: "+ str(len(group4_comp)/len(reactions)*100) + "\n")
            compartments_info.write("Percentage of reactions in group5: "+ str(len(group5_comp)/len(reactions)*100) + "\n")

        else:
            print("Percentage of reactions in group1: 0")
            print("Percentage of reactions in group2: 0")
            print("Percentage of reactions in group3: 0")
            print("Percentage of reactions in group4: 0")
            print("Percentage of reactions in group5: 0")

            compartments_info.write("Percentage of reactions in group1: 0" + "\n")
            compartments_info.write("Percentage of reactions in group2: 0" + "\n")
            compartments_info.write("Percentage of reactions in group3: 0" + "\n")
            compartments_info.write("Percentage of reactions in group4: 0" + "\n")
            compartments_info.write("Percentage of reactions in group5: 0" + "\n")

        print("")
        compartments_info.write("\n")

    compartments_info.close()

    #Create a model without the reactions in group1 (0 mean presence)
    model_reduced = model.copy()

    model_reduced.remove_reactions(group1) 

    print("Number of reactions in reduced model: ", len(model_reduced.reactions))
    print("Number of genes in reduced model: ", len(model_reduced.genes))
    print("Number of metabolites in reduced model: ", len(model_reduced.metabolites))

    #From this model, remove the genes and metabolites that are not associated with any reaction

    mets_to_remove = [met for met in model_reduced.metabolites if len(met.reactions) == 0]
    model_reduced.remove_metabolites(mets_to_remove)

    genes_to_remove = [gene for gene in model_reduced.genes if len(gene.reactions) == 0]
    for gene in genes_to_remove:
        model_reduced.genes.remove(gene)

    print("After removing genes and metabolites not associated with any reaction: ")
    print("Number of reactions in reduced model: ", len(model_reduced.reactions))
    print("Number of genes in reduced model: ", len(model_reduced.genes))
    print("Number of metabolites in reduced model: ", len(model_reduced.metabolites))


    #Save the reduced model
    write_sbml_model(model_reduced, output_model_path)


    #Now perform gapfilling on model_reduced using the complete model as the universal model

    tolerance = 0.001 #Same one used in PTR code

    #Try to optimize the model before gapfilling
    solution = model_reduced.optimize()
    
    print("Optimization status:", solution.status)
    print("Objective value:", solution.objective_value)
        
    try:
        gapfill_solution = gapfill(model_reduced, model, iterations=5, lower_bound=tolerance, penalties=None, exchange_reactions=False, demand_reactions=False)
        for solution in gapfill_solution:
            for reaction in solution:
                model_reduced.add_reactions([reaction])
                print("Reaction added: ", reaction.id)

        # Check the updated reaction count after gapfilling
        print("Number of reactions in the model for gapfilling after gapfilling: ", len(model_reduced.reactions))

        #Check if gapfill solution is just a list of empty lists (i.e. gapfilling was not successful)
        if gapfill_solution == [[]]*len(gapfill_solution):
            print("Gapfilling was not successful")
        else:
            print("Gapfilling was successful, creating gapfilled model")
            gapfilled_model_path = os.path.join(base_path, 'model_after_gapfilling.xml')
            write_sbml_model(model_reduced, gapfilled_model_path)

    except Exception as e:
        print("An error occurred during the gapfilling process:", str(e))

    
    

def main():

    base_path = os.path.join('transcriptomics') #Define base path for transcriptomics files
    model_path = os.path.join(base_path, 'model_full_THG_endoA.xml')

    output_model_path = os.path.join(base_path, 'model_reduced_endoA_gimme.xml')
    solutions = os.path.join(base_path, 'allsolutions_gimme_parallel.csv')

    model = read_sbml_model(model_path)
    
    model_reduce(model, solutions, output_model_path)

if __name__ == '__main__':
    main()



