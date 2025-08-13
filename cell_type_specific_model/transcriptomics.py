#Script for transcriptomics analysis
import numpy as np
import pdb
from cobra.io import read_sbml_model
import os
import re
import xml.etree.ElementTree as ET
import pandas as pd
from model_reduce import model_reduce
from gimme_parallel import gimme_parallel as gimme
import multiprocessing
import copy
import sys

# Determine the current file's directory and the project root.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(current_dir, "..")
# Add the project root to sys.path to access top-level folders like 'functions' and 'models'
if project_root not in sys.path:
    sys.path.append(project_root)

def biomass_fix(model):
    """
    Adjusts the biomass reaction lower bound to 80% of the optimal biomass flux.
    
    This function:
      - Prints the original bounds of the biomass reaction.
      - Optimizes the model (FBA) to determine the optimal biomass flux.
      - Sets the lower bound of the biomass reaction (MAR13082) to 80% of that optimal flux.
      - Re-optimizes the model and prints updated bounds and solution status.
      
    Parameters:
        model (cobra.Model): The metabolic model.
        
    Returns:
        cobra.Model: The updated model with the modified biomass reaction lower bound.
    """
     
    print("Original biomass reaction bounds: ", model.reactions.get_by_id("MAR13082").lower_bound, model.reactions.get_by_id("MAR13082").upper_bound) 
    
    # Run FBA to get the optimal biomass flux.
    solution = model.optimize()
    optimal_biomass_flux = solution.objective_value
    print("Optimal biomass flux: ", optimal_biomass_flux)
    print("Solution status: ", solution.status)

    # Set a lower bound for the biomass reaction to 80% of the optimal biomass flux
    model.reactions.get_by_id("MAR13082").lower_bound = 0.8*optimal_biomass_flux

    # Re-optimize the model after changing the biomass bounds.
    solution = model.optimize()
    optimal_biomass_flux = solution.objective_value
    print("After changing bounds: \nOptimal biomass flux: ", optimal_biomass_flux)
    print("Solution status: ", solution.status)
    print("Lower bound for biomass reaction: ", model.reactions.get_by_id("MAR13082").lower_bound)
    print("Final biomass reaction bounds: ", model.reactions.get_by_id("MAR13082").lower_bound, model.reactions.get_by_id("MAR13082").upper_bound) 

    return model

def sgpr_rules(model):
    """
    Extracts sGPR rules from reaction annotations in the model.
    
    For each reaction in the model, this function searches through its annotations
    for a key containing 'sGPR' and extracts the rule (assumed to be the last component
    of a slash-separated string).
    
    Parameters:
        model (cobra.Model): The metabolic model.
    
    Returns:
        list: A list of sGPR rules corresponding to the reactions in the model.
              If a reaction does not have an sGPR rule, its corresponding entry is None.
    """
    sgpr_rules = [None] * len(model.reactions)

    # Extract sGPR rules from the annotations of reactions
    for i, reaction in enumerate(model.reactions):
        if hasattr(reaction, 'annotation') and reaction.annotation:
            for key, resources in reaction.annotation.items():
                if 'sGPR' in key:
                    sgpr_rule = resources.split('/')[-1]
                    sgpr_rules[i] = sgpr_rule

    return sgpr_rules

def ensemble(model):
    """
    Extracts Ensembl IDs from the gene identifiers in the model.
    
    For each gene in the model, searches for an Ensembl ID using a regex pattern.
    If an Ensembl ID is found, it is added to the returned list; otherwise, an empty
    string is added.
    
    Parameters:
        model (cobra.Model): The metabolic model.
    
    Returns:
        list: A list of Ensembl IDs (or empty strings) corresponding to the model's genes.
    """
    ensembl_model = []
    pattern = r'ENSG\d+'

    # Extract Ensembl IDs from gene identifiers in the model
    for gene in model.genes:
        match = re.search(pattern, gene.id)
        if match:
            ensembl_model.append(match.group(0))
        else:
            ensembl_model.append('')

    return ensembl_model

def index_to_ensembl(model, ensembl_model):
    """
    Updates the gene-protein-reaction (GPR) rules in the model by replacing index tokens 
    with corresponding Ensembl IDs.
    
    The function searches for tokens in the form x(index) in each reaction's GPR rule,
    retrieves the corresponding Ensembl ID from the provided list, and replaces the token.
    
    Parameters:
        model (cobra.Model): The metabolic model.
        ensembl_model (list): A list of Ensembl IDs corresponding to gene indices.
    
    Returns:
        list: A list of updated GPR rules with indices replaced by Ensembl IDs.
    """

    new_rules3 = [None] * len(model.reactions)

    # Update GPR rules to replace indices with corresponding Ensembl IDs 
    for j, reaction in enumerate(model.reactions):
        rule = reaction.gene_reaction_rule
        tokens = re.findall(r'x\((\d+)\)', rule)
        new_rule = rule

        for token in tokens:
            index = int(token)
            ensg_id = ensembl_model[index] if index < len(ensembl_model) else ''
            new_rule = new_rule.replace(f'x({token})', f'({ensg_id})')

        # Clean up the GPR
        new_rule = new_rule.strip()

        new_rules3[j] = new_rule

    return new_rules3

def extract_pairs(model_path):
    """
    Parses an XML file (model) to extract pairs of HGNC gene symbols and Ensembl IDs.
    
    The function uses ElementTree to navigate through the XML structure, extracting 
    the HGNC symbol and Ensembl ID from specific resource links in the file. In cases 
    where duplicate symbols are encountered with different Ensembl IDs, it prints a warning.
    
    Parameters:
        model_path (str): Path to the XML file containing model information.
    
    Returns:
        dict: A dictionary mapping HGNC symbols to Ensembl IDs. If duplicates occur, 
              the value may be a list of Ensembl IDs.
    """
    tree = ET.parse(model_path)
    root = tree.getroot()
    symbol_and_ensembl = {}

    # Extract pairs of HGNC symbols and Ensembl IDs from an XML file

    # Iterate over all Description elements in the XML.

    for description in root.findall('.//{http://www.w3.org/1999/02/22-rdf-syntax-ns#}Description'):
        hgnc_symbol = ''
        ensg = ''

        # Iterate over list items in the Description.
        for li in description.findall('.//{http://www.w3.org/1999/02/22-rdf-syntax-ns#}li'):
            resource = li.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource')

            if resource:
                if 'hgnc.symbol/' in resource or 'hgcn.symbol/' in resource:
                    hgnc_symbol = resource.split('/')[-1]
                elif 'ensembl/ENSG' in resource:
                    ensg = resource.split('/')[-1]

        # If both HGNC and Ensembl ID are found, add them to the dictionary.
        if hgnc_symbol and ensg:
            if hgnc_symbol in symbol_and_ensembl:
                # If duplicate symbol, check whether it is the same Ensembl ID.
                if symbol_and_ensembl[hgnc_symbol] == ensg:
                    continue
                # Convert to list if needed and add the new Ensembl ID.
                if not isinstance(symbol_and_ensembl[hgnc_symbol], list):
                    symbol_and_ensembl[hgnc_symbol] = [symbol_and_ensembl[hgnc_symbol]]
                symbol_and_ensembl[hgnc_symbol].append(ensg)     
            else:
                symbol_and_ensembl[hgnc_symbol] = ensg

    return symbol_and_ensembl


def sgpr_to_ensembl(model, symbol_and_ensembl):
    """
    Converts GPR rules in the model from gene symbols/indices to Ensembl IDs using 
    the provided mapping.
    
    For each reaction, the function scans the gene reaction rule for gene tokens 
    (ignoring logical operators). It then replaces these tokens with the corresponding 
    Ensembl IDs from the symbol_and_ensembl dictionary. Unmatched genes are reported.
    
    Parameters:
        model (cobra.Model): The metabolic model.
        symbol_and_ensembl (dict): Dictionary mapping HGNC symbols to Ensembl IDs.
    
    Returns:
        list: A list of updated GPR rules with gene tokens replaced by Ensembl IDs.
    """

    new_rules3 = [None] * len(model.reactions)
    unmatched_indices = []       # Gene tokens not found in the mapping.
    unmatched_rules = []         # GPR rules with missing Ensembl IDs.
    unmatched_row_numbers = []   # Row numbers of reactions with unmatched genes.

    for j, reaction in enumerate(model.reactions):
        rule = reaction.gene_reaction_rule
        new_rule = rule

        # Find all words (potential gene tokens) in the rule.
        genes = re.findall(r'\w+', rule)
        # Remove logical operator tokens.
        genes = [gene for gene in genes if gene not in ['or', 'and']]
        for gene in genes:
            # Check if gene is already an Ensembl ID or is present in the mapping.
            if gene in symbol_and_ensembl.values():
                ensg_id = gene
            elif gene.startswith('ENSG'):
                ensg_id = gene # Already in Ensembl format.
            elif gene in symbol_and_ensembl.keys():
                # If multiple Ensembl IDs exist, choose the first.
                if isinstance(symbol_and_ensembl[gene], list):
                    print(f"Gene {gene} has multiple Ensembl IDs: {symbol_and_ensembl[gene]}")
                    ensg_id = symbol_and_ensembl[gene][0]
                else:   
                    ensg_id = symbol_and_ensembl[gene] 
            else:
                ensg_id = ''
                print(f"Gene {gene} not found in the Ensembl model")
                unmatched_indices.append(gene)
                unmatched_rules.append(rule)
                unmatched_row_numbers.append(j)
            # Replace the gene token in the rule with the Ensembl ID.
            new_rule = re.sub(r'\b{}\b'.format(re.escape(gene)), ensg_id, new_rule)

        new_rule = new_rule.strip()
        new_rules3[j] = new_rule

    return new_rules3


class CustomExpression:
    """
    Implements a custom expression object that allows arithmetic operations
    
    The class defines custom behavior for logical AND (&) and OR (|) operators.
    """
    def __init__(self, value):
        """
        Initialize with a numeric expression value.
        
        Parameters:
            value (float): The expression value.
        """
        self.value = value

    def __and__(self, other):
        """
        Implements the AND logic.
        
        Returns:
            CustomExpression: For AND, if either value is 0, returns the maximum;
                              otherwise, returns the minimum.
        """
        if self.value == 0 or other.value == 0:
            return CustomExpression(max(self.value, other.value))
        else:
            return CustomExpression(min(self.value, other.value))

    def __or__(self, other):
        """
        Implements the OR logic.
        
        Returns:
            CustomExpression: The maximum of the two values.
        """
        return CustomExpression(max(self.value, other.value))

    def __repr__(self):
        """Make the object printable for debugging."""
        return str(self.value)

def preprocess_rule(rule, ensembl_map):
    """
    Preprocesses a gene reaction rule to standardize Ensembl IDs and embed expression values.
    
    The function performs the following steps:
      - Finds all Ensembl IDs in the rule that match the pattern ENSG followed by 11 digits
        and optional extra characters, truncating any extra characters.
      - Replaces these Ensembl IDs with a call to CustomExpression using the corresponding
        expression value from ensembl_map.
      - Converts logical operator tokens ('[', ']', 'or', 'and') into standard symbols.
    
    Parameters:
        rule (str): The gene reaction rule as a string.
        ensembl_map (dict): Dictionary mapping Ensembl IDs to expression values.
    
    Returns:
        str: The preprocessed rule ready for evaluation.
    """

    # Find all Ensembl IDs and remove any extra characters beyond the standard 15 characters (ENSG + 11 digits)
    matches = re.findall(r'ENSG\d{11}[A-Za-z0-9]*', rule)
    # Process the matches by length
    matches.sort(key=len, reverse=True) 
    for match in matches:
        if len(match) > 15:
            new_match = match[:15]
            rule = rule.replace(match, new_match)

    # Replace gene IDs with their corresponding expression values
    matches = re.findall(r'ENSG\d+', rule)
    matches.sort(key=len, reverse=True)
    for match in matches:
        new_match = match[:15] 
        rule = rule.replace(match, new_match) 
        expression_value = ensembl_map.get(new_match, -1)  # Default to -1 if not found
        rule = rule.replace(new_match, f'CustomExpression({expression_value})')     
        
    # Perform regex replacements for logical operators
    rule = re.sub(r'\[', '(', rule)
    rule = re.sub(r'\]', ')', rule)
    rule = re.sub(r'\bor\b', '|', rule)
    rule = re.sub(r'\band\b', '&', rule)

    return rule

def process_expression_data(new_rules3, lung_data):
    """
    Processes gene expression data for multiple samples and evaluates the gene reaction rules.
    
    Steps:
      - Loads expression data from a CSV file.
      - For each sample, builds an Ensembl map from the expression data.
      - Preprocesses the gene reaction rules using the expression map.
      - Evaluates the preprocessed rules in a restricted environment (only CustomExpression available).
      - Stores the evaluated expression values in a 2D numpy array.
    
    Parameters:
        new_rules3 (list): List of preprocessed gene reaction rules.
        lung_data (str): Path to the CSV file containing expression data. The first column is gene IDs.
    
    Returns:
        - rules_samples (list): A list where each entry corresponds to a reaction's list of rules per sample.
        - expressionRxns (np.array): 2D array (reactions x samples) of evaluated expression values.
    """
    # Load the expression data
    lung_data = pd.read_csv(lung_data)
    lung_data_cell = lung_data.values
    num_samples = lung_data.shape[1] - 1  # Exclude the first column for gene IDs
    rules_samples = [[None] * num_samples for _ in range(len(new_rules3))]

    for sample in range(1, lung_data.shape[1]):  # Skip the first column
        print(f"Processing sample {sample} out of {num_samples}")
        # Build an expression mapping: gene ID -> expression value for the current sample.
        ensembl_map = {lung_data_cell[i, 0]: lung_data_cell[i, sample] for i in range(lung_data.shape[0])}

        updated_rules = new_rules3.copy()
        for i, rule in enumerate(updated_rules):
            if isinstance(rule, str):
                # Preprocess the rule
                rule = preprocess_rule(rule, ensembl_map)
                updated_rules[i] = rule

        # Store the processed rules for this sample
        for i in range(len(updated_rules)):
            rules_samples[i][sample - 1] = updated_rules[i]

    # Evaluate the preprocessed rules and calculate reaction expression scores.
    print("Evaluating gene expression levels using gene rules...")
    expressionRxns = np.zeros((len(rules_samples), num_samples))

    # Create a restricted evaluation environment with only CustomExpression available.
    local_env = {'CustomExpression': CustomExpression}
    for sample in range(num_samples):
        print(f"Evaluating sample {sample + 1} out of {num_samples}")
        rules3 = [rules_samples[i][sample] for i in range(len(rules_samples))]
        frxnscores = np.zeros(len(rules3))

        for k, rule in enumerate(rules3):
            if not rule or rule.strip() == '':
                frxnscores[k] = -1  # Default to -1 for empty rules, will be recognized by reconstruction algorithm
            else:
                try:
                    frxnscores[k] = eval(rule, {"__builtins__": {}}, local_env).value
                except Exception as e:
                    print(f"Error evaluating rule: {rule} ///, rule number {k}, setting score to -1. Error: {e}")
                    frxnscores[k] = -1

        expressionRxns[:, sample] = frxnscores

    print("Finished processing expression data.")
    return rules_samples, expressionRxns


#Main function

def main():
    """
    Main function to execute the transcriptomics analysis pipeline.
    
    Workflow:
      - Defines file paths for the model and expression data.
      - Loads the metabolic model (SBML format) and applies biomass_fix.
      - Extracts gene symbol and Ensembl ID pairs from the model XML.
      - Converts sGPR rules to use Ensembl IDs.
      - Processes the expression data (preprocess and evaluate gene rules).
      - Saves the processed expression values to a CSV file.
      - Reads exchange reaction IDs relevant in the specific cell type.
      - Runs the parallelized GIMME algorithm to obtain flux solutions.
      - Calls model_reduce to generate a tailored model based on the flux solutions.
      - Prints the output location of the tailored model.
    """
    # Define paths using os.path.join for cross-platform compatibility
    model_path = os.path.join(project_root, 'models','model_full_THG_round2.xml') # Output model after adding transport reactions
    lung_data= os.path.join(project_root, 'files', 'lung_data_gtex.csv') #GTEx expression data
    model = read_sbml_model(model_path)

    model = biomass_fix(model)
    # Extract gene symbol and Ensembl ID pairs from the model XML.
    symbol_and_ensembl = extract_pairs(model_path)
    # Convert sGPR rules to use Ensembl IDs.
    new_rules2 = sgpr_to_ensembl(model, symbol_and_ensembl)
    # Process the expression data to generate evaluated expression values.
    rules_samples, expressionRxns = process_expression_data(new_rules2, lung_data)

    # Convert rules_samples to a more memory-efficient data type 
    rules_samples = np.array(rules_samples, dtype=object)  # Use dtype=object for mixed types
    expressionRxns = np.array(expressionRxns, dtype=np.float32)  # Use float32 for smaller memory footprint

    expressionRxns_df = pd.DataFrame(expressionRxns)
    expression_rxns_path = os.path.join(project_root, 'files', 'expressionRxns.csv')
    expressionRxns_df.to_csv(expression_rxns_path, index=False, header=False)

    #Exchange reactions common between the general model and the cell type specific model (obtained from match_exch_rxns script)
    exchange_rs = pd.read_csv(os.path.join(project_root, 'files','common_rs.txt'), header=None).values.flatten().tolist()

    # Create a copy for flux computations so that the original model remains unmodified
    #model_flux = model.copy()
    model_flux = copy.deepcopy(model)
    all_solutions = gimme(model_flux, expressionRxns, exchange_rs, num_workers=1)
    
    #Tailor the model based on the flux solutions from the reconstruction algorithm
    output_path = os.path.join(project_root, 'models', 'model_THG_endoA_tailored.xml') # Output path for the tailored model
    
    model_reduce(model, all_solutions, output_path)

    print("Model tailored and saved to ", output_path)


if __name__ == "__main__":
    main()
