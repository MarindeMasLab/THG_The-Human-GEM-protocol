from cobra.io import read_sbml_model, load_json_model, save_json_model
from cobra.util import create_stoichiometric_matrix
from cobra.manipulation.delete import prune_unused_reactions, prune_unused_metabolites
import networkx as nx
from tqdm import tqdm
import json
import os


def find_network_components(
    model,
    cleanup=False,
    cleanup_save_path=None,
    verbose=True,
    visualize=False,
    viz_output_path=None,
):
    """
    Analyze the network connectivity of a metabolic model by finding weakly connected components.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to analyze.
    cleanup : bool, optional
        Whether to remove blocked reactions before analysis. Default is False.
    cleanup_save_path : str, optional
        Path to save the cleaned model (only used if cleanup=True).
    verbose : bool, optional
        Whether to print detailed component information. Default is True.
    visualize : bool, optional
        Whether to generate HTML visualization of components. Default is False.
    viz_output_path : str, optional
        Path for the HTML visualization file. Required if visualize=True.

    Returns
    -------
    dict
        Dictionary containing:
        - 'model': The (possibly cleaned) model used for analysis
        - 'graph': NetworkX DiGraph of the network
        - 'components': List of weakly connected components (sorted by size, descending)
        - 'largest_component': The largest component
        - 'is_fully_connected': Boolean indicating if network is fully connected
        - 'component_info': List of dicts with details about each component
        - 'viz_path': Path to HTML visualization (if visualize=True)
    """
    # Make a copy to avoid modifying the original
    model_copy = model.copy()

    if cleanup:
        from cobra.flux_analysis import find_blocked_reactions

        if verbose:
            print(f"Number of reactions before cleanup: {len(model_copy.reactions)}")
        blocked_reactions = find_blocked_reactions(model_copy)
        if verbose:
            print(f"Removing {len(blocked_reactions)} blocked reactions")
        model_copy.remove_reactions(blocked_reactions)

        if cleanup_save_path:
            save_json_model(model_copy, cleanup_save_path)

    model_copy, _ = prune_unused_metabolites(model_copy)

    if verbose:
        print(f"Number of metabolites: {len(model_copy.metabolites)}")
        print(f"Number of reactions: {len(model_copy.reactions)}")

    # Build a directed bipartite graph: metabolites <-> reactions
    G = nx.DiGraph()

    # Add metabolite nodes
    for met in model_copy.metabolites:
        G.add_node(met.id, bipartite=0)

    # Add reaction nodes
    for rxn in model_copy.reactions:
        G.add_node(rxn.id, bipartite=1)

    # Add edges: metabolite -> reaction (if consumed), reaction -> metabolite (if produced)
    for rxn in model_copy.reactions:
        for met, coeff in rxn.metabolites.items():
            if coeff < 0:
                G.add_edge(met.id, rxn.id)  # consumed
            elif coeff > 0:
                G.add_edge(rxn.id, met.id)  # produced

    # Precompute compartment mappings for metabolites and reactions
    met_compartments = {met.id: met.compartment for met in model_copy.metabolites}
    rxn_compartments = {rxn.id: set(rxn.compartments) for rxn in model_copy.reactions}

    # Get all weakly connected components
    components = list(nx.weakly_connected_components(G))
    components = sorted(components, key=len, reverse=True)  # Order by size descending

    if verbose:
        print(f"Number of weakly connected components: {len(components)}")

    # Create a set of reaction IDs for fast lookup
    reaction_ids = set(rxn_compartments.keys())

    # Collect component information
    component_info = []
    iterator = tqdm(enumerate(components)) if verbose else enumerate(components)

    for idx, comp in iterator:
        # Count reactions in this component by comparing node IDs
        reaction_count = sum(1 for node in comp if node in reaction_ids)

        # Efficiently collect compartments present in this component
        met_ids = comp & set(met_compartments.keys())
        rxn_ids = comp & reaction_ids
        compartments = set()
        compartments.update(met_compartments[mid] for mid in met_ids)
        for rid in rxn_ids:
            compartments.update(rxn_compartments[rid])

        comp_info = {
            "index": idx + 1,
            "total_nodes": len(comp),
            "reaction_count": reaction_count,
            "metabolite_count": len(comp) - reaction_count,
            "compartments": sorted(compartments),
            "nodes": comp,
        }
        component_info.append(comp_info)

        if verbose:
            print(
                f"Component {idx+1}: {len(comp)} nodes, {reaction_count} reactions, compartments: {sorted(compartments)}"
            )

    # Get the largest weakly connected component
    largest_component = components[0] if components else set()
    is_fully_connected = len(largest_component) == G.number_of_nodes()

    if verbose:
        print(f"Largest weakly connected component size: {len(largest_component)}")
        print(f"Total nodes in network: {G.number_of_nodes()}")

        if is_fully_connected:
            print("The network is fully connected.")
        else:
            print("The network is not fully connected.")

    results = {
        "model": model_copy,
        "graph": G,
        "components": components,
        "largest_component": largest_component,
        "is_fully_connected": is_fully_connected,
        "component_info": component_info,
    }

    # Generate visualization if requested
    if visualize:
        if viz_output_path is None:
            raise ValueError("viz_output_path must be provided when visualize=True")
        viz_path = visualize_small_components(results, viz_output_path, verbose=verbose)
        results["viz_path"] = viz_path

    return results


def visualize_small_components(
    results, output_path, max_metabolite_degree=10, exclude_largest=True, verbose=True
):
    """
    Create an interactive HTML visualization of network components using vis.js.

    Parameters
    ----------
    results : dict
        Output from find_network_components() function.
    output_path : str
        Path where the HTML file will be saved.
    max_metabolite_degree : int, optional
        Maximum degree for metabolites. Higher degree metabolites are filtered out. Default is 10.
    exclude_largest : bool, optional
        Whether to exclude the largest component from visualization. Default is True.
    verbose : bool, optional
        Whether to print information about filtered metabolites. Default is True.

    Returns
    -------
    str
        Path to the created HTML file.
    """
    model = results["model"]
    G = results["graph"]
    components = results["components"]

    # Filter components
    if exclude_largest and len(components) > 0:
        small_components = components[1:]  # Exclude largest
        if verbose:
            print(
                f"Number of small components (excluding largest): {len(small_components)}"
            )
    else:
        small_components = components
        if verbose:
            print(f"Number of components to visualize: {len(small_components)}")

    if not small_components:
        if verbose:
            print("No components to visualize.")
        return None

    # Collect all filtered metabolites across all components
    all_filtered_metabolites = set()

    # Build component_data directly from NetworkX graphs (nodes + edges)
    component_data = []
    for idx, comp in enumerate(small_components):
        nodes_to_plot = set(comp)
        G_small = G.subgraph(nodes_to_plot).copy()

        # Filter out metabolites with more than max_metabolite_degree connections
        nodes_to_remove = []
        for node in G_small.nodes():
            if node in {met.id for met in model.metabolites}:
                degree = G.degree(node)  # Use degree from original graph G
                if degree > max_metabolite_degree:
                    nodes_to_remove.append(node)

        # Collect filtered metabolites for this component
        for node_id in nodes_to_remove:
            met = model.metabolites.get_by_id(node_id)
            all_filtered_metabolites.add((met.name, met.id, G.degree(node_id)))

        # Remove high-degree metabolites from the subgraph
        G_small.remove_nodes_from(nodes_to_remove)
        # Update nodes_to_plot to exclude removed nodes
        nodes_to_plot = set(G_small.nodes())

        nodes_data = []
        # precompute metabolite and reaction id sets for speed
        met_ids = {met.id for met in model.metabolites}
        rxn_ids = {rxn.id for rxn in model.reactions}

        for node in G_small.nodes():
            if node in met_ids:
                met = model.metabolites.get_by_id(node)
                comp_name = met.compartment if hasattr(met, "compartment") else ""
                label_id = f"Met: {met.id} [{comp_name}]"
                label_name = f"Met: {met.name} [{comp_name}]"
                color = "#97C2FC"
                nodes_data.append(
                    {
                        "id": node,
                        "label_id": label_id,
                        "label_name": label_name,
                        "title_id": label_id,
                        "title_name": label_name,
                        "color": color,
                        "shape": "dot",
                        "size": 14,
                        "type": "metabolite",
                    }
                )
            elif node in rxn_ids:
                rxn = model.reactions.get_by_id(node)
                label = f"Rxn: {rxn.id}"
                nodes_data.append(
                    {
                        "id": node,
                        "label_id": label,
                        "label_name": label,
                        "title_id": label,
                        "title_name": label,
                        "color": "#FB7E81",
                        "shape": "dot",
                        "size": 14,
                        "type": "reaction",
                    }
                )
            else:
                label = str(node)
                nodes_data.append(
                    {
                        "id": node,
                        "label_id": label,
                        "label_name": label,
                        "title_id": label,
                        "title_name": label,
                        "color": "#DDDDDD",
                        "shape": "dot",
                        "size": 14,
                        "type": "other",
                    }
                )

        edges_data = []
        for u, v in G_small.edges():
            edges_data.append({"from": u, "to": v})

        component_data.append(
            {
                "nodes": nodes_data,
                "edges": edges_data,
                "n_nodes": len(G_small.nodes()),
                "n_edges": len(G_small.edges()),
            }
        )

    # Print the names of filtered metabolites
    if verbose:
        if all_filtered_metabolites:
            print(
                f"\nFiltered out {len(all_filtered_metabolites)} metabolites with >{max_metabolite_degree} connections:"
            )
            for name, met_id, degree in sorted(
                all_filtered_metabolites, key=lambda x: x[2], reverse=True
            ):
                print(f"  {name} ({met_id}) - {degree} connections")
        else:
            print(
                f"\nNo metabolites were filtered out (all have ≤{max_metabolite_degree} connections)"
            )

    # Load HTML template
    template_path = os.path.join(
        os.path.dirname(__file__), "component_visualization_template.html"
    )
    with open(template_path, "r") as f:
        html_template = f.read()

    # Build component buttons HTML
    start_idx = 2 if exclude_largest else 1
    buttons_html = ""
    for idx, comp in enumerate(small_components):
        comp_num = idx + start_idx
        buttons_html += f'  <button class="tab-button" data-target="Component{comp_num}">Component {comp_num} <br><small>{len(comp)} nodes</small></button>\n'

    # Build component divs HTML
    divs_html = ""
    for idx, comp in enumerate(small_components):
        comp_num = idx + start_idx
        div_id = f"network{comp_num}"
        divs_html += f'  <div id="Component{comp_num}" class="tabcontent">\n'
        divs_html += f"    <h3>Component {comp_num} ({len(comp)} nodes)</h3>\n"
        divs_html += f'    <div id="{div_id}" class="network-container"></div>\n'
        divs_html += f"  </div>\n"

    # Replace placeholders in template
    html = html_template.replace("{{COMPONENT_BUTTONS}}", buttons_html)
    html = html.replace("{{COMPONENT_DIVS}}", divs_html)
    html = html.replace("{{NETWORK_CONFIGS}}", json.dumps(component_data))
    html = html.replace("{{START_IDX}}", str(start_idx))

    with open(output_path, "w") as f:
        f.write(html)

    if verbose:
        print(f"\nTabbed HTML file created (toggleable id/name): {output_path}")

    return output_path


if __name__ == "__main__":
    # Example usage
    model_name = "endoA_251119"

    # model = read_sbml_model(f"../models/{model_name}.xml")
    model = load_json_model(f"../models/{model_name}.json")

    CLEANUP = True
    VISUALIZE = True

    # Run analysis
    results = find_network_components(
        model,
        cleanup=CLEANUP,
        cleanup_save_path=f"../models/{model_name}_clean.json" if CLEANUP else None,
        verbose=True,
        visualize=VISUALIZE,
        viz_output_path=f"{model_name}_components_tabs.html" if VISUALIZE else None,
    )
