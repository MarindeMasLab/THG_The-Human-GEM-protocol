#!/usr/bin/env python3
"""
Unified Pathway Network Visualization Tool

Generate different types of network visualizations for any pathway:
  0: Generate all visualizations (default)
  1: D3.js interactive HTML (force-directed layout)
  2: Plotly interactive network (spring layout with FVA)
  3: GraphViz static image (hierarchical layout, PNG/PDF)

Usage:
  python3 visualize.py [option] [--config CONFIG_PATH] [--model MODEL_PATH]

Examples:
  python3 visualize.py                                    # Generate all visualizations (default pathway)
  python3 visualize.py 0 --config config/config_Syndecan1_HS.json --model ../models/endoA_250917_3_with_Syndecan1_HS.json
  python3 visualize.py 1 --config config/config_Syndecan1_HS.json  # D3.js HTML only
  python3 visualize.py 2                                  # Plotly interactive only (default pathway)
"""

import sys
from pathlib import Path

# Add parent directory to path for shared functions module
sys.path.insert(0, str(Path(__file__).parent.parent))
import sys
import os
import json
import argparse
from collections import defaultdict

# Import config loader
from functions.config import load_config, get_model_paths

# Global variables (will be set by command line args or defaults)
PATHWAY_CONFIG = None
PATHWAY_NAME = "cytoskeleton"
MODEL_PATH = None
VIZ_OUTPUT_DIR = "figures"


def load_pathway_config(config_path):
    """Load pathway configuration from JSON file."""
    global PATHWAY_CONFIG, PATHWAY_NAME, MODEL_PATH, VIZ_OUTPUT_DIR

    with open(config_path, "r") as f:
        PATHWAY_CONFIG = json.load(f)

    # Extract pathway name from config
    PATHWAY_NAME = PATHWAY_CONFIG.get("pathway_name", "pathway")

    # Extract model path from config and resolve it
    output_model = PATHWAY_CONFIG.get("model", {}).get("output")
    if output_model:
        output_path = Path(output_model)
        # If relative path, resolve relative to project root (parent of implement_pathway)
        if not output_path.is_absolute():
            project_root = Path(__file__).parent.parent
            MODEL_PATH = str(project_root / output_path)
        else:
            MODEL_PATH = str(output_path)

    # Extract visualization output directory from config
    viz_config = PATHWAY_CONFIG.get("visualizations", {})
    if "output_directory" in viz_config:
        VIZ_OUTPUT_DIR = viz_config["output_directory"]

    return PATHWAY_CONFIG


def get_pathway_reaction_ids(config):
    """Extract all reaction IDs from the pathway configuration."""
    reaction_ids = set()

    if "reactions" in config:
        for rxn in config["reactions"]:
            if "id" in rxn:
                reaction_ids.add(rxn["id"])

    return reaction_ids


# Load default configuration if available
try:
    CONFIG = load_config()
    BASE_MODEL_PATH, OUTPUT_MODEL_PATH = get_model_paths(CONFIG)

    # Get visualization settings from config
    viz_config = CONFIG.get("visualizations", {})
    VIZ_OUTPUT_DIR = viz_config.get("output_directory", "figures")

except Exception:
    # Silently use defaults - will be overridden by --config if provided
    BASE_MODEL_PATH = "../models/endoA_250917_3.json"
    OUTPUT_MODEL_PATH = "../models/endoA_250917_3_with_cytoskeleton.json"
    VIZ_OUTPUT_DIR = "figures"

# ============================================================================
# VISUALIZATION 1: D3.js Interactive HTML
# ============================================================================


def visualize_d3js(model_path=None, output_path=None, pathway_config=None):
    """Generate D3.js interactive HTML visualization."""
    global PATHWAY_NAME

    if model_path is None:
        model_path = MODEL_PATH if MODEL_PATH else OUTPUT_MODEL_PATH
    if output_path is None:
        output_path = os.path.join(VIZ_OUTPUT_DIR, f"{PATHWAY_NAME}_network_d3.html")

    print("\n" + "=" * 80)
    print("VISUALIZATION 1: D3.js Interactive HTML")
    print("=" * 80)

    def load_model(model_path):
        """Load the model JSON file."""
        print(f"Loading model from: {model_path}")
        with open(model_path, "r") as f:
            model = json.load(f)
        return model

    def extract_pathway_network(model, pathway_config):
        """Extract pathway reactions and metabolites."""
        print(f"Extracting {PATHWAY_NAME} network...")

        # Get reaction IDs from config if available
        if pathway_config:
            target_reaction_ids = get_pathway_reaction_ids(pathway_config)
            print(f"  Target reaction IDs: {len(target_reaction_ids)}")
        else:
            # Fallback: Use cytoskeleton compartment
            target_reaction_ids = None

        # Find matching reactions
        pathway_reactions = []
        for rxn in model["reactions"]:
            if target_reaction_ids:
                # Use config-based filtering
                if rxn["id"] in target_reaction_ids:
                    pathway_reactions.append(rxn)
            else:
                # Fallback: filter by 'ck' compartment
                metabolites = rxn.get("metabolites", {})
                has_ck = any(met_id.endswith("ck") for met_id in metabolites.keys())
                if has_ck:
                    pathway_reactions.append(rxn)

        print(f"  Found {len(pathway_reactions)} pathway reactions")

        # Find all involved metabolites
        metabolite_ids = set()
        for rxn in pathway_reactions:
            metabolite_ids.update(rxn.get("metabolites", {}).keys())

        # Get metabolite details
        metabolite_map = {}
        for met in model["metabolites"]:
            if met["id"] in metabolite_ids:
                metabolite_map[met["id"]] = {
                    "id": met["id"],
                    "name": met.get("name", "Unknown"),
                    "formula": met.get("formula", ""),
                    "compartment": met.get("compartment", ""),
                    "charge": met.get("charge", 0),
                }

        print(f"  Found {len(metabolite_map)} metabolites")

        # Categorize reactions by subsystem
        subsystems = defaultdict(list)
        for rxn in pathway_reactions:
            subsystem = rxn.get("subsystem", "Other")
            subsystems[subsystem].append(rxn)

        return pathway_reactions, metabolite_map, subsystems

    def create_network_data(reactions, metabolite_map, subsystems):
        """Create nodes and edges for network visualization."""
        print("Creating network data structure...")

        nodes = []
        edges = []

        # Color scheme for compartments
        compartment_colors = {
            "c": "#3498db",  # Cytosol - Blue
            "r": "#e74c3c",  # ER - Red
            "g": "#2ecc71",  # Golgi - Green
            "v": "#9b59b6",  # Vesicle - Purple
            "m": "#f39c12",  # Mitochondria - Orange
            "gl": "#1abc9c",  # Glycocalyx - Turquoise
            "e": "#95a5a6",  # Extracellular - Gray
            "ck": "#e67e22",  # Cytoskeleton - Dark Orange
            "l": "#34495e",  # Lysosome - Dark Gray
            "n": "#16a085",  # Nucleus - Dark Turquoise
            "x": "#c0392b",  # Peroxisome - Dark Red
        }

        # Compartment full names for legend
        compartment_names = {
            "c": "Cytosol",
            "r": "Endoplasmic Reticulum",
            "g": "Golgi Apparatus",
            "v": "Vesicle",
            "m": "Mitochondria",
            "gl": "Glycocalyx",
            "e": "Extracellular",
            "ck": "Cytoskeleton",
            "l": "Lysosome",
            "n": "Nucleus",
            "x": "Peroxisome",
        }

        # Color scheme for subsystems (for reactions)
        subsystem_colors = {
            "Cytoskeleton transport": "#3498db",
            "Cytoskeleton synthesis": "#e74c3c",
            "Cytoskeleton regulation": "#f39c12",
            "Cytoskeleton dynamics": "#9b59b6",
            "Cytoskeleton turnover": "#1abc9c",
            "Other": "#95a5a6",
        }

        # Add metabolite nodes
        metabolite_node_map = {}
        for met_id, met_info in metabolite_map.items():
            node_id = f"met_{met_id}"

            compartment = met_info["compartment"]
            color = compartment_colors.get(compartment, "#95a5a6")

            # Base size
            size = 12

            # Increase size for key metabolites
            key_metabolites = [
                "ATP",
                "ADP",
                "GTP",
                "GDP",
                "H2O",
                "Pi",
                "H+",
                "Phosphate",
            ]
            if any(key in met_info["name"] for key in key_metabolites):
                size = 18

            # Increase size for proteins and large molecules
            if "protein" in met_info["name"].lower() or met_info["id"].startswith(
                "MAM100"
            ):
                size = 25

            nodes.append(
                {
                    "id": node_id,
                    "label": met_info["name"][:30],
                    "full_name": met_info["name"],
                    "formula": met_info["formula"],
                    "compartment": compartment,
                    "compartment_name": compartment_names.get(compartment, compartment),
                    "charge": met_info["charge"],
                    "type": "metabolite",
                    "color": color,
                    "size": size,
                }
            )

            metabolite_node_map[met_id] = node_id

        # Add reaction nodes
        reaction_node_map = {}
        for rxn in reactions:
            node_id = f"rxn_{rxn['id']}"
            subsystem = rxn.get("subsystem", "Other")

            nodes.append(
                {
                    "id": node_id,
                    "label": "",
                    "full_name": rxn.get("name", rxn["id"]),
                    "subsystem": subsystem,
                    "type": "reaction",
                    "color": subsystem_colors.get(subsystem, "#95a5a6"),
                    "size": 8,
                    "shape": "square",
                }
            )

            reaction_node_map[rxn["id"]] = node_id

        # Add edges
        for rxn in reactions:
            rxn_node = reaction_node_map[rxn["id"]]
            reversible = rxn.get("lower_bound", 0) < 0

            for met_id, coeff in rxn.get("metabolites", {}).items():
                if met_id not in metabolite_node_map:
                    continue

                met_node = metabolite_node_map[met_id]

                if coeff < 0:  # Substrate
                    edges.append(
                        {
                            "source": met_node,
                            "target": rxn_node,
                            "stoichiometry": abs(coeff),
                            "type": "substrate",
                            "reversible": reversible,
                        }
                    )
                else:  # Product
                    edges.append(
                        {
                            "source": rxn_node,
                            "target": met_node,
                            "stoichiometry": abs(coeff),
                            "type": "product",
                            "reversible": reversible,
                        }
                    )

        print(f"  Created {len(nodes)} nodes and {len(edges)} edges")

        return nodes, edges, subsystem_colors, compartment_colors, compartment_names

    # Generate HTML content (using the D3.js template from original script)
    # ... [HTML content omitted for brevity - same as original visualize_cytoskeleton_network.py]

    def generate_html(
        nodes,
        edges,
        subsystem_colors,
        compartment_colors,
        compartment_names,
        output_path,
    ):
        """Generate interactive HTML visualization using D3.js."""
        print(f"Generating HTML visualization: {output_path}")

        nodes_json = json.dumps(nodes, indent=2)
        edges_json = json.dumps(edges, indent=2)
        colors_json = json.dumps(subsystem_colors, indent=2)
        comp_colors_json = json.dumps(compartment_colors, indent=2)
        comp_names_json = json.dumps(compartment_names, indent=2)

        # Get unique compartments present in the data
        unique_compartments = sorted(
            set(n["compartment"] for n in nodes if n["type"] == "metabolite")
        )

        # Generate legend HTML for compartments
        legend_html = ""
        for comp in unique_compartments:
            color = compartment_colors.get(comp, "#95a5a6")
            name = compartment_names.get(comp, comp)
            legend_html += f'            <div class="legend-item"><span class="legend-color" style="background: {color}; border-radius: 50%;"></span>{name} [{comp}]</div>\n'

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{PATHWAY_NAME} Metabolic Network</title>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <style>
        body {{ margin: 0; padding: 10px; font-family: 'Arial', sans-serif; background: #f5f5f5; }}
        #header {{ background: white; padding: 10px 15px; margin-bottom: 10px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        h1 {{ margin: 0 0 5px 0; color: #2c3e50; font-size: 24px; }}
        button {{ padding: 6px 12px; margin-right: 8px; background: #3498db; color: white; border: none; border-radius: 4px; cursor: pointer; font-size: 13px; }}
        button:hover {{ background: #2980b9; }}
        #legend {{ background: white; padding: 10px 15px; margin-bottom: 10px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .legend-item {{ display: inline-block; margin-right: 20px; margin-bottom: 5px; }}
        .legend-color {{ display: inline-block; width: 20px; height: 20px; margin-right: 5px; vertical-align: middle; border-radius: 3px; }}
        #visualization {{ background: white; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); width: 100%; overflow: hidden; }}
        #tooltip {{ position: absolute; padding: 12px; background: rgba(0, 0, 0, 0.9); color: white; border-radius: 4px; 
                    pointer-events: none; opacity: 0; transition: opacity 0.3s; font-size: 12px; max-width: 300px; z-index: 1000; }}
        .link {{ fill: none; stroke: #999; stroke-opacity: 0.6; stroke-width: 2px; }}
        .node {{ cursor: pointer; stroke: #fff; stroke-width: 2px; }}
        .node:hover {{ stroke: #000; stroke-width: 3px; }}
        .node-label {{ font-size: 10px; pointer-events: none; text-anchor: middle; fill: #2c3e50; }}
        .stats {{ background: #ecf0f1; padding: 8px; border-radius: 4px; margin-top: 8px; font-size: 13px; display: inline-block; }}
        p {{ margin: 5px 0; font-size: 14px; }}
    </style>
</head>
<body>
    <div id="header">
        <h1>🧬 {PATHWAY_NAME} Metabolic Network</h1>
        <p>Interactive visualization of {PATHWAY_NAME.lower().replace('_', ' ')} pathway</p>
        <div class="stats">
            <strong>Network Statistics:</strong>
            <span id="stats-metabolites"></span> metabolites, 
            <span id="stats-reactions"></span> reactions, 
            <span id="stats-edges"></span> connections
        </div>
        <div id="controls">
            <button onclick="resetZoom()">Reset View</button>
            <button onclick="toggleLabels()">Toggle Labels</button>
            <button onclick="colorByCompartment()">Color by Compartment</button>
        </div>
    </div>
    
    <div id="legend">
        <strong>Compartments:</strong><br>
        <div style="margin: 10px 0;">
{legend_html}
            <div class="legend-item"><span class="legend-color" style="background: #666; width: 16px; height: 16px;"></span>Reactions</div>
        </div>
    </div>
    
    <div id="visualization"></div>
    <div id="tooltip"></div>
    
    <script>
        const nodes = {nodes_json};
        const edges = {edges_json};
        const subsystemColors = {colors_json};
        
        const metaboliteCount = nodes.filter(n => n.type !== 'reaction').length;
        const reactionCount = nodes.filter(n => n.type === 'reaction').length;
        document.getElementById('stats-metabolites').textContent = metaboliteCount;
        document.getElementById('stats-reactions').textContent = reactionCount;
        document.getElementById('stats-edges').textContent = edges.length;
        
        const width = window.innerWidth - 40;
        const height = window.innerHeight - 200;
        
        const svg = d3.select("#visualization")
            .append("svg")
            .attr("width", width)
            .attr("height", height);
        
        // Define arrow markers for directed edges
        svg.append("defs").selectAll("marker")
            .data(["arrow"])
            .join("marker")
            .attr("id", "arrow")
            .attr("viewBox", "0 -5 10 10")
            .attr("refX", 20)
            .attr("refY", 0)
            .attr("markerWidth", 6)
            .attr("markerHeight", 6)
            .attr("orient", "auto")
            .append("path")
            .attr("d", "M0,-5L10,0L0,5")
            .attr("fill", "#999");
        
        const g = svg.append("g");
        
        const zoom = d3.zoom()
            .scaleExtent([0.1, 10])
            .on("zoom", (event) => {{ g.attr("transform", event.transform); }});
        
        svg.call(zoom);
        
        const tooltip = d3.select("#tooltip");
        
        // Define compartment centers for spatial clustering
        const compartmentCenters = {{
            'c': {{ x: width * 0.3, y: height * 0.5 }},      // Cytosol - Left center
            'r': {{ x: width * 0.5, y: height * 0.3 }},      // ER - Top center
            'g': {{ x: width * 0.7, y: height * 0.3 }},      // Golgi - Top right
            'v': {{ x: width * 0.7, y: height * 0.5 }},      // Vesicle - Right center
            'm': {{ x: width * 0.3, y: height * 0.7 }},      // Mitochondria - Bottom left
            'gl': {{ x: width * 0.5, y: height * 0.7 }},     // Glycocalyx - Bottom center
            'e': {{ x: width * 0.7, y: height * 0.7 }},      // Extracellular - Bottom right
            'ck': {{ x: width * 0.5, y: height * 0.5 }},     // Cytoskeleton - Center
        }};
        
        // Custom force to cluster nodes by compartment
        function forceCompartment(alpha) {{
            for (let node of nodes) {{
                if (node.type === 'metabolite' && node.compartment) {{
                    const center = compartmentCenters[node.compartment];
                    if (center) {{
                        node.vx += (center.x - node.x) * alpha * 0.1;
                        node.vy += (center.y - node.y) * alpha * 0.1;
                    }}
                }}
            }}
        }}
        
        const simulation = d3.forceSimulation(nodes)
            .force("link", d3.forceLink(edges).id(d => d.id).distance(80))
            .force("charge", d3.forceManyBody().strength(d => d.type === 'reaction' ? -300 : -500))
            .force("center", d3.forceCenter(width / 2, height / 2))
            .force("collision", d3.forceCollide().radius(d => d.size + 5))
            .force("compartment", forceCompartment);
        
        const link = g.append("g")
            .selectAll("line")
            .data(edges)
            .join("line")
            .attr("class", "link")
            .attr("marker-end", "url(#arrow)")
            .attr("stroke-width", d => Math.sqrt(d.stoichiometry) * 1.5);
        
        const node = g.append("g")
            .selectAll("circle,rect")
            .data(nodes)
            .join(enter => {{
                return enter.append(d => {{
                    if (d.shape === 'square') {{
                        return document.createElementNS("http://www.w3.org/2000/svg", "rect");
                    }}
                    return document.createElementNS("http://www.w3.org/2000/svg", "circle");
                }});
            }})
            .attr("class", "node")
            .attr("fill", d => d.color)
            .each(function(d) {{
                if (d.shape === 'square') {{
                    d3.select(this).attr("width", d.size * 2).attr("height", d.size * 2);
                }} else {{
                    d3.select(this).attr("r", d.size);
                }}
            }})
            .call(drag(simulation))
            .on("mouseover", showTooltip)
            .on("mouseout", hideTooltip);
        
        let showLabels = true;
        const labels = g.append("g")
            .selectAll("text")
            .data(nodes.filter(d => d.type !== 'reaction' && d.size >= 15))
            .join("text")
            .attr("class", "node-label")
            .text(d => d.label);
        
        simulation.on("tick", () => {{
            link
                .attr("x1", d => d.source.x)
                .attr("y1", d => d.source.y)
                .attr("x2", d => d.target.x)
                .attr("y2", d => d.target.y);
            
            node.each(function(d) {{
                if (d.shape === 'square') {{
                    d3.select(this).attr("x", d.x - d.size).attr("y", d.y - d.size);
                }} else {{
                    d3.select(this).attr("cx", d.x).attr("cy", d.y);
                }}
            }});
            
            labels.attr("x", d => d.x).attr("y", d => d.y + d.size + 12);
        }});
        
        function drag(simulation) {{
            function dragstarted(event, d) {{
                if (!event.active) simulation.alphaTarget(0.3).restart();
                d.fx = d.x; d.fy = d.y;
            }}
            function dragged(event, d) {{ d.fx = event.x; d.fy = event.y; }}
            function dragended(event, d) {{
                if (!event.active) simulation.alphaTarget(0);
                d.fx = null; d.fy = null;
            }}
            return d3.drag()
                .on("start", dragstarted)
                .on("drag", dragged)
                .on("end", dragended);
        }}
        
        function showTooltip(event, d) {{
            let content = `<strong>${{d.full_name || d.label}}</strong><br>`;
            if (d.type === 'reaction') {{
                content += `Type: Reaction<br>Subsystem: ${{d.subsystem}}`;
            }} else {{
                content += `Type: Metabolite<br>`;
                content += `Compartment: ${{d.compartment_name || d.compartment}}<br>`;
                content += `Formula: ${{d.formula}}<br>`;
                content += `Charge: ${{d.charge}}`;
            }}
            tooltip.style("opacity", 1)
                .style("left", (event.pageX + 10) + "px")
                .style("top", (event.pageY - 10) + "px")
                .html(content);
        }}
        
        function hideTooltip() {{ tooltip.style("opacity", 0); }}
        function resetZoom() {{ svg.transition().duration(750).call(zoom.transform, d3.zoomIdentity); }}
        function toggleLabels() {{ showLabels = !showLabels; labels.style("opacity", showLabels ? 1 : 0); }}
        function colorByCompartment() {{
            node.transition().duration(500).attr("fill", d => {{
                if (d.type === 'protein') return '#9b59b6';
                if (d.compartment === 'ck') return '#e74c3c';
                if (d.compartment === 'c') return '#3498db';
                if (d.type === 'reaction') return subsystemColors[d.subsystem] || '#95a5a6';
                return '#95a5a6';
            }});
        }}
    </script>
</body>
</html>
"""

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "w") as f:
            f.write(html_content)

        print(f"✓ Visualization saved!")

    model = load_model(model_path)
    reactions, metabolite_map, subsystems = extract_pathway_network(
        model, pathway_config
    )
    nodes, edges, subsystem_colors, compartment_colors, compartment_names = (
        create_network_data(reactions, metabolite_map, subsystems)
    )
    generate_html(
        nodes,
        edges,
        subsystem_colors,
        compartment_colors,
        compartment_names,
        output_path,
    )

    print(f"✓ D3.js visualization saved: {output_path}")
    print(f"  Open in browser: file://{os.path.abspath(output_path)}")


# ============================================================================
# VISUALIZATION 2: Plotly Interactive Network
# ============================================================================


def visualize_plotly(model_path=None, output_path=None, pathway_config=None):
    """Generate Plotly interactive network visualization with FVA."""
    global PATHWAY_NAME

    if model_path is None:
        model_path = MODEL_PATH if MODEL_PATH else OUTPUT_MODEL_PATH
    if output_path is None:
        output_path = os.path.join(
            VIZ_OUTPUT_DIR, f"{PATHWAY_NAME}_network_plotly.html"
        )

    print("\n" + "=" * 80)
    print("VISUALIZATION 2: Plotly Interactive Network")
    print("=" * 80)

    try:
        import cobra
        import plotly.graph_objects as go
        import networkx as nx
    except ImportError as e:
        print(f"✗ Error: Required package not installed: {e}")
        print("  Install with: pip install cobra plotly networkx")
        return

    # Load model
    print(f"Loading model from: {model_path}")
    model = cobra.io.load_json_model(model_path)

    # Run FBA and FVA
    print("Running FBA...")
    solution = model.optimize()
    optimal_biomass = solution.objective_value

    print("Getting pathway reactions...")

    # Get reaction IDs from config if available
    if pathway_config:
        target_reaction_ids = get_pathway_reaction_ids(pathway_config)
        ck_reactions = [rxn for rxn in model.reactions if rxn.id in target_reaction_ids]
    else:
        # Fallback: use cytoskeleton compartment
        ck_reactions = [
            rxn
            for rxn in model.reactions
            if any(met.compartment == "ck" for met in rxn.metabolites)
        ]

    print("Running FVA analysis...")
    fva_result = cobra.flux_analysis.flux_variability_analysis(
        model, reaction_list=[rxn.id for rxn in ck_reactions], fraction_of_optimum=0.9
    )
    fva_result["range"] = fva_result["maximum"] - fva_result["minimum"]

    # Build network graph
    print("Building network graph...")
    G = nx.Graph()

    # Add metabolite nodes (from all compartments involved in pathway reactions)
    metabolites = set()
    for rxn in ck_reactions:
        for met in rxn.metabolites:
            metabolites.add(met)
            G.add_node(
                met.id,
                type="metabolite",
                name=met.name,
                formula=met.formula if met.formula else "N/A",
                compartment=met.compartment,
            )

    # Add reaction nodes and edges
    for rxn in ck_reactions:
        is_flexible = fva_result.loc[rxn.id, "range"] > 1e-6
        flux_range = fva_result.loc[rxn.id, "range"]
        min_flux = fva_result.loc[rxn.id, "minimum"]
        max_flux = fva_result.loc[rxn.id, "maximum"]

        G.add_node(
            rxn.id,
            type="reaction",
            name=rxn.name,
            flexible=is_flexible,
            flux_range=flux_range,
            min_flux=min_flux,
            max_flux=max_flux,
        )

        for met, coeff in rxn.metabolites.items():
            if coeff < 0:  # Substrate
                G.add_edge(met.id, rxn.id, type="substrate", coeff=abs(coeff))
            elif coeff > 0:  # Product
                G.add_edge(rxn.id, met.id, type="product", coeff=coeff)

    # Calculate layout
    print("Calculating node positions...")
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # Create traces for edges and nodes (same as original script)
    # ... [Plotly traces code omitted for brevity]

    # Separate nodes
    met_nodes = [node for node in G.nodes() if G.nodes[node]["type"] == "metabolite"]
    rxn_nodes = [node for node in G.nodes() if G.nodes[node]["type"] == "reaction"]
    flexible_rxn_nodes = [n for n in rxn_nodes if G.nodes[n]["flexible"]]
    blocked_rxn_nodes = [n for n in rxn_nodes if not G.nodes[n]["flexible"]]

    # Create edge traces
    edge_traces = []

    substrate_x, substrate_y = [], []
    for edge in G.edges():
        if G.edges[edge].get("type") == "substrate":
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            substrate_x.extend([x0, x1, None])
            substrate_y.extend([y0, y1, None])

    edge_trace_substrate = go.Scatter(
        x=substrate_x,
        y=substrate_y,
        line=dict(width=1, color="#888"),
        hoverinfo="none",
        mode="lines",
        name="Substrate",
        showlegend=True,
    )
    edge_traces.append(edge_trace_substrate)

    product_x, product_y = [], []
    for edge in G.edges():
        if G.edges[edge].get("type") == "product":
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            product_x.extend([x0, x1, None])
            product_y.extend([y0, y1, None])

    edge_trace_product = go.Scatter(
        x=product_x,
        y=product_y,
        line=dict(width=1, color="#BBB"),
        hoverinfo="none",
        mode="lines",
        name="Product",
        showlegend=True,
    )
    edge_traces.append(edge_trace_product)

    # Metabolite trace
    met_x, met_y, met_text, met_hover = [], [], [], []
    for node in met_nodes:
        x, y = pos[node]
        met_x.append(x)
        met_y.append(y)
        node_data = G.nodes[node]
        met_text.append(node_data["name"][:20])
        met_hover.append(
            f"<b>{node_data['name']}</b><br>"
            + f"ID: {node}<br>"
            + f"Formula: {node_data['formula']}<br>"
            + f"Type: Metabolite"
        )

    met_trace = go.Scatter(
        x=met_x,
        y=met_y,
        mode="markers+text",
        name="Metabolites",
        text=met_text,
        textposition="top center",
        textfont=dict(size=8),
        hovertext=met_hover,
        hoverinfo="text",
        marker=dict(
            size=15,
            color="lightblue",
            symbol="circle",
            line=dict(width=2, color="darkblue"),
        ),
        showlegend=True,
    )

    # Flexible reactions trace
    flex_x, flex_y, flex_text, flex_hover, flex_colors = [], [], [], [], []
    for node in flexible_rxn_nodes:
        x, y = pos[node]
        flex_x.append(x)
        flex_y.append(y)
        node_data = G.nodes[node]

        label = (
            node_data["name"]
            .replace("Transport of ", "")
            .replace(" from cytosol to cytoskeleton", "")
        )
        label = label.replace("Export of ", "").replace(
            " from cytoskeleton to cytosol", ""
        )
        flex_text.append(label[:25])

        flex_hover.append(
            f"<b>{node_data['name']}</b><br>"
            + f"ID: {node}<br>"
            + f"Status: FLEXIBLE<br>"
            + f"Min flux: {node_data['min_flux']:.2f}<br>"
            + f"Max flux: {node_data['max_flux']:.2f}<br>"
            + f"Range: {node_data['flux_range']:.2f}"
        )

        flex_colors.append(node_data["flux_range"])

    flex_trace = go.Scatter(
        x=flex_x,
        y=flex_y,
        mode="markers+text",
        name="Flexible Reactions",
        text=flex_text,
        textposition="bottom center",
        textfont=dict(size=7, color="darkred"),
        hovertext=flex_hover,
        hoverinfo="text",
        marker=dict(
            size=20,
            color=flex_colors,
            colorscale="Reds",
            symbol="square",
            line=dict(width=2, color="darkred"),
            colorbar=dict(title="Flux Range", thickness=15, len=0.5, x=1.1),
            showscale=True,
        ),
        showlegend=True,
    )

    # Blocked reactions trace
    blocked_x, blocked_y, blocked_text, blocked_hover = [], [], [], []
    for node in blocked_rxn_nodes:
        x, y = pos[node]
        blocked_x.append(x)
        blocked_y.append(y)
        node_data = G.nodes[node]

        label = (
            node_data["name"]
            .replace("Transport of ", "")
            .replace(" from cytosol to cytoskeleton", "")
        )
        label = label.replace("Export of ", "").replace(
            " from cytoskeleton to cytosol", ""
        )
        blocked_text.append(label[:25])

        blocked_hover.append(
            f"<b>{node_data['name']}</b><br>"
            + f"ID: {node}<br>"
            + f"Status: BLOCKED<br>"
            + f"Min flux: {node_data['min_flux']:.2f}<br>"
            + f"Max flux: {node_data['max_flux']:.2f}<br>"
            + f"Range: {node_data['flux_range']:.2f}"
        )

    blocked_trace = go.Scatter(
        x=blocked_x,
        y=blocked_y,
        mode="markers+text",
        name="Blocked Reactions",
        text=blocked_text,
        textposition="bottom center",
        textfont=dict(size=7, color="gray"),
        hovertext=blocked_hover,
        hoverinfo="text",
        marker=dict(
            size=20,
            color="lightgray",
            symbol="square",
            line=dict(width=2, color="gray"),
        ),
        showlegend=True,
    )

    # Create figure
    print("Creating interactive figure...")
    fig = go.Figure(data=edge_traces + [met_trace, flex_trace, blocked_trace])

    fig.update_layout(
        title=dict(
            text=f"<b>{PATHWAY_NAME} Metabolic Network</b><br>"
            + f"<sub>FVA at 90% optimal biomass ({0.9*optimal_biomass:.2f})</sub><br>"
            + f"<sub>{len(flexible_rxn_nodes)} flexible reactions (red) | {len(blocked_rxn_nodes)} blocked reactions (gray)</sub>",
            x=0.5,
            xanchor="center",
        ),
        showlegend=True,
        hovermode="closest",
        margin=dict(b=20, l=5, r=5, t=100),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor="white",
        width=1400,
        height=900,
    )

    # Save
    fig.write_html(output_path)
    print(f"✓ Plotly visualization saved: {output_path}")
    print(f"\n  Network Statistics:")
    print(f"    Total metabolites: {len(metabolites)}")
    print(f"    Total reactions: {len(ck_reactions)}")
    print(
        f"    Flexible reactions: {len(flexible_rxn_nodes)} ({100*len(flexible_rxn_nodes)/len(ck_reactions):.1f}%)"
    )
    print(
        f"    Blocked reactions: {len(blocked_rxn_nodes)} ({100*len(blocked_rxn_nodes)/len(ck_reactions):.1f}%)"
    )


# ============================================================================
# VISUALIZATION 3: GraphViz Static Image
# ============================================================================


def visualize_graphviz(model_path=None, output_prefix=None, pathway_config=None):
    """Generate GraphViz static visualization (PNG/PDF/SVG)."""
    global PATHWAY_NAME

    if model_path is None:
        model_path = MODEL_PATH if MODEL_PATH else OUTPUT_MODEL_PATH
    if output_prefix is None:
        output_prefix = os.path.join(VIZ_OUTPUT_DIR, f"{PATHWAY_NAME}_network")

    print("\n" + "=" * 80)
    print("VISUALIZATION 3: GraphViz Static Image")
    print("=" * 80)

    try:
        import cobra
        from graphviz import Digraph
    except ImportError as e:
        print(f"✗ Error: Required package not installed: {e}")
        print("  Install with: pip install cobra graphviz")
        return

    # Load model
    print(f"Loading model from: {model_path}")
    model = cobra.io.load_json_model(model_path)

    # Run FBA and FVA
    print("Running FBA...")
    solution = model.optimize()

    # Get reaction IDs from config if available
    if pathway_config:
        target_reaction_ids = get_pathway_reaction_ids(pathway_config)
        ck_reactions = [rxn for rxn in model.reactions if rxn.id in target_reaction_ids]
    else:
        # Fallback: use cytoskeleton compartment
        ck_reactions = [
            rxn
            for rxn in model.reactions
            if any(met.compartment == "ck" for met in rxn.metabolites)
        ]

    print("Running FVA analysis...")
    fva_result = cobra.flux_analysis.flux_variability_analysis(
        model, reaction_list=[rxn.id for rxn in ck_reactions], fraction_of_optimum=0.9
    )
    fva_result["range"] = fva_result["maximum"] - fva_result["minimum"]

    # Create digraph
    dot = Digraph(comment=f"{PATHWAY_NAME} Metabolic Network")
    # Left-to-right layout for widescreen monitors, compact spacing
    dot.attr(
        rankdir="LR",
        size="24,12",
        ratio="compress",
        dpi="300",
        ranksep="0.5",
        nodesep="0.3",
    )
    dot.attr(
        "node",
        shape="box",
        style="rounded,filled",
        fontname="Arial",
        fontsize="12",
        width="1.2",
        height="0.5",
    )

    # Metabolite categories
    energy_mets = ["atp_ck", "adp_ck", "gtp_ck", "gdp_ck", "pi_ck", "h2o_ck", "h_ck"]
    amino_acids = [
        "ala__L_ck",
        "arg__L_ck",
        "asn__L_ck",
        "asp__L_ck",
        "cys__L_ck",
        "gln__L_ck",
        "glu__L_ck",
        "gly_ck",
        "his__L_ck",
        "ile__L_ck",
        "leu__L_ck",
        "lys__L_ck",
        "met__L_ck",
        "phe__L_ck",
        "pro__L_ck",
        "ser__L_ck",
        "thr__L_ck",
        "trp__L_ck",
        "tyr__L_ck",
        "val__L_ck",
    ]
    proteins = [
        "actin_protein_ck",
        "tubulin_dimer_ck",
        "vimentin_protein_ck",
        "myosin_light_chain_ck",
        "actin_complex_ck",
        "vimentin_glcnac_ck",
        "myosin_p_ck",
    ]
    modifications = ["udpacgal_ck", "udp_ck", "pail4p_ck", "pail45p_ck"]

    # Categorize reactions
    reaction_categories = {}
    for rxn in ck_reactions:
        name_lower = rxn.name.lower()
        rxn_id_lower = rxn.id.lower()
        is_flexible = fva_result.loc[rxn.id, "range"] > 1e-6

        # Check for specific patterns - order matters!
        if "demand" in name_lower or rxn.id.startswith("DM_"):
            reaction_categories[rxn.id] = ("Demand", is_flexible)
        elif "synthesis" in name_lower or "_synth" in rxn_id_lower:
            reaction_categories[rxn.id] = ("Synthesis", is_flexible)
        elif "export" in name_lower:
            reaction_categories[rxn.id] = ("Export", is_flexible)
        elif "transport" in name_lower or rxn.id.startswith("T_"):
            # Subcategorize transport reactions
            if any(
                aa in name_lower
                for aa in [
                    "alanine",
                    "arginine",
                    "asparagine",
                    "aspartate",
                    "cysteine",
                    "glutamine",
                    "glutamate",
                    "glycine",
                    "histidine",
                    "isoleucine",
                    "leucine",
                    "lysine",
                    "methionine",
                    "phenylalanine",
                    "proline",
                    "serine",
                    "threonine",
                    "tryptophan",
                    "tyrosine",
                    "valine",
                ]
            ):
                reaction_categories[rxn.id] = ("AA Transport", is_flexible)
            elif "atp" in name_lower or "gtp" in name_lower:
                reaction_categories[rxn.id] = ("Energy Transport", is_flexible)
            elif (
                "h2o" in name_lower
                or "h+" in name_lower
                or "phosphate" in name_lower
                or "adp" in name_lower
                or "gdp" in name_lower
                or "udp" in name_lower
            ):
                reaction_categories[rxn.id] = ("Small Molecule Transport", is_flexible)
            else:
                reaction_categories[rxn.id] = ("Other Transport", is_flexible)
        elif any(kw in name_lower for kw in ["hydrolysis", "regeneration"]):
            reaction_categories[rxn.id] = ("Energy", is_flexible)
        elif any(
            kw in name_lower
            for kw in [
                "remodel",
                "phosphorylation",
                "dephosphorylation",
                "modification",
                "glycosylation",
                "glyco",
            ]
        ):
            reaction_categories[rxn.id] = ("Modification", is_flexible)
        else:
            reaction_categories[rxn.id] = ("Other", is_flexible)

    # Add metabolite nodes
    added_mets = set()
    for rxn in ck_reactions:
        for met in rxn.metabolites:
            if met.id in added_mets:
                continue

            # Color by type
            if met.id in energy_mets:
                color = "#FFE6E6"
                label = met.name.replace(" [cytoskeleton]", "")
            elif met.id in amino_acids:
                color = "#E6F3FF"
                label = met.name.replace(" [cytoskeleton]", "").replace("L-", "")
            elif met.id in proteins:
                color = "#E6FFE6"
                label = met.name.replace(" [cytoskeleton]", "")
            elif met.id in modifications:
                color = "#FFF5E6"
                label = met.name.replace(" [cytoskeleton]", "")
            else:
                color = "#F0F0F0"
                label = met.name.replace(" [cytoskeleton]", "")

            dot.node(met.id, label, fillcolor=color, shape="ellipse", fontsize="11")
            added_mets.add(met.id)

    # Add reaction nodes and edges
    for rxn in ck_reactions:
        cat, is_flexible = reaction_categories[rxn.id]

        # Color by category - distinct colors for each
        if cat == "Synthesis":
            color = "#FF6B6B"  # Red
        elif cat == "Demand":
            color = "#9B59B6"  # Purple
        elif cat == "AA Transport":
            color = "#3498DB"  # Blue
        elif cat == "Energy Transport":
            color = "#E67E22"  # Orange
        elif cat == "Small Molecule Transport":
            color = "#1ABC9C"  # Teal
        elif cat == "Other Transport":
            color = "#16A085"  # Dark Teal
        elif cat == "Energy":
            color = "#F1C40F"  # Yellow
        elif cat == "Modification":
            color = "#E91E63"  # Pink
        elif cat == "Export":
            color = "#95A5A6"  # Gray
        else:  # Other
            color = "#BDC3C7"  # Light Gray

        # Style by flexibility (bold=flexible, dashed=blocked)
        if is_flexible:
            style = "filled,bold"
        else:
            style = "filled,dashed"

        # Simplified label
        label = rxn.name.replace("Transport of ", "").replace(
            " from cytosol to cytoskeleton", ""
        )
        label = label.replace("Export of ", "").replace(
            " from cytoskeleton to cytosol", ""
        )
        label = label.replace(" [cytoskeleton]", "")
        if len(label) > 40:
            label = label[:37] + "..."

        dot.node(
            rxn.id, label, fillcolor=color, shape="box", style=style, fontsize="11"
        )

        # Add edges with larger labels
        for met, coeff in rxn.metabolites.items():
            if coeff < 0:  # Substrate
                dot.edge(
                    met.id,
                    rxn.id,
                    label=f"{abs(coeff):.1f}" if abs(coeff) != 1.0 else "",
                    fontsize="10",
                    len="1.0",
                )
            elif coeff > 0:  # Product
                dot.edge(
                    rxn.id,
                    met.id,
                    label=f"{coeff:.1f}" if coeff != 1.0 else "",
                    fontsize="10",
                    len="1.0",
                )

    # Add legend - dynamically include only categories present in the network
    with dot.subgraph(name="cluster_legend") as legend:
        legend.attr(label="Legend", fontsize="12", style="filled", color="lightgrey")

        # Metabolites
        legend.node(
            "leg_title", "Metabolites:", shape="plaintext", fontsize="10", style=""
        )
        legend.node(
            "leg_met",
            "Metabolite",
            shape="ellipse",
            fillcolor="#E6F3FF",
            style="filled",
            fontsize="10",
        )

        # Get unique categories present in the network
        categories_present = sorted(set(cat for cat, _ in reaction_categories.values()))

        # Define colors for all possible categories
        category_colors = {
            "Synthesis": "#FF6B6B",
            "Demand": "#9B59B6",
            "AA Transport": "#3498DB",
            "Energy Transport": "#E67E22",
            "Small Molecule Transport": "#1ABC9C",
            "Other Transport": "#16A085",
            "Energy": "#F1C40F",
            "Modification": "#E91E63",
            "Export": "#95A5A6",
            "Other": "#BDC3C7",
        }

        # Add reaction categories to legend
        legend.node(
            "leg_title2",
            "Reactions (by category):",
            shape="plaintext",
            fontsize="10",
            style="",
        )
        for i, cat in enumerate(categories_present):
            color = category_colors.get(cat, "#BDC3C7")
            legend.node(
                f"leg_rxn_{i}",
                cat,
                shape="box",
                fillcolor=color,
                style="filled,bold",
                fontsize="10",
            )

        # Border style
        legend.node(
            "leg_title3", "Border style:", shape="plaintext", fontsize="10", style=""
        )
        legend.node(
            "leg_style_flex",
            "Bold = Flexible",
            shape="box",
            fillcolor="#3498DB",
            style="filled,bold",
            fontsize="10",
        )
        legend.node(
            "leg_style_blocked",
            "Dashed = Blocked",
            shape="box",
            fillcolor="#3498DB",
            style="filled,dashed",
            fontsize="10",
        )

    # Render multiple formats
    print("Generating network visualizations...")

    output_png = dot.render(output_prefix, format="png", cleanup=True)
    print(f"✓ PNG saved: {output_png}")

    output_svg = dot.render(output_prefix, format="svg", cleanup=True)
    print(f"✓ SVG saved: {output_svg}")

    output_pdf = dot.render(output_prefix, format="pdf", cleanup=True)
    print(f"✓ PDF saved: {output_pdf}")

    print("\n  Network Statistics:")
    flexible_count = sum(1 for _, (_, flex) in reaction_categories.items() if flex)
    blocked_count = sum(1 for _, (_, flex) in reaction_categories.items() if not flex)
    print(f"    Total metabolites: {len(added_mets)}")
    print(f"    Total reactions: {len(ck_reactions)}")
    print(f"    Flexible reactions: {flexible_count}")
    print(f"    Blocked reactions: {blocked_count}")


# ============================================================================
# MAIN
# ============================================================================


def print_usage():
    """Print usage instructions."""
    print(__doc__)


def main():
    """Main function."""
    global MODEL_PATH, PATHWAY_CONFIG, PATHWAY_NAME

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Unified Pathway Network Visualization Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 visualize.py 0
  python3 visualize.py 1 --config config/config_Syndecan1_HS.json
  python3 visualize.py 2 --model ../models/endoA_250917_3_with_Syndecan1_HS.json
  python3 visualize.py --config config/config_Syndecan1_HS.json --model ../models/endoA_250917_3_with_Syndecan1_HS.json
        """,
    )

    parser.add_argument(
        "option",
        type=int,
        nargs="?",
        default=0,
        help="Visualization option (0=all, 1=D3.js, 2=Plotly, 3=GraphViz)",
    )
    parser.add_argument(
        "--config",
        "-c",
        dest="config_path",
        help="Path to pathway configuration JSON file",
    )
    parser.add_argument(
        "--model", "-m", dest="model_path", help="Path to model JSON file"
    )

    args = parser.parse_args()
    option = args.option

    # Load pathway config if provided
    if args.config_path:
        if not os.path.exists(args.config_path):
            print(f"✗ Error: Config file not found: {args.config_path}")
            return
        print(f"Loading pathway config: {args.config_path}")
        PATHWAY_CONFIG = load_pathway_config(args.config_path)

    # Set model path if provided
    if args.model_path:
        if not os.path.exists(args.model_path):
            print(f"✗ Error: Model file not found: {args.model_path}")
            return
        MODEL_PATH = args.model_path

    # Validate option
    if option not in [0, 1, 2, 3]:
        print(f"Error: Invalid option {option}. Must be 0, 1, 2, or 3.")
        print_usage()
        return

    print("\n" + "=" * 80)
    print(f"{PATHWAY_NAME.upper()} NETWORK VISUALIZATION TOOL")
    print("=" * 80)

    # Create output directory
    os.makedirs(VIZ_OUTPUT_DIR, exist_ok=True)

    # Generate visualizations based on option
    if option == 0:
        print("\nGenerating ALL visualizations (1, 2, and 3)...\n")

        # Track which visualizations succeeded
        success_count = 0
        generated_files = []

        # Track files (using dynamic names)
        d3_file = os.path.join(VIZ_OUTPUT_DIR, f"{PATHWAY_NAME}_network_d3.html")
        plotly_file = os.path.join(
            VIZ_OUTPUT_DIR, f"{PATHWAY_NAME}_network_plotly.html"
        )
        png_file = os.path.join(VIZ_OUTPUT_DIR, f"{PATHWAY_NAME}_network.png")
        svg_file = os.path.join(VIZ_OUTPUT_DIR, f"{PATHWAY_NAME}_network.svg")
        pdf_file = os.path.join(VIZ_OUTPUT_DIR, f"{PATHWAY_NAME}_network.pdf")

        # Clear old timestamps to detect new files
        d3_exists_before = os.path.exists(d3_file)
        d3_mtime_before = os.path.getmtime(d3_file) if d3_exists_before else 0

        plotly_exists_before = os.path.exists(plotly_file)
        plotly_mtime_before = (
            os.path.getmtime(plotly_file) if plotly_exists_before else 0
        )

        png_exists_before = os.path.exists(png_file)
        png_mtime_before = os.path.getmtime(png_file) if png_exists_before else 0

        # Try D3.js (no dependencies)
        try:
            visualize_d3js(pathway_config=PATHWAY_CONFIG)
            # Check if file was created/updated
            if os.path.exists(d3_file) and os.path.getmtime(d3_file) > d3_mtime_before:
                success_count += 1
                generated_files.append(f"  ✓ {d3_file} (D3.js interactive)")
        except Exception as e:
            print(f"✗ D3.js failed: {e}")

        # Try Plotly (requires cobra, plotly, networkx)
        try:
            visualize_plotly(pathway_config=PATHWAY_CONFIG)
            # Check if file was created/updated
            if (
                os.path.exists(plotly_file)
                and os.path.getmtime(plotly_file) > plotly_mtime_before
            ):
                success_count += 1
                generated_files.append(f"  ✓ {plotly_file} (Plotly interactive)")
        except Exception:
            pass  # Error already printed by visualize_plotly

        # Try GraphViz (requires cobra, graphviz)
        try:
            visualize_graphviz(pathway_config=PATHWAY_CONFIG)
            # Check if files were created/updated
            graphviz_success = False
            if (
                os.path.exists(png_file)
                and os.path.getmtime(png_file) > png_mtime_before
            ):
                generated_files.append(f"  ✓ {png_file} (GraphViz PNG)")
                graphviz_success = True
            if os.path.exists(svg_file):
                generated_files.append(f"  ✓ {svg_file} (GraphViz SVG)")
            if os.path.exists(pdf_file):
                generated_files.append(f"  ✓ {pdf_file} (GraphViz PDF)")
            if graphviz_success:
                success_count += 1
        except Exception:
            pass  # Error already printed by visualize_graphviz

        print("\n" + "=" * 80)
        print(f"VISUALIZATION SUMMARY: {success_count}/3 completed")
        print("=" * 80)

        if generated_files:
            print("\nSuccessfully generated files:")
            for f in generated_files:
                print(f)
        else:
            print("\nNo files were generated.")

        if success_count < 3:
            print("\nℹ To enable all visualizations, install missing dependencies:")
            print("  pip install cobra plotly networkx graphviz")
            print("  sudo apt-get install graphviz  # System package")

        print("=" * 80 + "\n")

    elif option == 1:
        visualize_d3js(pathway_config=PATHWAY_CONFIG)
        print("\n✓ D3.js visualization complete!")

    elif option == 2:
        visualize_plotly(pathway_config=PATHWAY_CONFIG)
        print("\n✓ Plotly visualization complete!")

    elif option == 3:
        visualize_graphviz(pathway_config=PATHWAY_CONFIG)
        print("\n✓ GraphViz visualization complete!")

    else:
        print(f"Error: Invalid option {option}. Must be 0, 1, 2, or 3.")
        print_usage()
        return


if __name__ == "__main__":
    main()
