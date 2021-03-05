#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
From a standard MODSIM-style FEP output (both bound and solvated) text output file, generates molecule 
images and computes a networkx graph. This graph is written to JSON format in a style that is 
parsable by networkgen.html
"""

# internal imports
import helper_functions

# general imports

# process specific imports
import networkx as nx
from networkx.readwrite import json_graph
import json


if __name__ == "__main__":

	##! change paths to universal naming at time of integration.
	protein_fep_file = "input/protein_cdk2_opls.txt"
	water_fep_file = "input/water_cdk2_opls.txt"

	# retrieve a dictionary of perturbations with corresponding
	# relative free energies of binding (ddG), standard errors of mean
	# (SEMs) and number of crashed simulations (crashes).
	ddG_dict = helper_functions.computeDDGs(
										bound_path=protein_fep_file, 
										solvated_path=water_fep_file
										)

	# get unique lists of ligands and edges
	nodes = [] 
	edges = []

	for pert in ddG_dict.keys():
		node0 = pert.split("-")[0]
		node1 = pert.split("-")[1]

		nodes.append(node0)
		nodes.append(node1)
		edges.append((node0,node1))

	# get unique nodes:
	nodes = set(nodes)

	################ NETWORKX GRAPH GENERATION ######################
	G = nx.Graph()
	G.add_nodes_from(nodes)
	G.add_edges_from(edges)

	# add information to nodes.
	for node in G.nodes:

		# keyword for visjs to recognise that node should be a circular image.
		G.nodes[node]["shape"] = "circularImage"

		# add image path.
		##! change path to web path on integration! github link is for dev only.
		G.nodes[node]["image"] = "https://raw.githubusercontent.com/GPCR-ModSim/qfepweb"+\
				"/networkgen/networkgen_playground/visjs/imgs/{}.png".format(node)

	# add information to edges.

	for edge in G.edges:

		# get edge data from initial dict.
		# networkx sometimes changes the order of nodes in edge naming, so try both ways. 
		# this is ok because edges are unique in FEP anyway.
		try:
			pert_name = "{}-{}".format(edge[0],edge[1])
			freenrg, sem, crashes = ddG_dict[pert_name]
		except KeyError:
			pert_name = "{}-{}".format(edge[1],edge[0])
			freenrg, sem, crashes = ddG_dict[pert_name]			

		# similar to method with nodes.
		G.edges[edge]["label"] = str(round(freenrg, 2)) + "\n±" + str(round(sem, 1))
		G.edges[edge]["freenrg"] = str(round(freenrg, 2))
		G.edges[edge]["sem"] = str(round(sem, 2))
		G.edges[edge]["crashes"] = str(crashes)

	# networkx to json dict; use keywords that conform to Visjs terminology.
	graph_data = json_graph.node_link_data(G, {"link": "edges", "source": "from", "target": "to"})

	# we only need node + edge information, remove redundant networkx keys.
	graph_data = {k: graph_data[k] for k in ("nodes", "edges")}

	with open("data/tmp.json", "w") as jsonfile:
		json.dump(graph_data, jsonfile, indent=4)

	############# RDKIT MOLECULE IMAGE GENERATION ###################

	# compile nested list of ligand names so that we can check if there are >1 ligands.
	ligands = [[node] for node in G.nodes]
	helper_functions.writeLigandImages(ligands)

	# load the SD supplier in input, 

	# use notebook functions to write image files.

	# write CIRCULAR PNGs

	# adjust visjs to load node shape as image shape.












