#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Helper functions for graph.py to declutter. 
"""

import csv
import os
import numpy as np 
import pandas as pd
from collections import defaultdict

# rdkit
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.warning')
from rdkit.Chem import Draw, rdFMCS

from rdkit.Chem import rdRGroupDecomposition
from rdkit.Chem import rdqueries
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Geometry
rdDepictor.SetPreferCoordGen(True)





# from IPython.display import SVG,Image

def importFEPFile(input_file):
	"""
	Reads a MODSIM-style FEP output file and returns parsed values. 
	The file should follow this format:

	FEP                           dG        sem       crashes   
	FEP_a-b                       float     float     int
	FEP_...                       ...       ...       ...
	FEP_...                       float     float     int

	This function will throw an error if the file diverges from format or
	if the matrix contains empty fields.    
	"""
	# checks.
	if not input_file:
		raise Exception("Path to input file is not defined.")
	elif not os.path.exists(input_file):
		raise Exception("File does not exist: {}".format(input_file))

	perts_dict = {}		# dict with pertname : freenrg, sem, crash

	# parse input file.
	with open(input_file, "r") as inputfile:
		reader = csv.reader(inputfile)

		# skip header.
		next(reader)
		for row in reader:

			# input file is written as space sep. values.
			try:
				pert, freenrg, sem, crash = row[0].rsplit()
			except ValueError:
				raise Exception("Input data contains empty field(s) on row {}".format(row))

			# set key as perturbation while removing the FEP regex.
			key = pert.replace("FEP_","")

			perts_dict[key] = freenrg, sem, crash

	# with check in row loop values shouldn't be None at this point.
	return perts_dict

def computeDDGs(bound_path, solvated_path):
	"""
	Given two FEP output files, generates a dictionary with ddG, SEM and number of crashes per 
	perturbation. 
	"""

	# daughter function contains file tests.
	p_perts = importFEPFile(bound_path)
	w_perts = importFEPFile(solvated_path)

	if not len(p_perts) == len(w_perts):
		raise Exception("Bound and solvated FEP output files an contain uneven number"
		" of perturbations: {}!={}".format(len(p_perts), len(w_perts)))

	ddG_dict = {}	# dict into which we collect ddG information with perturbation names as keys.

	# iterate over perturbations.
	for p_pert, (p_freenrg, p_sem, p_crashes) in p_perts.items():

		# grab the solvated information to go with this bound information by perturbation name.
		try:
			(w_freenrg, w_sem, w_crashes) = w_perts[p_pert]
		except KeyError:
			raise Exception("Bound perturbation \"{}\" not found in solvated FEP output"
			" file.".format(p_pert))

		# Define our final values. 
		pert_name = p_pert
		freenrg_ddG = float(w_freenrg) - float(p_freenrg)

		# propagate SEM.
		sem_ddG = (float(p_sem)+float(w_sem))/np.sqrt(2)

		# sum number of crashes
		crashes_ddG = int(p_crashes) + int(w_crashes)

		# add to dict.
		ddG_dict[p_pert] = freenrg_ddG, sem_ddG, crashes_ddG
	return ddG_dict

def generateImages(mol,row,core,width=350,height=200,
                      fillRings=False,legend="",
                      sourceIdxProperty="SourceAtomIdx",
                      lbls=None):
    # copy the molecule and core
    mol = Chem.Mol(mol)
    core = Chem.Mol(core)

    # -------------------------------------------
    # include the atom map numbers in the substructure search in order to 
    # try to ensure a good alignment of the molecule to symmetric cores
    for at in core.GetAtoms():
        if at.GetAtomMapNum():
            at.ExpandQuery(rdqueries.IsotopeEqualsQueryAtom(200+at.GetAtomMapNum()))
            
    for lbl in row:
        if lbl=='Core':
            continue
        rg = row[lbl]
        for at in rg.GetAtoms():
            if not at.GetAtomicNum() and at.GetAtomMapNum() and \
            at.HasProp('dummyLabel') and at.GetProp('dummyLabel')==lbl:
                # attachment point. the atoms connected to this
                # should be from the molecule
                for nbr in at.GetNeighbors():
                    if nbr.HasProp(sourceIdxProperty):
                        mAt = mol.GetAtomWithIdx(nbr.GetIntProp(sourceIdxProperty))
                        if mAt.GetIsotope():
                            mAt.SetIntProp('_OrigIsotope',mAt.GetIsotope())
                        mAt.SetIsotope(200+at.GetAtomMapNum())
    # remove unmapped hs so that they don't mess up the depiction
    rhps = Chem.RemoveHsParameters()
    rhps.removeMapped = False
    tmol = Chem.RemoveHs(mol,rhps)
    rdDepictor.GenerateDepictionMatching2DStructure(tmol,core)

    oldNewAtomMap={}
    # reset the original isotope values and account for the fact that
    # removing the Hs changed atom indices
    for i,at in enumerate(tmol.GetAtoms()):
        if at.HasProp(sourceIdxProperty):
            oldNewAtomMap[at.GetIntProp(sourceIdxProperty)] = i
            if at.HasProp("_OrigIsotope"):
                at.SetIsotope(at.GetIntProp("_OrigIsotope"))
                at.ClearProp("_OrigIsotope")
            else:
                at.SetIsotope(0)
      
    # ------------------
    #  set up our colormap
    #   the three choices here are all "colorblind" colormaps
    
    # "Tol" colormap from https://davidmathlogic.com/colorblind
    colors = [(51,34,136),(17,119,51),(68,170,153),(136,204,238),(221,204,119),(204,102,119),(170,68,153),(136,34,85)]
    # "IBM" colormap from https://davidmathlogic.com/colorblind
    colors = [(100,143,255),(120,94,240),(220,38,127),(254,97,0),(255,176,0)]
    # Okabe_Ito colormap from https://jfly.uni-koeln.de/color/
    colors = [(230,159,0),(86,180,233),(0,158,115),(240,228,66),(0,114,178),(213,94,0),(204,121,167)]
    for i,x in enumerate(colors):
        colors[i] = tuple(y/255 for y in x)
  
    #----------------------
    # Identify and store which atoms, bonds, and rings we'll be highlighting
    highlightatoms = defaultdict(list)
    highlightbonds = defaultdict(list)
    atomrads = {}
    widthmults = {}

    rings = []
    # loop over R groups.
    for i,lbl in enumerate(lbls):    
        color = colors[i%len(colors)]
        try:
            rquery = row[lbl]
        # we don't know the number of R-groups, so just quit this loop if we've reached the end.
        except KeyError:
        	continue

        Chem.GetSSSR(rquery)
        rinfo = rquery.GetRingInfo()
        for at in rquery.GetAtoms():
            if at.HasProp(sourceIdxProperty):
                origIdx = oldNewAtomMap[at.GetIntProp(sourceIdxProperty)]
                highlightatoms[origIdx].append(color)
                atomrads[origIdx] = 0.4
        if fillRings:
            for aring in rinfo.AtomRings():
                tring = []
                allFound = True
                for aid in aring:
                    at = rquery.GetAtomWithIdx(aid)
                    if not at.HasProp(sourceIdxProperty):
                        allFound = False
                        break
                    tring.append(oldNewAtomMap[at.GetIntProp(sourceIdxProperty)])
                if allFound:
                    rings.append((tring,color))
        for qbnd in rquery.GetBonds():
            batom = qbnd.GetBeginAtom()
            eatom = qbnd.GetEndAtom()
            if batom.HasProp(sourceIdxProperty) and eatom.HasProp(sourceIdxProperty):
                origBnd = tmol.GetBondBetweenAtoms(oldNewAtomMap[batom.GetIntProp(sourceIdxProperty)],
                                                 oldNewAtomMap[eatom.GetIntProp(sourceIdxProperty)])
                bndIdx = origBnd.GetIdx()
                highlightbonds[bndIdx].append(color)
                widthmults[bndIdx] = 2

    d2d = rdMolDraw2D.MolDraw2DCairo(width,height)
    dos = d2d.drawOptions()
    dos.useBWAtomPalette()
                
    #----------------------
    # if we are filling rings, go ahead and do that first so that we draw
    # the molecule on top of the filled rings
    if fillRings and rings:
        # a hack to set the molecule scale
        d2d.DrawMoleculeWithHighlights(tmol,legend,dict(highlightatoms),
                                       dict(highlightbonds),
                                       atomrads,widthmults)
        d2d.ClearDrawing()
        conf = tmol.GetConformer()
        for (aring,color) in rings:
            ps = []
            for aidx in aring:
                pos = Geometry.Point2D(conf.GetAtomPosition(aidx))
                ps.append(pos)
            d2d.SetFillPolys(True)
            d2d.SetColour(color)
            d2d.DrawPolygon(ps)
        dos.clearBackground = False

    #----------------------
    # now draw the molecule, with highlights:
    d2d.DrawMoleculeWithHighlights(tmol,"",dict(highlightatoms),dict(highlightbonds),
                                   atomrads,widthmults)
    d2d.FinishDrawing()
    png = d2d.GetDrawingText()

    # save to file:
    d2d.WriteDrawingText("data/mol_images/"+legend+".png")
        
    return png

def writeLigandImages(ligand_names):
	"""
	Takes a list of ligand names and generates .png images in ./data/mol_images/ using RDKit.
	Workflow mostly based on https://rdkit.blogspot.com/2020/10/molecule-highlighting-and-r-group.html
	"""

	##! Change ligand SD file path during integration!
	ligands_path = "input/CDK2_ligands.sdf"

	if len(ligand_names) < 2:
		raise Exception("Number of ligands ({}) should be > 1".format(len(ligand_names)))

	# unnest the ligand names:
	ligand_names = [lig[0] for lig in ligand_names]

	# we work with multimol SD files. Use a supplier to generate a list of ligand names and 
	# a list of molecules (while preserving order).

	ligand_names = []		# list of ligand names (str)
	rdkit_mols = []			# list of rdkit objects: <rdkit.Chem.rdchem.Mol object at ...>

	suppl = Chem.SDMolSupplier(ligands_path)
	for mol in suppl:
		ligand_names.append(mol.GetProp("_Name"))
		rdkit_mols.append(mol)

	################## RDKIT IMAGE GENERATION ################
	# flatten all ligands into 2D space. Note that most rdkit mol operations are in-place.
	for m in rdkit_mols:
		rdDepictor.Compute2DCoords(m)
		m.UpdatePropertyCache()


	# find the core using standard MCS algorithm. Flatten the core again just to be sure.
	mcs = rdFMCS.FindMCS(rdkit_mols, matchValences=False,
                                    ringMatchesRingOnly=True,
                                    completeRingsOnly=True,
                                    matchChiralTag=False)

	core = Chem.MolFromSmarts(mcs.smartsString)
	rdDepictor.Compute2DCoords(core)

	# find subtructure matches with MCS per ligand, then tag matching atom indices in each ligand.
	ps = Chem.AdjustQueryParameters.NoAdjustments()
	ps.makeDummiesQueries=True
	qcore = Chem.AdjustQueryProperties(core,ps)
	#mhs = [Chem.AddHs(x,addCoords=True) for x in ms]
	mms = [x for x in rdkit_mols if x.HasSubstructMatch(qcore)]
	for m in mms:
		for atom in m.GetAtoms():
			atom.SetIntProp("SourceAtomIdx",atom.GetIdx())

	# do an RDKit R-group decomposition.	
	groups,_ = rdRGroupDecomposition.RGroupDecompose([qcore],mms,asSmiles=False,asRows=True)

	# call the writer function with each molecule.
	for i, m in enumerate(rdkit_mols):

		png = generateImages(m,groups[i],qcore,lbls=('R1','R2','R3','R4', 'R5', 'R6', 'R7'),
								legend=ligand_names[i],
								width=400,height=400)








































