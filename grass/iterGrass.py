import sys
import os
import argparse
import itertools
import statistics
import re
import csv
import fileinput
import numpy as np
from collections import defaultdict
from subprocess import call
import pandas as pd
from collections import Counter

alpha = 0.8 #ratio alloted to the probability value for labels
newEdgeProb = 0.9 #add a new edge if the shared label prob is greater than this


def flattenClusters(infile, outfile):
    with open(outfile, 'w') as ofile:
        with open(infile) as ifile:
            for i,l in enumerate(ifile):
                toks = l.rstrip().split()
                cname = "cluster{}".format(i)
                for t in toks:
                    ofile.write("{}\t{}\n".format(cname, t))
                    
def getMedWeight(graph, node1, node2):
	weights = []
	for (x, weight) in graph[node1]:
		if weight != 1.1:
			weights.append(weight)
		else:
			weights.append(1)
	for (x, weight) in graph[node2]:
		if weight != 1.1:
			weights.append(weight)
		else:
			weights.append(1)
	
	if not weights:
		return(0)
	else:
		return(statistics.median(weights))

def readLabels(keys):
	#below we get labelToContigs = labels -> contigs that map to this label
	#and contigToLabels = contigs -> all labels and their corresponding probs
	labelToContigs = defaultdict(set)
	contigToLabels = defaultdict(set)
	with open(keys["seed_file"], 'r') as ifile:
		for line in ifile:
			data = (line.strip('\n')).split('\t')
			contigName = data[0]
			data = data[3].split()
			curLabel = ""
			for i in range(1, len(data), 2):
				curLabel = data[i-1]
				if (curLabel != "" and curLabel != "__DUMMY__" and data[i].lower() != "nan"):
					contigToLabels[contigName].add((curLabel, float(data[i])))
					labelToContigs[curLabel].add((contigName, float(data[i])))
	return(labelToContigs, contigToLabels)

#get the sum of the probabilities that 2 nodes have the same label
def getProb(contigToLabels, node1, node2):
	labels1 = []
	probs1 = []
	totalProbability = 0
	for (label, prob) in contigToLabels[node1]:
		labels1.append(label)
		probs1.append(prob)
	for (label, prob) in contigToLabels[node2]:
		if label in labels1:
			totalProbability += (prob * probs1[labels1.index(label)])
	return totalProbability

def changeEdgeWeights(orgGraph, graph, contigToLabels, ofile):
	changesMade = 0
	noChange = 0
	weightCalc = 0;
	for (node1, node2), weight in graph.iteritems():
		if node1 == node2:
			ofile.write(node1 + "\t" + node2 + "\t" + str(weight) + "\n")
			noChange += 1
		else:
			prob = 0
			if (node1 in contigToLabels) and (node2 in contigToLabels):
				for (label1, prob1) in  contigToLabels[node1]:
					for (label2, prob2) in contigToLabels[node2]:
						if label1 == label2:
							prob += (prob1 * prob2)
			if prob > 0:
				for (x, w) in orgGraph[node1]:
					if x == node2:
						orgWeight = w
				for (x, w) in orgGraph[node2]:
					if x == node1:
						orgWeight = w
				newWeight = ((1 - alpha) * orgWeight) + (prob * alpha)
				changesMade += 1
			else:
				newWeight = weight
				noChange += 1
			weightCalc += newWeight
			ofile.write(node1 + "\t" + node2 + "\t" + str(newWeight) + "\n")
	weightCalc /= (changesMade+noChange)
	return (weightCalc, (changesMade + noChange))

def addNewEdges(orgGraph, graph, contigToLabels, labelToContigs, ofile):
	changesMade = 0
	weightCalc = 0
	for label in labelToContigs.iterkeys():
		contigs = []
		probs = []
		for (contig, prob) in labelToContigs[label]:
			contigs.append(contig)
			probs.append(prob)
		edgeExists = 0
		for node1, node2 in itertools.combinations(contigs, 2):
			if (node1, node2) in graph or (node2, node1) in graph:
				edgeExists = 1
			if edgeExists == 0:
				totalProb = getProb(contigToLabels, node1, node2)
				orgWeight = getMedWeight(orgGraph, node1, node2)
				newWeight = ((1 - alpha) * orgWeight) + (totalProb * alpha)
				if (totalProb) > newEdgeProb:
					ofile.write(node1 + "\t" + node2 + "\t" + str(newWeight) + "\n")
					orgGraph[node1].add((node2, newWeight))
					changesMade += 1
					weightCalc += newWeight
	weightCalc /= changesMade
	return (weightCalc, changesMade)

def run(keys, labelFile, juntoConfigFile, outdir):
	netFile = os.path.sep.join([outdir, "mag.filt.net"])
	
	if (netFile == keys["graph_file"]):
		print ("ERROR: junto graph file should be different from the net file otherwise it will be overwritten")
		return 0

	orgGraph = defaultdict(set)
	orgGraphSize = 0
	with open(netFile, 'r') as ifile:
		for line in ifile:
			edge = (line.strip('\n')).split('\t')
			orgGraph[edge[0]].add((edge[1], float(edge[2])))
			orgGraphSize += 1

	diff = 10000
	orgDiff = 0
	iters = 1
	i = 0
	call(["cp", netFile, keys["graph_file"]])

	while (diff > (0.05*orgDiff)):
		i += 1
		print ("Started iteration number: " + str(i) + "\n")

		if (i == 2):
			orgDiff = diff;

		call(["junto", "config", juntoConfigFile])
		call(["cp", keys["output_file"], keys["seed_file"]])

		graph = {}
		graphSize = 0
		with open(keys["graph_file"], 'r') as ifile:
			for line in ifile:
				edge = (line.strip('\n')).split('\t')
				graph[(edge[0], edge[1])] = float(edge[2])
				graphSize += 1

		labelToContigs, contigToLabels = readLabels(keys)

		sizeNewGraph = 0

		with open(keys["graph_file"], 'w') as ofile:
			#write previous edges with new weights
			(avgOldWeight, temp) = changeEdgeWeights(orgGraph, graph, contigToLabels, ofile)
			sizeNewGraph += temp
			#add new edges between contigs with same labels
			(avgNewWeight, temp) = addNewEdges(orgGraph, graph, contigToLabels, labelToContigs, ofile)
			sizeNewGraph += temp
		
		with open(keys["seed_file"], 'w') as ofile:
			for contig in contigToLabels.iterkeys():
				for (label, prob) in contigToLabels[contig]:
					ofile.write(contig + "\t" + label + "\t" + str(prob) + "\n")
		
		diff = abs(sizeNewGraph - graphSize)
		print("Number of edges added/removed: " + str(diff) + "\n")
		
		if ((avgOldWeight/avgNewWeight) > 0.5): 
			iters+=1
		
		with open(juntoConfigFile, "w") as configFile:
                	for x in keys.iterkeys():
                        	configFile.write(x + " = " + keys[x] + "\n")
                	configFile.write("data_format = edge_factored\n")
                	configFile.write("iters = " + str(iters) + "\n")
               		configFile.write("prune_threshold = 0\n")
                	configFile.write("algo = adsorption\n")