#Take a table of ENCODE psi values.
#First, change the name of every sample to something more readable by relating
#to the sample ids in SamplesTable
#Second, calculate delta psis, using the replicate and kd/control relationships in SamplesTable.
#There are two replicates per kd for each cell line.  Can either keep cell lines separate or use them as a covariate.

import os
import pandas as pd
import sys
import itertools


def readSamplesTable(samplestable):
	samples = {} #{sampleid: [KdGene, cellline, controlID]}
	s = []
	with open(samplestable, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			#Skip header
			if line[0] == 'ID':
				continue
			if len(line) == 8:
				sampleid = line[1]
				kdgene = line[3]
				cellline = line[4]
				controlID = line[7].split('/')[2]
			elif len(line) == 7:
				#If this is a control sample line
				sampleid = line[1]
				kdgene = 'controlkd'
				cellline = line[4]
				controlID = sampleid

			samples[sampleid] = [kdgene, cellline, controlID]

	return samples

def loadpsis(psis, samples):
	#Load psi table and rename columns
	#Split into a HepG2 table and a K562 table
	renamedict = {} #{oldname : newname}
	psidf = pd.read_table(psis, sep = '\t', header = 0, index_col = False)
	for column in psidf:
		if column == 'Gene':
			continue
		sampleID = column.split('_')[0]
		rep = column.split('_')[1]
		kdgene = samples[sampleID][0]
		cellline = samples[sampleID][1]
		controlID = samples[sampleID][2]
		renamedict[column] = kdgene + '_' + cellline + '_' + rep + '_' + controlID
	psidf.rename(columns = renamedict, inplace = True)

	#Split into HepG2 and K562 dfs
	hepg2df = psidf.select(lambda col: 'Gene' in col or 'HepG2' in col, axis = 1)
	k562df = psidf.select(lambda col: 'Gene' in col or 'K562' in col, axis = 1)

	return hepg2df, k562df

samples = readSamplesTable(sys.argv[1])
loadpsis(sys.argv[2], samples)

