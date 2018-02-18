#Take a table of ENCODE psi values.
#First, change the name of every sample to something more readable by relating
#to the sample ids in SamplesTable
#Second, calculate delta psis, using the replicate and kd/control relationships in SamplesTable.
#There are two replicates per kd for each cell line.  Can either keep cell lines separate or use them as a covariate.

import os
import pandas as pd
import sys
import itertools
import numpy as np
from collections import OrderedDict
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, Formula, FloatVector
from rpy2.rinterface import RRuntimeError


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

	#Remove any row with an NA value
	psidf.dropna(axis = 0, how = 'any', inplace = True)
	hepg2df.dropna(axis = 0, how = 'any', inplace = True)
	k562df.dropna(axis = 0, how = 'any', inplace = True)

	psidf.to_csv('~/Desktop/psidf.txt', sep = '\t', index = False)

	return psidf, hepg2df, k562df

def getdpsis(df, samples):
	#Use the control/kd relationships in samples to get delta psis 
	
	#Loop through columns, finding its control column and calculating dpsi
	#Put this dpsi in a new df

	#First column in dpsi dataframe is gene column of psi dataframe
	genes = list(df['Gene'])
	#Remove that annoying dot in the ensembl gene id
	genes = [gene.split('.')[0] for gene in genes]
	d = {'Gene' : genes}
	dpsidf = pd.DataFrame(data = d)
	
	analyzedsamples = []
	for column in df:
		#If this is a control kd column, skip it
		if column == 'Gene' or column.split('_')[0] == 'controlkd':
			continue
		#If this is a replicate of a kd we've already seen, skip it
		samplename = column.split('_')[0] + '_' + column.split('_')[1]
		if samplename in analyzedsamples:
			continue

		controlid = column.split('_')[3]
		samplecolumn = df[column]
		
		#Get kd sample columns
		kddf = psidf.select(lambda col: samplename in col, axis = 1)
		#Get control sample columns
		controldf = psidf.select(lambda col: controlid in col and 'controlkd' in col, axis = 1)
		#Take the mean of the two kd samples
		kddf['kdmean'] = (kddf.iloc[:, 0] + kddf.iloc[:, 1]) / 2
		#Take the mean of the two control samples
		controldf['controlmean'] = (controldf.iloc[:, 0] + controldf.iloc[:, 1]) / 2

		dpsicolumnname = samplename + '_dpsi'
		#dpsi is the difference of the means
		dpsicolumn = pd.Series(kddf['kdmean'] - controldf['controlmean']).reset_index(drop = True) #have to reset index to add column to dpsidf
		dpsidf[dpsicolumnname] = dpsicolumn

		#Add this sample to analyzedsamples
		analyzedsamples.append(samplename)

	
	return dpsidf

def getlmep_singlecellline(df, samples):
	#Use the control/kd relationships to get LME p values.
	#This is for a single cell line df so it will be 2 kd samples vs 2 control samples

	genes = list(df['Gene'])
	#Remove that annoying dot in the ensembl gene id
	genes = [gene.split('.')[0] for gene in genes]
	d = {'Gene' : genes}
	pvaluedf = pd.DataFrame(data = d)

	#Loop through columns, finding its matching replicate and control kd samples
	#Use R nlme package to calculate linear model p values

	pandas2ri.activate() #allow conversion between pandas dataframes and r dataframes

	#Define R packages
	nlme = importr('nlme')
	base = importr('base')
	stats = importr('stats')
	qv = importr('qvalue')

	#define formulae
	fmla = Formula('value ~ 1 + cond1')
	rndm = Formula('~ 1 | samples')
	nullfmla = Formula('value ~ 1')
	nullrndm = Formula('~1 | samples')

	analyzedsamples = []
	for column in df:
		samplename = column.split('_')[0] + '_' + column.split('_')[1]
		controlid = column.split('_')[3]
		#If this is a control kd column, or we've already done this sample (this is the 2nd rep), skip it
		if column == 'Gene' or column.split('_')[0] == 'controlkd' or samplename in analyzedsamples:
			continue

		#Get 2 replicate columns
		kddf = df.select(lambda col: samplename in col, axis = 1)
		#Get 2 control knockdown columns

		#Put kd and control dfs together



samples = readSamplesTable(sys.argv[1])
psidf, hepg2df, k562df = loadpsis(sys.argv[2], samples)
dpsidf = getdpsis(psidf, samples)
print(dpsidf.head())

