#Take a table of ENCODE psi values.
#First, change the name of every sample to something more readable by relating
#to the sample ids in SamplesTable
#Second, calculate delta psis, using the replicate and kd/control relationships in SamplesTable.
#There are two replicates per kd for each cell line.  Can either keep cell lines separate or use them as a covariate.
#Third, calculate nlme pvalues.

import os
import pandas as pd
import sys
import itertools
import numpy as np
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, Formula, FloatVector
from rpy2.rinterface import RRuntimeError
import warnings
import argparse


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

	return psidf, hepg2df, k562df

def getdpsis_onecellline(df):
	#Use the control/kd relationships in samples to get delta psis 
	
	#Loop through columns, finding its control column and calculating dpsi
	#Put this dpsi in a new df

	#First column in dpsi dataframe is gene column of psi dataframe
	genes = list(df['Gene'])
	print(len(genes))
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
		kddf = df.select(lambda col: samplename in col, axis = 1)
		print(kddf.shape)
		#Get control sample columns
		controldf = df.select(lambda col: controlid in col and 'controlkd' in col, axis = 1)
		#print(controldf.shape)
		#Take the mean of the two kd samples
		kddf['kdmean'] = kddf.mean(axis = 1)
		#Take the mean of the two control samples
		controldf['controlmean'] = controldf.mean(axis = 1)

		dpsicolumnname = samplename + '_dpsi'
		#dpsi is the difference of the means
		dpsicolumn = pd.Series(kddf['kdmean'] - controldf['controlmean']).reset_index(drop = True) #have to reset index to add column to dpsidf
		dpsidf[dpsicolumnname] = dpsicolumn

		#Add this sample to analyzedsamples
		analyzedsamples.append(samplename)

	
	return dpsidf

def getdpsis_twocelllines(df):
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
		samplename = column.split('_')[0]
		if samplename in analyzedsamples:
			continue
		
		#Get kd sample columns
		kddf = df.select(lambda col: samplename in col, axis = 1).reset_index(drop = True)
		
		#Sometimes this kd wasn't done in both HepG2 and K562
		#If it wasn't, we have to skip it
		kdsamples = kddf.columns.values.tolist()
		hepg2kds = [kd for kd in kdsamples if 'HepG2' in kd]
		k562kds = [kd for kd in kdsamples if 'K562' in kd]
		if len(hepg2kds) < 2 or len(k562kds) < 2:
			continue

		#Get control sample IDs for the 4 knockdown samples
		controlIDs = []
		for column in kddf.columns.values.tolist():
			controlID = column.split('_')[3]
			controlIDs.append(controlID)
		controlIDs = list(set(controlIDs))

		#Get control sample columns
		controldf = df.select(lambda col: 'controlkd' in col and col.split('_')[3] in controlIDs, axis = 1).reset_index(drop = True)
		
		#Take the mean of the kd samples
		for cellline in ['HepG2', 'K562']:
			kd_celllinedf = kddf.select(lambda col: cellline in col, axis = 1)
			kddf['{0}_kdmean'.format(cellline)] = kd_celllinedf.mean(axis = 1)
			#Take the mean of the control samples
			control_celllinedf = controldf.select(lambda col: cellline in col, axis = 1)
			controldf['{0}_controlmean'.format(cellline)] = control_celllinedf.mean(axis = 1)

		#Take the mean of the two cell line means
		kddf['kdmean'] = (kddf.loc[:, 'HepG2_kdmean'] + kddf.loc[:, 'K562_kdmean']) / 2
		controldf['controlmean'] = (controldf.loc[:, 'HepG2_controlmean'] + controldf.loc[:, 'K562_controlmean']) / 2

		dpsicolumnname = samplename + '_dpsi'
		#dpsi is the difference of the means
		dpsicolumn = pd.Series(kddf['kdmean'] - controldf['controlmean']).reset_index(drop = True) #have to reset index to add column to dpsidf
		dpsidf[dpsicolumnname] = dpsicolumn

		#Add this sample to analyzedsamples
		analyzedsamples.append(samplename)
	
	return dpsidf

def getlmep_singlecellline(df):
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
	samplecounter = 0
	for column in df:
		#Going to make a dataframe of kd and control samples.
		sampledf = {'Gene' : genes}
		sampledf = pd.DataFrame(data = sampledf)

		
		#If this is a control kd column, skip it
		if column == 'Gene' or column.split('_')[0] == 'controlkd':
			continue

		#If we've already seen this sample (it's the second rep), skip it
		samplename = column.split('_')[0] + '_' + column.split('_')[1]
		controlid = column.split('_')[3]
		if samplename in analyzedsamples:
			continue

		samplecounter +=1
		print('Calculating pvalues for {0} (number {1})...'.format(samplename, samplecounter))
		
		#Get 2 replicate columns
		kddf = df.select(lambda col: samplename in col, axis = 1).reset_index(drop = True)
		kddfcolumns = kddf.columns.values.tolist()
		kddfcolumns = [colname + '_kd' for colname in kddfcolumns]
		kddf.columns = kddfcolumns
		#Get 2 control knockdown columns
		controldf = df.select(lambda col: controlid in col and 'controlkd' in col, axis = 1).reset_index(drop = True)
		#Get the number of total samples (kd + control)
		nsamps = len(kddf.columns.tolist()) + len(controldf.columns.tolist())
		#Put kd and control dfs together
		sampledf = pd.concat([sampledf, kddf, controldf], axis = 1)
		
		#List of pvalues so that we can later turn them into qvalues
		pvalues = []
		#Loop through rows (genes) of sampledf and calculate pvalues
		for index, row in sampledf.iterrows():
			d = {}
			d['Gene'] = [row['Gene']] * nsamps
			d['variable'] = sampledf.columns.values.tolist()[1:] #first column is 'Gene'

			psivalues = [] #get psi values for this gene
			for sample in sampledf.columns.values.tolist()[1:]:
				psivalue = row[sample]
				psivalues.append(psivalue)
			d['value'] = psivalues

			#Get the condition (kd or control) for the samples for this gene
			samps = sampledf.columns.values.tolist()[1:]
			cond1 = [1 if 'control' in samp else 0 for samp in samps] #0 if kd, 1 if control
			cond2 = [0 if 'control' in samp else 1 for samp in samps] #1 if kd, 0 if control
			d['cond1'] = cond1
			d['cond2'] = cond2

			d['samples'] = [x + 1 for x in range(nsamps)]

			#Turn this dictionary into a DataFrame
			rowdf = pd.DataFrame.from_dict(d)
			
			#Get LME p value
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				try:
					lm_alt = nlme.lme(fmla, random = rndm, data = rowdf, method = 'ML') #test
					lm_null = nlme.lme(nullfmla, random = nullrndm, data = rowdf, method = 'ML') #control
					logratio = (stats.logLik(lm_alt)[0] - stats.logLik(lm_null)[0]) * 2
					pvalue = stats.pchisq(logratio, df = 1, lower_tail = False)[0]
					#format decimal
					pvalue = float('{:.2e}'.format(pvalue))
				except RRuntimeError:
					print('RRuntime error for {0}/{1}!'.format(row['Gene'], samplename))
					pvalue = 1.0

			pvalues.append(pvalue)

		#Turn pvalues into qvalues
		#Turn list of pvalues into R vector
		pvec = FloatVector(pvalues)
		#Get qvalues object
		qobj = qv.qvalue(p = pvec)
		#qvalues are index 2 of qvalue object
		qvalues = list(qobj[2])
		#format decimal
		qvalues = [float('{:.2e}'.format(qvalue)) for qvalue in qvalues]
		#Add pvalues and qvalues to pvaluedf
		pvalues = pd.Series(pvalues).reset_index(drop = True)
		qvalues = pd.Series(qvalues).reset_index(drop = True)
		pvaluedf['{0}_pval'.format(samplename)] = pvalues
		pvaluedf['{0}_qval'.format(samplename)] = qvalues
		#Add this sample to analyzedsamples so that we don't do the same thing all over again for replicate 2
		analyzedsamples.append(samplename)

	return pvaluedf

def getlmepdf_twocelllines(df):
	#This is just like getlmep_singlecellline, except that we are calculating one p value per sample across both (HepG2 and K562)
	#cell lines. The cell line will be used as a covariate. It only makes sense to give this function the combined psidf,
	#which has both cell lines.

	#Use the control/kd relationships to get LME p values.
	#This is for a combined cell line df so it will be 4 kd samples vs 4 control samples

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
	fmla = Formula('value ~ 1 + cond1 + cellline1')
	rndm = Formula('~ 1 | samples')
	nullfmla = Formula('value ~ 1 + cellline1')
	nullrndm = Formula('~1 | samples')

	samplecounter = 0
	analyzedsamples = []
	for column in df:
		#Going to make a dataframe of kd and control samples.
		sampledf = {'Gene' : genes}
		sampledf = pd.DataFrame(data = sampledf)

		
		#If this is a control kd column, skip it
		if column == 'Gene' or column.split('_')[0] == 'controlkd':
			continue

		#If we've already seen this sample (it's the second rep), skip it
		samplename = column.split('_')[0]
		cellline = column.split('_')[1]
		if samplename in analyzedsamples:
			continue

		samplecounter +=1
		print('Calculating pvalues for {0} (number {1})...'.format(samplename, samplecounter))
		
		#Get 2 replicate columns
		kddf = df.select(lambda col: samplename in col, axis = 1).reset_index(drop = True)
		kddfcolumns = kddf.columns.values.tolist()
		kddfcolumns = [colname + '_kd' for colname in kddfcolumns]
		kddf.columns = kddfcolumns
		
		#Get control sample IDs for the 4 knockdown samples
		controlIDs = []
		for column in kddf.columns.values.tolist():
			controlID = column.split('_')[3]
			controlIDs.append(controlID)
		controlIDs = list(set(controlIDs))

		#Get 2 control knockdown columns
		controldf = df.select(lambda col: 'controlkd' in col and col.split('_')[3] in controlIDs, axis = 1).reset_index(drop = True)
		#Get the number of total samples (kd + control)
		nsamps = len(kddf.columns.tolist()) + len(controldf.columns.tolist())

		#Put kd and control dfs together
		sampledf = pd.concat([sampledf, kddf, controldf], axis = 1)
		
		#List of pvalues so that we can later turn them into qvalues
		pvalues = []
		#Loop through rows (genes) of sampledf and calculate pvalues
		for index, row in sampledf.iterrows():
			d = {}
			d['Gene'] = [row['Gene']] * nsamps
			d['variable'] = sampledf.columns.values.tolist()[1:] #first column is 'Gene'

			psivalues = [] #get psi values for this gene
			for sample in sampledf.columns.values.tolist()[1:]:
				psivalue = row[sample]
				psivalues.append(psivalue)
			d['value'] = psivalues

			#Get the condition (kd or control) for the samples for this gene
			samps = sampledf.columns.values.tolist()[1:]
			cond1 = [1 if 'control' in samp else 0 for samp in samps] #0 if kd, 1 if control
			cond2 = [0 if 'control' in samp else 1 for samp in samps] #1 if kd, 0 if control
			d['cond1'] = cond1
			d['cond2'] = cond2

			#Get the cell line (HepG2 or K562) for the samples for this gene
			samps = sampledf.columns.values.tolist()[1:]
			cellline1 = [1 if 'HepG2' in samp else 0 for samp in samps] #1 if HepG2, 0 if K562
			cellline2 = [0 if 'HepG2' in samp else 1 for samp in samps] #0 if HepG2, 1 if K562

			d['cellline1'] = cellline1
			d['cellline2'] = cellline2

			d['samples'] = [x + 1 for x in range(nsamps)]

			#Turn this dictionary into a DataFrame
			rowdf = pd.DataFrame.from_dict(d)

			#Get LME p value
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				try:
					lm_alt = nlme.lme(fmla, random = rndm, data = rowdf, method = 'ML') #test
					lm_null = nlme.lme(nullfmla, random = nullrndm, data = rowdf, method = 'ML') #control
					logratio = (stats.logLik(lm_alt)[0] - stats.logLik(lm_null)[0]) * 2
					pvalue = stats.pchisq(logratio, df = 1, lower_tail = False)[0]
					#format decimal
					pvalue = float('{:.2e}'.format(pvalue))
				except RRuntimeError:
					print('RRuntime error for {0}/{1}!'.format(row['Gene'], samplename))
					pvalue = 1.0

			pvalues.append(pvalue)

		#Turn pvalues into qvalues
		#Turn list of pvalues into R vector
		pvec = FloatVector(pvalues)
		#Get qvalues object
		qobj = qv.qvalue(p = pvec)
		#qvalues are index 2 of qvalue object
		qvalues = list(qobj[2])
		#format decimal
		qvalues = [float('{:.2e}'.format(qvalue)) for qvalue in qvalues]
		#Add pvalues and qvalues to pvaluedf
		pvalues = pd.Series(pvalues).reset_index(drop = True)
		qvalues = pd.Series(qvalues).reset_index(drop = True)
		pvaluedf['{0}_pval'.format(samplename)] = pvalues
		pvaluedf['{0}_qval'.format(samplename)] = qvalues
		#Add this sample to analyzedsamples so that we don't do the same thing all over again for the other samples for this kd
		analyzedsamples.append(samplename)

	return pvaluedf
			

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--mode', type = str, choices = ['dpsi_oneline', 'dpsi_twolines', 'pval_oneline', 'pval_twolines'])
	parser.add_argument('--samplestable', type = str, help = 'File relating sample IDs to their controls. Usually SamplesTable.txt.')
	parser.add_argument('--psitable', type = str, help = 'Table of PSI values for all samples and all genes. Output of TXpsi.py.')
	args = parser.parse_args()

	samples = readSamplesTable(args.samplestable)
	psidf, hepg2df, k562df = loadpsis(args.psitable, samples)
	if args.mode == 'dpsi_oneline':
		hepg2dpsidf = getdpsis_onecellline(hepg2df)
		hepg2dpsidf.to_csv('HepG2.dpsi.txt', sep = '\t', header = True, index = False)
		k562dpsidf = getdpsis_onecellline(k562df)
		k562dpsidf.to_csv('K562.dpsi.txt', sep = '\t', header = True, index = False)

	elif args.mode == 'dpsi_twolines':
		dpsidf = getdpsis_twocelllines(psidf)
		dpsidf.to_csv('HepG2andK562.dpsi.txt', sep = '\t', header = True, index = False)

	elif args.mode == 'pval_oneline':
		hepg2lmepdf = getlmep_singlecellline(hepg2df)
		hepg2lmepdf.to_csv('HepG2.psi.pval.txt', sep = '\t', header = True, index = False)
		k562lmepdf = getlmep_singlecellline(k562df)
		k562lmepdf.to_csv('K562.psi.pval.txt', sep = '\t', header = True, index = False)

	elif args.mode == 'pval_twolines':
		pvaluedf = getlmepdf_twocelllines(psidf)
		pvaluedf.to_csv('HepG2andK562.psi.pval.txt', sep = '\t', header = True, index = False)

