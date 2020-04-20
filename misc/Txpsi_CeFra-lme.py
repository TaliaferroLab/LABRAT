import gffutils
import os
import sys
import itertools
import gzip
import numpy as np
from Bio import SeqIO
import argparse
import subprocess
import pandas as pd
from collections import OrderedDict
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, Formula, FloatVector
from rpy2.rinterface import RRuntimeError
import math
from itertools import combinations
import statsmodels.api as sm

def getrowdict(index, row, samp_conds, fraca, fracb):
	#Get a list of all samples
	samps = []
	for cond in samp_conds:
		for libprep in samp_conds[cond]:
			samps += samp_conds[cond][libprep]

	d = {}
	d['Gene'] = [row['Gene']] * len(samps)
	d['variable'] = samps
	values = [] #psi values
	for cond in samp_conds:
		for libprep in samp_conds[cond]:
			for sample in samp_conds[cond][libprep]:
				value = row[sample]
				values.append(value)
	d['value'] = values

	conds = []
	for cond in samp_conds:
		conds += [cond] * sum(len(libprep) for libprep in samp_conds.itervalues())
	condas = []
	condbs = []
	for cond in conds:
		if cond == fraca:
			condas.append(1)
			condbs.append(0)
		elif cond == fracb:
			condas.append(0)
			condbs.append(1)
	d['conda'] = condas #e.g. [1, 1, 1, 1, 0, 0, 0, 0]
	d['condb'] = condbs #e.g. [0, 0, 0, 0, 1, 1, 1, 1]

	libpreps = []
	for cond in samp_conds:
		for libprep in samp_conds[cond]:
			libpreps += [libprep] * len(samp_conds[cond][libprep])
	polyAs = []
	ribodeps = []
	for libprep in libpreps:
		if libprep == 'polyA':
			polyAs.append(1)
			ribodeps.append(0)
		elif libprep == 'ribodep':
			polyAs.append(0)
			ribodeps.append(1)
	d['polyA'] = polyAs #e.g. [1, 1, 1, 1, 0, 0, 0, 0]
	d['ribodep'] = ribodeps #e.g. [0, 0, 0, 0, 1, 1, 1, 1]

	d['samples'] = [x + 1 for x in range(len(samps))]

	return d


def getdpsis(psifile):
	#Given a table of psis, calculate LME-based p values
	psidf = pd.read_table(psifile, sep = '\t', header = 0, index_col = False)

	pandas2ri.activate() #allow conversion between pandas dataframes and r dataframes

	#define R packages
	nlme = importr('nlme')
	base = importr('base')
	stats = importr('stats')
	qv = importr('qvalue')

	#define formulae
	fmla = Formula('value ~ 1 + conda + polyA')
	rndm = Formula('~ 1 | samples')
	nullfmla = Formula('value ~ 1 + polyA')
	nullrndm = Formula('~1 | samples')

	#Remove any gene that has a psi of NA in any sample
	psidf = psidf.dropna(axis=0)

	#Store relationships of conditions and the samples in that condition
	#It's important that this dictionary be ordered because we are going to be iterating through it
	fracs = ['cytosol', 'membrane', 'insoluble', 'total']

	for combination in combinations(fracs, 2):
		fraca = combination[0]
		fracb = combination[1]
		pvalues = []

		samp_conds = OrderedDict({fraca : {'polyA' : [fraca + '_polyA_Rep1', fraca + '_polyA_Rep2'], 
			'ribodep' : [fraca + '_ribodep_Rep1', fraca + '_ribodep_Rep2']},
			fracb : {'polyA' : [fracb + '_polyA_Rep1', fracb + '_polyA_Rep2'],
			'ribodep' : [fracb + '_ribodep_Rep1', fracb + '_ribodep_Rep2']}})

		#Get a list of all samples
		samps = []
		for cond in samp_conds:
			for libprep in samp_conds[cond]:
				samps += samp_conds[cond][libprep]

		#Iterate through rows, making a dictionary from every row, turning it into a dataframe, then calculating p value
		genecounter = 0
		for index, row in psidf.iterrows():
			genecounter +=1
			if genecounter % 1000 == 0:
				print 'Gene {0}...'.format(genecounter)

			d = getrowdict(index, row, samp_conds, fraca, fracb)

			#Turn this dictionary into a dataframe
			rowdf = pd.DataFrame.from_dict(d)

			#Get lme p value
			try:
				lm_alt = nlme.lme(fmla, random = rndm, data = rowdf, method = 'ML') #test
				lm_null = nlme.lme(nullfmla, random = nullrndm, data = rowdf, method = 'ML') #control
				logratio = (stats.logLik(lm_alt)[0] - stats.logLik(lm_null)[0]) * 2
				pvalue = stats.pchisq(logratio, df = 1, lower_tail = False)[0]
				#format decimal
				pvalue = float('{:.2e}'.format(pvalue))
			except RRuntimeError:
				print 'RRuntime error for {0}!'.format(row['Gene'])
				pvalue = 1.0

			pvalues.append(pvalue)

		#Turn list of pvalues into qvalues
		pvec = FloatVector(pvalues)
		#Get qvalues object
		qobj = qv.qvalue(p = pvec)
		#qvalues are index 2 of qvalue object
		qvalues = list(qobj[2])
		#format decimal
		qvalues = [float('{:.2e}'.format(qvalue)) for qvalue in qvalues]

		#Add pvalues and qvalues to df
		pvalcolname = fraca + '_vs_' + fracb + '_pval'
		qvalcolname = fraca + '_vs_' + fracb + '_qval'
		psidf[pvalcolname] = pvalues
		psidf[qvalcolname] = qvalues
	psidf.to_csv('Drosophilapsi.pval.txt', sep = '\t', header = True, index = False)

getdpsis(sys.argv[1])