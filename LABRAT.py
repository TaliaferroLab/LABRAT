#!/usr/bin/env python

import gffutils
import os
import sys
import itertools
import gzip
import numpy as np
from Bio import SeqIO
import argparse
import subprocess
from functools import reduce
import pandas as pd
from collections import OrderedDict
import math
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy.stats.distributions import chi2
import warnings
import time

#For every gene, the goal is to define what the overall usage of proximal/distal polyA sites is.  This is done by defining 
#a "psi" value for every gene.  For each transcript in a gene, the "terminal fragment" (TF) is the last two
#exons of the transcript concatenated together.  For each transcript, a position factor (PF) is calculated as (m / n + 1; will range 
#between 0 and 1).
#The tpm for its TF is multiplied by its PF.  This is then done for all TFs and their scaled values are summed and divided
#by the unscaled sum of tpms for all TFs to give psi, which will also range between 0 (exculsive usage of the proximal polyA site)
#and 1 (exclusive usage of the distal polyA site).  

#Does a gene pass filters?
def genefilters(gene, db):
	proteincoding = False

	if 'protein_coding' in gene.attributes['gene_type']:
		proteincoding = True

	if proteincoding == True:
		return True
	else:
		return False

#Does a transcript pass filters?
def transcriptfilters(transcript, db):
	exonnumberpass = False
	TFlengthpass = False
	proteincoding = False 
	mrnaendpass = False
	#How many exons does it have
	if len(list(db.children(transcript, featuretype = 'exon'))) >= 2:
		exonnumberpass = True
	else:
		return False
	
	#What is the length of the terminal fragment
	exons = []
	if transcript.strand == '+':
		for exon in db.children(transcript, featuretype = 'exon', order_by = 'start'):
			exons.append([exon.start, exon.end + 1])
	elif transcript.strand == '-':
		for exon in db.children(transcript, featuretype = 'exon', order_by = 'start', reverse = True):
			exons.append([exon.start, exon.end + 1])
	penultimateexonlen = len(range(exons[-2][0], exons[-2][1]))
	lastexonlen = len(range(exons[-1][0], exons[-1][1]))
	TFlength = penultimateexonlen + lastexonlen
	if TFlength > 200:
		TFlengthpass = True

	#Is this transcript protein coding
	if 'protein_coding' in transcript.attributes['transcript_type']:
		proteincoding = True

	#Are we confident in the 3' end of this mrnaendpass
	if 'tag' not in transcript.attributes or 'mRNA_end_NF' not in transcript.attributes['tag']:
		mrnaendpass = True

	if exonnumberpass and TFlengthpass and proteincoding and mrnaendpass:
		return True
	else:
		return False



#Given an annotation in gff format, get the position factors for all transcripts.
#It merges any two transcript ends that are less than <lengthfilter> away from each other into a single end.
#This is so that you dont end up with unique regions that are like 4 nt long.
#They might causes issues when it comes to counting kmers or reads that map to a given region.

def getpositionfactors(gff, lengthfilter):
	lengthfilter = int(lengthfilter)
	genecount = 0
	txends = {} #{ENSMUSG : [strand, [list of distinct transcript end coords]]}
	posfactors = {} #{ENSMUSG : {ENSMUST : positionfactor}}

	#Make gff database
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	#Get number of distinct transcript ends for each gene
	genes = db.features_of_type('gene')
	for gene in genes:
		#Only protein coding genes
		passgenefilters = genefilters(gene, db)
		if passgenefilters == False:
			continue
		genename = str(gene.id).replace('gene:', '')
		ends = []
		if gene.strand == '+':
			for transcript in db.children(gene, featuretype = 'transcript', level = 1, order_by = 'end'):
				#Skip transcripts that do not pass filters
				passtranscriptfilters = transcriptfilters(transcript, db)
				if passtranscriptfilters == False:
					continue
				if transcript.end not in ends:
					ends.append(transcript.end)
		elif gene.strand == '-':
			for transcript in db.children(gene, featuretype = 'transcript', level = 1, order_by = 'start', reverse = True):
				#Skip transcripts that do not pass filters
				passtranscriptfilters = transcriptfilters(transcript, db)
				if passtranscriptfilters == False:
					continue
				if transcript.start not in ends:
					ends.append(transcript.start)

		if ends: #Sometimes there are no 'transcripts' for a gene, like with pseudogenes, etc.
			txends[genename] = [gene.strand, ends]

	#Sort transcript end coords
	s_txends = {} #{ENSMUSG : [sorted (most upstream to most downstream) tx end coords]}
	for gene in txends:
		strand = txends[gene][0]
		coords = txends[gene][1]
		if strand == '+':
			sortedcoords = sorted(coords)
		elif strand == '-':
			sortedcoords = sorted(coords, reverse = True)
		s_txends[gene] = sortedcoords

	#Get m values (the numerator of the position factor fraction), combining an end that is less than <lengthfilter> nt away from
	#the previous utr into the same m value as the previous utr

	mvalues = {} #{ENSMUSG : {txendcoord : mvalue}}
	for gene in s_txends:
		mvalues[gene] = {}
		currentendcoord = s_txends[gene][0]
		currentmvalue = 0
		mvalues[gene][currentendcoord] = 0 #the first one has to have m = 0
		for endcoord in s_txends[gene][1:]:
			#If this endcoord is too close to the last one
			if abs(endcoord - currentendcoord) <= lengthfilter:
				#this end gets the current m value
				mvalues[gene][endcoord] = currentmvalue
				#we stay on this m value for the next round
				currentmvalue = currentmvalue
				#update currentendcoord
				currentendcoord = endcoord
			#If this endcoord is sufficiently far away from the last one
			elif abs(endcoord - currentendcoord) > lengthfilter:
				#this end coord gets the next m value
				mvalues[gene][endcoord] = currentmvalue + 1
				#we move on to the next m value for the next round
				currentmvalue = currentmvalue + 1
				#update currentendcoord
				currentendcoord = endcoord


	#Figure out postion scaling factor for each transcript (position / (number of total positions - 1)) (m / (n - 1))
	genes = db.features_of_type('gene')
	for gene in genes:
		genecount +=1
		genename = str(gene.id).replace('gene:', '')
		#Only protein coding genes
		passgenefilters = genefilters(gene, db)
		if passgenefilters == False:
			continue
		#If this gene isnt in mvalues or there is only one m value for the entire gene, skip it
		if genename not in mvalues:
			continue
		if len(set(mvalues[genename].values())) == 1:
			continue
		#Get number of different m values for this gene
		n = len(set(mvalues[genename].values()))
		posfactors[genename] = {}
		for transcript in db.children(gene, featuretype = 'transcript', level = 1, order_by = 'end'):
			#Skip transcripts that do not pass filters
			passtranscriptfilters = transcriptfilters(transcript, db)
			if passtranscriptfilters == False:
				continue
			txname = str(transcript.id).replace('transcript:', '')
			if gene.strand == '+':
				m = mvalues[genename][transcript.end]
			elif gene.strand == '-':
				m = mvalues[genename][transcript.start]
			posfactor = m / float(n - 1)
			posfactors[genename][txname] = posfactor

	#Output file of the number of posfactors for each gene
	with open('numberofposfactors.txt', 'w') as outfh:
		outfh.write(('\t').join(['Gene', 'numberofposfactors', 'txids', 'interpolyAdist']) + '\n')
		for gene in posfactors: #{ENSMUSG : {ENSMUST : positionfactor}}
			pfs = []
			for tx in posfactors[gene]:
				pfs.append(posfactors[gene][tx])
			pfs = list(set(pfs))
			
			#write distance between polyA sites for those genes that only have 2 pfs
			if len(pfs) == 2:
				g = db[gene]
				for tx in posfactors[gene]:
					if posfactors[gene][tx] == 0:
						txpf1 = tx
					elif posfactors[gene][tx] == 1:
						txpf2 = tx
				t1 = db[txpf1]
				t2 = db[txpf2]
				if g.strand == '+':
					interpolyAdist = t2.end - t1.end
				elif g.strand == '-':
					interpolyAdist = t1.start - t2.start
			elif len(pfs) != 2:
				interpolyAdist = 'NA'

			#Get list of txs that belong to each pf
			txids = {} #{positionfactor : [list of transcriptIDs]}
			for pf in sorted(pfs):
				txids[pf] = []
				for tx in posfactors[gene]:
					if posfactors[gene][tx] == pf:
						txids[float(pf)].append(tx)

			alltxids = []
			for pf in sorted(list(txids.keys())):
				alltxids.append((',').join(txids[pf]))
			outfh.write(('\t').join([gene, str(len(pfs)), ('_').join(alltxids), str(interpolyAdist)]) + '\n')


	return posfactors


#Make a fasta file containing the "terminal fragments" of all transcripts
def makeTFfasta(gff, genomefasta, lasttwoexons, librarytype):
	TFs = {} #{transcriptid: [chrm, strand, [next to last exon1start, next to last exon1stop], [last exon2start, last exon2stop]]}
	
	#Make gff database
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	print('Indexing genome sequence...')
	if os.path.basename(genomefasta).endswith('.gz'):
		seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta, 'rt'), 'fasta'))
	else:
		seq_dict = SeqIO.to_dict(SeqIO.parse(genomefasta, 'fasta'))
	print('Done indexing!')

	genes = db.features_of_type('gene')

	#Get last two exons
	genecounter = 0
	for gene in genes:
		genecounter +=1
		if genecounter % 5000 == 0:
			print('Gene {0}...'.format(genecounter))
		#Only protein coding genes
		passgenefilters = genefilters(gene, db)
		if passgenefilters == False:
			continue
		for transcript in db.children(gene, featuretype = 'transcript'):
			#Skip transcripts that do not pass filters
			passtranscriptfilters = transcriptfilters(transcript, db)
			if passtranscriptfilters == False:
				continue
			exons = [] #[[exon1start, exon1stop], ..., [exonnstart, exonnstop]]
			txname = str(transcript.id).replace('transcript:', '').split('.')[0]
			if gene.strand == '+':
				for exon in db.children(transcript, featuretype = 'exon', order_by = 'start'):
					exons.append([exon.start, exon.end])
			elif gene.strand == '-':
				for exon in db.children(transcript, featuretype = 'exon', order_by = 'start', reverse = True):
					exons.append([exon.start, exon.end])

			#For last two exons
			if lasttwoexons:
				TFs[txname] = [transcript.chrom, transcript.strand, exons[-2], exons[-1]]

			#For all exons
			elif not lasttwoexons:
				TFs[txname] = [transcript.chrom, transcript.strand]
				for exon in exons:
					TFs[txname].append(exon)

	
	#Get sequences of TFs
	if lasttwoexons:
		with open('TFseqs.fasta', 'w') as outfh:
			for TF in TFs:
				chrm, strand, exon1, exon2 = TFs[TF]
				if strand == '+':
					exon1seq = seq_dict[chrm].seq[exon1[0] - 1:exon1[1]].upper()
					exon2seq = seq_dict[chrm].seq[exon2[0] - 1:exon2[1]].upper()
				elif strand == '-':
					exon1seq = seq_dict[chrm].seq[exon1[0] - 1:exon1[1]].reverse_complement().upper()
					exon2seq = seq_dict[chrm].seq[exon2[0] - 1:exon2[1]].reverse_complement().upper()

				TFseq = exon1seq + exon2seq

				#If this is 3' end data we are only going to quantify the last 300 nt of every transcript
				if librarytype == '3pseq':
					if len(TFseq) < 300:
						outfh.write('>' + TF + '\n' + str(TFseq) + '\n')
					else:
						outfh.write('>' + TF + '\n' + str(TFseq)[-300:] + '\n')

				elif librarytype == 'RNAseq':
					outfh.write('>' + TF + '\n' + str(TFseq) + '\n')

	#Get sequences of all exons
	elif not lasttwoexons:
		with open('wholetranscriptseqs.fasta', 'w') as outfh:
			for TF in TFs:
				chrm, strand = TFs[TF][0], TFs[TF][1]
				exons = TFs[TF][2:]
				seq = ''
				for exon in exons:
					if strand == '+':
						exonseq = seq_dict[chrm].seq[exon[0] - 1:exon[1]].upper()
					elif strand == '-':
						exonseq = seq_dict[chrm].seq[exon[0] - 1:exon[1]].reverse_complement().upper()

					seq += exonseq

				#If this is 3' end data we are only going to quantify the last 300 nt of every transcript
				if librarytype == '3pseq':
					if len(seq) < 300:
						outfh.write('>' + TF + '\n' + str(seq) + '\n')
					else:
						outfh.write('>' + TF + '\n' + str(seq)[-300:] + '\n')

				elif librarytype == 'RNAseq':
					outfh.write('>' + TF + '\n' + str(seq) + '\n')


def runSalmon(threads, reads1, reads2, samplename):
	#paired end
	if reads2:
		print('Running salmon for {0}...'.format(samplename))
		command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--seqBias', '--gcBias', '-1', reads1, '-2', reads2, '-o', samplename, '--index', 'txfasta.idx', '--validateMappings']

	#Single end
	elif not reads2:
		fldMean = '250' #fragment length distribution mean
		fldSD = '20' #fragment length distribution standard deviation	
		command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--fldMean', fldMean, '--fldSD', fldSD, '--seqBias', '-r', reads1, '-o', samplename, '--index', 'txfasta.idx', '--validateMappings']

	subprocess.call(command)


def calculatepsi(positionfactors, salmondir, librarytype):
	txtpms = {} #{transcriptid : tpm} this is actually tpm or counts depending on library type
	genetpms = {} #{geneid : [transcript tpms]}
	posfactorgenetpms = {} #{geneid : [transcripttpms scaled by posfactor]}
	psis = {} #{geneid : psi}

	#Read in salmonout file
	#salmondir is the salmon output directory
	print('Calculating psi values for {0}...'.format(salmondir))
	salmonoutfile = os.path.join(os.path.abspath(salmondir), 'quant.sf')
	with open(salmonoutfile, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			#Skip header
			if line[0] == 'Name':
				continue
			transcriptid = str(line[0]).split('.')[0]
			if librarytype == 'RNAseq':
				tpm = float(line[3])
			#For 3p seq data, we don't need to worry about length normalization, and actually we shouldn't
			#because it would unfairly penalize long transcripts. So take the counts from the salmon quant
			#instead of the tpm.
			elif librarytype == '3pseq':
				tpm = float(line[4]) #this is really counts, but it's simpler to keep the variable name the same
			txtpms[transcriptid] = tpm

	#Put the transcript tpms for every gene together
	for gene in positionfactors:
		genetpms[gene] = []
		posfactorgenetpms[gene] = []
		for transcript in positionfactors[gene]:
			txtpm = txtpms[transcript.split('.')[0]]
			genetpms[gene].append(txtpm)
			posfactor = positionfactors[gene][transcript]
			posfactorgenetpms[gene].append(txtpm * posfactor)

	#Calculate psi for every gene (can apply filters for min gene tpm here)
	for gene in posfactorgenetpms:
		scaledtpm = sum(posfactorgenetpms[gene])
		totaltpm = sum(genetpms[gene])
		if librarytype == 'RNAseq':
			if totaltpm >= 5: #gene level tpm filter
				psi = scaledtpm / float(totaltpm)
				psis[gene] = float(format(psi, '.3f'))
			else:
				#psis[gene] = 'NA'
				psis[gene] = np.nan
		elif librarytype == '3pseq':
			if totaltpm >= 100: #gene-level count filter
				psi = scaledtpm / float(totaltpm)
				psis[gene] = float(format(psi, '.3f'))
			else:
				psis[gene] = np.nan

	return psis


def getdpsis_covariate(psifile, samp_conds_file, conditionA, conditionB):
	#Given a table of psis, calculate dpsis and LME-based p values
	deltapsidict = OrderedDict() #{genename : pvalue} ordered so it's easy to match it up with df
	pvaluedict = OrderedDict() #{genename : pvalue} ordered so it's easy to match it up with q values

	psidf = pd.read_csv(psifile, sep = '\t', header = 0, index_col = False)

	sampconddf = pd.read_csv(samp_conds_file, sep = '\t', index_col = False, header = 0)
	colnames = list(sampconddf.columns)
	covariate_columns = [c for c in colnames if 'covariate' in c]

	condasamps = sampconddf.loc[sampconddf['condition'] == conditionA, 'sample'].tolist()
	condbsamps = sampconddf.loc[sampconddf['condition'] == conditionB, 'sample'].tolist()

	#Reorder columns of psidf
	colorder = ['Gene'] + condasamps + condbsamps + ['genetype']
	psidf = psidf[colorder]

	print('Condition A samples: ' + (', ').join(condasamps))
	print('Condition B samples: ' + (', ').join(condbsamps))

	#Store relationships of conditions and the samples in that condition
	#It's important that this dictionary be ordered because we are going to be iterating through it
	samp_conds = OrderedDict({'cond1' : condasamps, 'cond2' : condbsamps})

	#Get a list of all samples
	samps = []
	for cond in samp_conds:
		samps += samp_conds[cond]

	#Iterate through rows, making a dictionary from every row, turning it into a dataframe, then calculating p value
	genecounter = 0
	for index, row in psidf.iterrows():
		genecounter +=1
		if genecounter % 1000 == 0:
			print('Calculating pvalue for gene {0}...'.format(genecounter))
		
		d = {}
		d['Gene'] = [row['Gene']] * len(samps)
		d['variable'] = samps
		
		values = [] #psi values
		for cond in samp_conds:
			for sample in samp_conds[cond]:
				value = row[sample]
				values.append(value)
		d['value'] = values

		if covariate_columns:
			for covcol in covariate_columns:
				covs= [] #classifications for this covariate 
				for samp in samps:
					cov = sampconddf.loc[sampconddf['sample'] == samp, covcol].tolist()[0]
					covs.append(cov)
				d[covcol] = covs

			covariatestring = '+'.join(covariate_columns)

		#If there is an NA psi value, we are not going to calculate a pvalue for this gene
		p = None
		if True in np.isnan(values):
			p = np.nan
		
		conds = []
		for cond in samp_conds:
			conds += [cond] * len(samp_conds[cond])
		cond1s = []
		cond2s = []
		for cond in conds:
			if cond == 'cond1':
				cond1s.append(1)
				cond2s.append(0)
			elif cond == 'cond2':
				cond1s.append(0)
				cond2s.append(1)
		d['cond1'] = cond1s #e.g. [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
		d['cond2'] = cond2s #e.g. [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]

		d['samples'] = [x + 1 for x in range(len(samps))]

		#Turn this dictionary into a DataFrame
		rowdf = pd.DataFrame.from_dict(d)

		#delta psi is difference between mean psi of two conditions (cond2 - cond1)
		cond1meanpsi = float(format(np.mean(rowdf.query('cond1 == 1').value.dropna()), '.3f'))
		cond2meanpsi = float(format(np.mean(rowdf.query('cond2 == 1').value.dropna()), '.3f'))
		deltapsi = cond2meanpsi - cond1meanpsi
		deltapsi = float(format(deltapsi, '.3f'))
		deltapsidict[row['Gene']] = deltapsi

		covariatestring = '+'.join(covariate_columns)

		#Get LME pvalue, but only if we haven't already determined that the pvalue is NA because we are missing one or more psi values
		#Lots of warnings about convergence, etc. Suppress them.
		if not p:
			with warnings.catch_warnings():
				warnings.filterwarnings('ignore')

				#So apparently, some combinations of psi values will give nan p values due to a LinAlgError that arises from a singular
				#hessian matrix during the fit of the model.  However, the equivalent code in R (nlme::lme) never gives this error, even with
				#the same data. It's not clear from just looking at the psi values why this is.  However, I found that by varying the 
				#start_params in the fit, this can be avoided. If this is done, the resulting p value always matches what is given in R.
				#Further, the p value is the same regardless of the start_param.
				#But it's not clear to me why changing the start_param matters, or what the default is here or with nlme.
				#So let's try a few starting paramters.  Regardless, this seems to affect a small number of genes (<1%), and it is causing 
				#false negatives because genes that should get p values (may or may not be sig) are getting NA.
				possible_start_params = [0, 0, 1, -1, 2, -2]
				numberoftries = -1
				for param in possible_start_params:
					#if we already have a pvalue, don't try again
					if p != None and not np.isnan(p):
						break
					#First time through, numberoftries = 0, and we are just using a placeholder startparam (0) here because we aren't even using it.
					#Gonna use whatever the default is
					numberoftries +=1
					try:
						#actual model
						if covariate_columns:
							md = smf.mixedlm('value ~ cond1' + '+' + covariatestring, data = rowdf, groups = 'samples', missing = 'drop')
						else:
							md = smf.mixedlm('value ~ cond1', data = rowdf, groups = 'samples', missing = 'drop')
						if numberoftries == 0:
							mdf = md.fit(reml = False) #REML needs to be false in order to use log-likelihood for pvalue calculation
						elif numberoftries > 0:
							mdf = md.fit(reml = False, start_params = [param])

						#null model
						if covariate_columns:
							nullmd = smf.mixedlm('value ~ ' + covariatestring, data = rowdf, groups = 'samples', missing = 'drop')
						else:
							nullmd = smf.mixedlm('value ~ 1', data = rowdf, groups = 'samples', missing = 'drop')
						if numberoftries == 0:
							nullmdf = nullmd.fit(reml = False)
						elif numberoftries > 0:
							nullmdf = nullmd.fit(reml = False, start_params = [param])

						#Likelihood ratio
						LR = 2 * (mdf.llf - nullmdf.llf)
						p = chi2.sf(LR, df = 1)

					#These exceptions are needed to catch cases where either all psi values are nan (valueerror) or all psi values for one condition are nan (linalgerror)
					except (ValueError, np.linalg.LinAlgError):
						p = np.nan

		pvaluedict[row['Gene']] = float('{:.2e}'.format(p))

	#Correct pvalues using BH method, but only using pvalues that are not NA
	pvalues = list(pvaluedict.values())
	pvaluesformultitest = [pval for pval in pvalues if str(pval) != 'nan']
	fdrs = list(multipletests(pvaluesformultitest, method = 'fdr_bh')[1])
	fdrs = [float('{:.2e}'.format(fdr)) for fdr in fdrs]

	#Now we need to incorporate the places where p = NA into the list of FDRs (also as NA)
	fdrswithnas = []
	fdrindex = 0
	for pvalue in pvalues:
		#print(pvalue)
		if str(pvalue) != 'nan':
			fdrswithnas.append(fdrs[fdrindex])
			fdrindex +=1
		elif str(pvalue) == 'nan':
			fdrswithnas.append(np.nan)

	#Add deltapsis, pvalues, and FDRs to df
	deltapsis = list(deltapsidict.values())
	psidf = psidf.assign(deltapsi = deltapsis)
	psidf = psidf.assign(pval = pvalues)
	psidf = psidf.assign(FDR = fdrswithnas)

	#Write df
	fn = os.path.abspath(psifile) + '.pval'
	psidf.to_csv(fn, sep = '\t', index = False, float_format = '%.3g', na_rep = 'NA')		

def getdpsis(psifile, samp_conds_file):
	#Old function, superceded by getdpsis_covariate
	#Given a table of psis, calculate dpsis and LME-based p values
	deltapsidict = OrderedDict() #{genename : pvalue} ordered so it's easy to match it up with df
	pvaluedict = OrderedDict() #{genename : pvalue} ordered so it's easy to match it up with q values

	psidf = pd.read_csv(psifile, sep = '\t', header = 0, index_col = False)

	sampconddf = pd.read_csv(samp_conds_file, sep = '\t', index_col = False, header = None)
	cond1samps = sampconddf[0].dropna().tolist()
	cond2samps = sampconddf[1].dropna().tolist()

	#Reorder columns of psidf
	colorder = ['Gene'] + cond1samps + cond2samps + ['genetype']
	psidf = psidf[colorder]

	print('Condition 1 samples: ' + (', ').join(cond1samps))
	print('Condition 2 samples: ' + (', ').join(cond2samps))

	#Store relationships of conditions and the samples in that condition
	#It's important that this dictionary be ordered because we are going to be iterating through it
	samp_conds = OrderedDict({'cond1' : cond1samps, 'cond2' : cond2samps})

	#Get a list of all samples
	samps = []
	for cond in samp_conds:
		samps += samp_conds[cond]

	#Iterate through rows, making a dictionary from every row, turning it into a dataframe, then calculating p value
	genecounter = 0
	for index, row in psidf.iterrows():
		genecounter +=1
		if genecounter % 1000 == 0:
			print('Calculating pvalue for gene {0}...'.format(genecounter))
		
		d = {}
		d['Gene'] = [row['Gene']] * len(samps)
		d['variable'] = samps
		
		values = [] #psi values
		for cond in samp_conds:
			for sample in samp_conds[cond]:
				value = row[sample]
				values.append(value)
		d['value'] = values

		#If there is an NA psi value, we are not going to calculate a pvalue for this gene
		p = None
		if True in np.isnan(values):
			p = np.nan
		
		conds = []
		for cond in samp_conds:
			conds += [cond] * len(samp_conds[cond])
		cond1s = []
		cond2s = []
		for cond in conds:
			if cond == 'cond1':
				cond1s.append(1)
				cond2s.append(0)
			elif cond == 'cond2':
				cond1s.append(0)
				cond2s.append(1)
		d['cond1'] = cond1s #e.g. [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
		d['cond2'] = cond2s #e.g. [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]

		d['samples'] = [x + 1 for x in range(len(samps))]

		#Turn this dictionary into a DataFrame
		rowdf = pd.DataFrame.from_dict(d)

		#delta psi is difference between mean psi of two conditions (cond2 - cond1)
		cond1meanpsi = float(format(np.mean(rowdf.query('cond1 == 1').value.dropna()), '.3f'))
		cond2meanpsi = float(format(np.mean(rowdf.query('cond2 == 1').value.dropna()), '.3f'))
		deltapsi = cond2meanpsi - cond1meanpsi
		deltapsi = float(format(deltapsi, '.3f'))
		deltapsidict[row['Gene']] = deltapsi

		#Get LME pvalue, but only if we haven't already determined that the pvalue is NA because we are missing one or more psi values
		#Lots of warnings about convergence, etc. Suppress them.
		if not p:
			with warnings.catch_warnings():
				warnings.filterwarnings('ignore')

				try:
					#actual model
					md = smf.mixedlm('value ~ cond1', data = rowdf, groups = 'cond1', missing = 'drop')
					mdf = md.fit(reml = False) #REML needs to be false in order to use log-likelihood for pvalue calculation

					#null model
					nullmd = smf.mixedlm('value ~ 1', data = rowdf, groups = 'samples', missing = 'drop')
					nullmdf = nullmd.fit(reml = False)

					#Likelihood ratio
					LR = 2 * (mdf.llf - nullmdf.llf)
					p = chi2.sf(LR, df = 1)

				#These exceptions are needed to catch cases where either all psi values are nan (valueerror) or all psi values for one condition are nan (linalgerror)
				except (ValueError, np.linalg.LinAlgError):
					p = np.nan

		pvaluedict[row['Gene']] = float('{:.2e}'.format(p))

	#Correct pvalues using BH method, but only using pvalues that are not NA
	pvalues = list(pvaluedict.values())
	pvaluesformultitest = [pval for pval in pvalues if str(pval) != 'nan']
	fdrs = list(multipletests(pvaluesformultitest, method = 'fdr_bh')[1])
	fdrs = [float('{:.2e}'.format(fdr)) for fdr in fdrs]

	#Now we need to incorporate the places where p = NA into the list of FDRs (also as NA)
	fdrswithnas = []
	fdrindex = 0
	for pvalue in pvalues:
		#print(pvalue)
		if str(pvalue) != 'nan':
			fdrswithnas.append(fdrs[fdrindex])
			fdrindex +=1
		elif str(pvalue) == 'nan':
			fdrswithnas.append(np.nan)

	#Add deltapsis, pvalues, and FDRs to df
	deltapsis = list(deltapsidict.values())
	psidf = psidf.assign(deltapsi = deltapsis)
	psidf = psidf.assign(pval = pvalues)
	psidf = psidf.assign(FDR = fdrswithnas)

	#Write df
	fn = os.path.abspath(psifile) + '.pval'
	psidf.to_csv(fn, sep = '\t', index = False, float_format = '%.3g', na_rep = 'NA')


#We want to know if the alternative 3' ends for this gene are either all of ALE forms or all of tandem UTR forms.
#If it's a mixture, we can't really deal with that cleanly, so we are going to ignore it
def getexoniccoords(posfactors, gff):
	#Take posfactors dictionary
	#Take the transcripts for each gene, get their exonic coords.
	#if txend of a transcript is exonic for every transcript that has a higher PF, then this gene is all TUTR
	#if all txends are not exonic in any other transcript, then this gene is all ALE
	#if neither of these are true, then the gene is somehow mixed
	#posfactors = {ENSMUSG : {ENSMUST : positionfactor}}
	#Make gff database
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	exoniccoords = {} #{geneid : {txid : [positionfactor, [set of exonic coords]]}}
	for gene in posfactors:
		txs = posfactors[gene].keys()
		geneexons = {} #{txid : [positionfactor, [set of exonic coords]]}
		for tx in txs:
			txexons = []
			pf = posfactors[gene][tx]
			tx = db[tx]
			if tx.strand == '+':
				for exon in db.children(tx, featuretype = 'exon', order_by = 'start'):
					txexons += list(range(exon.start, exon.end + 1))
			elif tx.strand == '-':
				for exon in db.children(tx, featuretype = 'exon', order_by = 'start', reverse = True):
					txexons += list(range(exon.start, exon.end + 1))

			geneexons[tx.id] = [pf, set(txexons)]

		exoniccoords[gene] = geneexons

	return exoniccoords



#Could this gene contain only ALEs for 3' ends?
def isitALE(geneexoniccoords, db):
	couldbeALE = True
	for txid in geneexoniccoords:
		tx = db[txid]
		if tx.strand == '+':
			txend = int(tx.end)
		elif tx.strand == '-':
			txend = int(tx.start)

		#See if this txend is exonic for any other transcript with a different positionfactor
		pf = geneexoniccoords[txid][0]
		for othertx in geneexoniccoords:
			othertxpf = geneexoniccoords[othertx][0]
			#If it has the same PF, move on
			if pf == othertxpf:
				continue
			othertxexoncoords = geneexoniccoords[othertx][1]
			if txend in othertxexoncoords:
				couldbeALE = False

	return couldbeALE

#Could this gene contain only TUTRs for 3' ends?
def isitTUTR(geneexoniccoords, db):
	couldbeTUTR = True
	for txid in geneexoniccoords:
		tx = db[txid]
		if tx.strand == '+':
			txend = int(tx.end)
		elif tx.strand == '-':
			txend = int(tx.start)

		#See if this txend is exonic for every transcript that has a higher PF
		pf = geneexoniccoords[txid][0]
		for othertx in geneexoniccoords:
			othertxpf = geneexoniccoords[othertx][0]
			#If the transcript of interest has a greater of equal PF to the transcript we are comparing to, move on
			if pf >= othertxpf:
				continue
			othertxexoncoords = geneexoniccoords[othertx][1]
			if txend not in othertxexoncoords:
				couldbeTUTR = False

	return couldbeTUTR


def classifygenes(exoniccoords, gff):
	#exoniccoords = {} #{geneid : {txid : [positionfactor, [set of exonic coords]]}}
	genetypes = {} #{geneid : type of 3' UTR end}
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	for gene in exoniccoords:
		geneexoniccoords = exoniccoords[gene]
		coulditbeALE = isitALE(geneexoniccoords, db)
		coulditbeTUTR = isitTUTR(geneexoniccoords, db)
		if coulditbeALE == True and coulditbeTUTR == False:
			genetype = 'ALE'
		elif coulditbeALE == False and coulditbeTUTR == True:
			genetype = 'TUTR'
		elif coulditbeALE == False and coulditbeTUTR == False:
			genetype = 'mixed'
		elif coulditbeALE == True and coulditbeTUTR == True:
			genetype = 'ERROR'

		genetypes[gene] = genetype

	return genetypes


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--mode', type = str, choices = ['makeTFfasta', 'runSalmon', 'calculatepsi', 'test'])
	parser.add_argument('--librarytype', type = str, choices = ['RNAseq', '3pseq'], help = 'Is this RNAseq data or 3\' seq data? Needed for makeTFfasta and calculatepsi.')
	parser.add_argument('--gff', type = str, help = 'GFF of transcript annotation. Needed for makeTFfasta and calculatepsi.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format. Needed for makeTFfasta.')
	parser.add_argument('--lasttwoexons', action = 'store_true', help = 'Used for makeTFfasta. Do you want to only use the last two exons?')
	parser.add_argument('--txfasta', type = str, help = 'Fasta file of sequences to quantitate with salmon. Often the output of makeTFfasta mode.')
	parser.add_argument('--reads1', type = str, help = 'Comma separated list of forward read fastq files. Needed for runSalmon.')
	parser.add_argument('--reads2', type = str, help = 'Comma separated list of reverse read fastq files. Needed for runSalmon. Omit for single end data.')
	parser.add_argument('--samplename', type = str, help = 'Comma separated list of samplenames.  Needed for runSalmon.')
	parser.add_argument('--threads', type = str, help = 'Number of threads to use.  Needed for runSalmon.')
	parser.add_argument('--salmondir', type = str, help = 'Salmon output directory. Needed for calculatepsi.')
	parser.add_argument('--sampconds', type = str, help = 'File relating sample names and conditions. See README for details.')
	parser.add_argument('--conditionA', type = str, help = 'Condition A. deltapsi will be calculated as B-A. Must be a value in the \'condition\' column of the sampconds file.')
	parser.add_argument('--conditionB', type = str, help = 'Condition B. deltapsi will be calculated as B-A. Must be a value in the \'condition\' column of the sampconds file.')
	args = parser.parse_args()


	if args.mode == 'makeTFfasta':
		if not args.gff or not args.genomefasta or not args.librarytype:
			print('You have not supplied all the required arguments! See the --help for more info.')
			sys.exit()
		makeTFfasta(args.gff, args.genomefasta, args.lasttwoexons, args.librarytype)

	elif args.mode == 'runSalmon':
		if not args.txfasta or not args.reads1 or not args.samplename or not args.threads:
			print('You have not supplied all the required arguments! See the --help for more info.')
			sys.exit()
		forreads = args.reads1.split(',')
		if args.reads2:
			revreads = args.reads2.split(',')
		samplenames = args.samplename.split(',')

		command = ['salmon', 'index', '-t', args.txfasta, '-i', 'txfasta.idx', '--type', 'quasi', '-k', '31', '--keepDuplicates']
		print('Indexing transcripts...')
		subprocess.call(command)
		print('Done indexing!')

		if args.reads2:
			if len(forreads) != len(revreads) or len(forreads) != len(samplenames) or len(revreads) != len(samplenames):
				print('ERROR: The number of forward read files, reverse read files, and sample names must match!')
				sys.exit()

		elif not args.reads2:
			if len(forreads) != len(samplenames):
				print('ERROR: The number of forward read files and sample names must match!')
				sys.exit()

		if args.reads2:
			for i in range(len(forreads)):
				freads = forreads[i]
				rreads = revreads[i]
				samplename = samplenames[i]
				runSalmon(args.threads, freads, rreads, samplename)

		elif not args.reads2:
			for i in range(len(forreads)):
				freads = forreads[i]
				samplename = samplenames[i]
				runSalmon(args.threads, freads, None, samplename)

	elif args.mode == 'calculatepsi':
		if not args.gff or not args.salmondir or not args.librarytype:
			print('You have not supplied all the required arguments! See the --help for more info.')
			sys.exit()
		print('Calculating position factors for every transcript...')
		positionfactors = getpositionfactors(args.gff, 25)
		print('Done with position factors!')
		salmondirs = [os.path.join(os.path.abspath(args.salmondir), d) for d in os.listdir(args.salmondir) if os.path.isdir(os.path.join(os.path.abspath(args.salmondir), d)) and d != 'txfasta.idx']
		psidfs = []
		for sd in salmondirs:
			samplename = os.path.basename(sd)
			
			#For ENCODE
			#samplename = os.path.basename(sd).split('_')
			#samplename = ('_rep').join([samplename[0], samplename[1]])
			
			psis = calculatepsi(positionfactors, sd, args.librarytype)
			psidf = pd.DataFrame.from_dict(psis, orient = 'index')
			psidf.reset_index(level = 0, inplace = True)
			psidf.columns = ['Gene', samplename]
			psidfs.append(psidf)
		bigpsidf = reduce(lambda x, y: pd.merge(x, y, on = 'Gene'), psidfs)
		
		#Add in genetypes (ALE or TUTR or mixed)
		exoniccoords = getexoniccoords(positionfactors, args.gff)
		genetypes = classifygenes(exoniccoords, args.gff)
		genetypedf = pd.DataFrame.from_dict(genetypes, orient = 'index')
		genetypedf.reset_index(level = 0, inplace = True)
		genetypedf.columns = ['Gene', 'genetype']
		finalpsidf = reduce(lambda x, y: pd.merge(x, y, on = 'Gene'), [bigpsidf, genetypedf])
		finalpsidf.to_csv('LABRAT.psis', sep = '\t', index = False, na_rep = 'NA')
		getdpsis_covariate('LABRAT.psis', args.sampconds, args.conditionA, args.conditionB)
