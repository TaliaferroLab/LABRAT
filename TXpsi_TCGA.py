#Given a directory of TCGA outputs, get psis and pvalues

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

#Does a gene pass filters?
def genefilters(gene, db):
	return True #not using this filter at all
	proteincoding = True #don't want this filter

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
	proteincoding = True #don't want this filter
	mrnaendpass = False #don't want this filter (for drosophila samples)
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
#THIS IS AN UPDATED GETPOSITIONFACTORS FUNCTION
#It merges any two transcript ends that are less than <lengthfilter> away from each other into a single end.
#This is so that you dont end up with unique regions that are like 4 nt long.
#They might causes issues when it comes to counting kmers or reads that map to a given region.

def getpositionfactors(gff, lengthfilter):
	lengthfilter = int(lengthfilter)
	genecount = 0
	txends = {} #{ENSMUSG : [strand, [list of distinct transcript end coords]]}
	posfactors = {} #{ENSMUSG : {ENSMUST : positionfactor}}

	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

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
				#m = possibleends.index(transcript.end)
				m = mvalues[genename][transcript.end]
			elif gene.strand == '-':
				#m = possibleends.index(transcript.start)
				m = mvalues[genename][transcript.start]
			posfactor = m / float(n - 1)
			posfactors[genename][txname] = posfactor
			#posfactors[genename][txname] = m

	return posfactors


#Calculate psi values in a single salmon output.
#Used by calculatepsi_multisample
def calculatepsi_singlesample(positionfactors, salmondir):
	txtpms = {} #{transcriptid : tpm}
	genetpms = {} #{geneid : [transcript tpms]}
	posfactorgenetpms = {} #{geneid : [transcripttpms scaled by posfactor]}
	psis = {} #{geneid : psi}

	#Read in salmonout file
	#salmondir is the salmon output directory
	print 'Calculating psi values for {0}...'.format(salmondir)
	salmonoutfile = os.path.join(os.path.abspath(salmondir), 'quant.sf')
	with open(salmonoutfile, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			#Skip header
			if line[0] == 'Name':
				continue
			transcriptid = str(line[0])
			tpm = float(line[3])
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
		if totaltpm > 5:
			psi = scaledtpm / float(totaltpm)
			psis[gene] = psi
		else:
			psis[gene] = 'NA'

	return psis

#Calculate psis for a samples in a directory containing salmon outputs
#Output a psi value table
def calculatepsi_multisample(gff, salmondir):
	print 'Calculating position factors for every transcript...'
	positionfactors = getpositionfactors(gff, 25)
	print 'Done with position factors!'
	salmondirs = [os.path.join(os.path.abspath(salmondir), d) for d in os.listdir(salmondir) if os.path.isdir(os.path.join(os.path.abspath(salmondir), d)) and d.endswith('_salmon')]
	psidfs = []
	samplenames = []
	for sd in salmondirs:
		samplename = os.path.basename(sd).split('_')[0]
		samplenames.append(samplename)
		psis = calculatepsi_singlesample(positionfactors, sd)
		psidf = pd.DataFrame.from_dict(psis, orient = 'index')
		psidf.reset_index(level = 0, inplace = True)
		psidf.columns = ['Gene', samplename]
		psidfs.append(psidf)
	bigpsidf = reduce(lambda x, y: pd.merge(x, y, on = 'Gene'), psidfs)
	bigpsidf.to_csv('{0}psis.txt'.format(os.path.basename(os.path.abspath(salmondir))), sep = '\t', index = False)
	#why the fuck we have to write the table then reread it, I don't know, but it won't work unless you do.
	bigpsidf = pd.read_table('{0}psis.txt'.format(os.path.basename(os.path.abspath(salmondir))), sep = '\t', header = 0, index_col = False)

	#Calculate dpsis
	samplebasenames = [s.replace('-norm', '').replace('-tumor', '') for s in samplenames]
	samplebasenames = list(set(samplebasenames))
	for samplebasename in samplebasenames:
		tumorsample = samplebasename + '-tumor'
		normalsample = samplebasename + '-norm'
		dpsicolname = samplebasename + '_dpsi'
		bigpsidf[dpsicolname] = bigpsidf[tumorsample] - bigpsidf[normalsample]

	bigpsidf.to_csv('{0}psis.txt'.format(os.path.basename(os.path.abspath(salmondir))), sep = '\t', index = False)


def getdpsis(salmondir):
	#Given a table of psis, calculate dpsis and LME-based p values
	pvaluedict = OrderedDict() #{genename : pvalue} ordered so it's easy to match it up with q values

	psifile = os.path.basename(os.path.abspath(salmondir)) + 'psis.txt'
	psidf = pd.read_table(psifile, sep = '\t', header = 0, index_col = False)

	pandas2ri.activate() #allow conversion between pandas dataframes and r dataframes

	#define R packages
	nlme = importr('nlme')
	base = importr('base')
	stats = importr('stats')
	qv = importr('qvalue')

	#define formulae
	fmla = Formula('value ~ 1 + cond1')
	rndm = Formula('~ 1 | samples')
	nullfmla = Formula('value ~ 1')
	nullrndm = Formula('~1 | samples')

	#Remove any gene that has a psi of NA in any sample
	psidf = psidf.dropna(axis=0)

	#Get name of normal sample names and tumor sample names
	salmondirs = [os.path.basename(d) for d in os.listdir(salmondir) if os.path.isdir(os.path.join(os.path.abspath(salmondir), d)) and d.endswith('_salmon')]
	samplenames = [s.split('_')[0] for s in salmondirs]
	normalnames = [sample for sample in samplenames if 'norm' in sample]
	tumornames = [sample for sample in samplenames if 'tumor' in sample]

	#Store relationships of conditions and the samples in that condition
	#It's important that this dictionary be ordered because we are going to be iterating through it
	samp_conds = OrderedDict({'cond1' : normalnames,
		'cond2' : tumornames})

	#Get a list of all samples
	samps = []
	for cond in samp_conds:
		samps += samp_conds[cond]

	#Iterate through rows, making a dictionary from every row, turning it into a dataframe, then calculating p value
	genecounter = 0
	for index, row, in psidf.iterrows():
		genecounter +=1
		if genecounter % 1000 == 0:
			print 'Gene {0}...'.format(genecounter)

		d = {}
		d['Gene'] = [row['Gene']] * len(samps)
		d['variable'] = samps
		
		values = [] #psi values
		for cond in samp_conds:
			for sample in samp_conds[cond]:
				value = row[sample]
				values.append(value)
		d['value'] = values
		
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

		#Get LME pvalue
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

		pvaluedict[row['Gene']] = pvalue

	#Get qvalues
	#Turn list of pvalues into R vector
	pvalues = pvaluedict.values()
	pvec = FloatVector(pvalues)
	#Get qvalues object
	qobj = qv.qvalue(p = pvec)
	#qvalues are index 2 of qvalue object
	qvalues = list(qobj[2])
	#format decimal
	qvalues = [float('{:.2e}'.format(qvalue)) for qvalue in qvalues]

	#Add pvalues and qvalues to df
	psidf = psidf.assign(pval = pvalues)
	psidf = psidf.assign(qval = qvalues)

	#Write df
	fn = os.path.abspath(psifile) + '.pval'
	psidf.to_csv(fn, sep = '\t', index = False)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'Gff of annotations.')
	parser.add_argument('--cancerdir', type = str, help = 'Directory containing all salmon directories for a particular cancer.')
	args = parser.parse_args()

	calculatepsi_multisample(args.gff, args.cancerdir)
	getdpsis(args.cancerdir)

	#calculatepsi_multisample(sys.argv[1], sys.argv[2])
	#getdpsis(sys.argv[2])
















