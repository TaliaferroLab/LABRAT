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
#from rpy2.robjects.packages import importr
#from rpy2.robjects import pandas2ri, Formula, FloatVector
#from rpy2.rinterface import RRuntimeError
import math
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy.stats.distributions import chi2
import warnings



#For every gene, the goal is to define what the overall usage of proximal/distal polyA sites is.  This is done by defining 
#a "psi" value for every gene.  For each transcript in a gene, the "terminal fragment" (TF) is the last two
#exons of the transcript concatenated together.  For each transcript, a position factor (PF) is calculated as (m / n + 1; will range 
#between 0 and 1).
#The tpm for its TF is multiplied by its PF.  This is then done for all TFs and their scaled values are summed and divided
#by the unscaled sum of tpms for all TFs to give psi, which will also range between 0 (exculsive usage of the proximal polyA site)
#and 1 (exclusive usage of the distal polyA site).  

#Does a gene pass filters?
def genefilters(gene, db):
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

	#Output file of the number of posfactors for each gene
	with open('numberofposfactors.txt', 'w') as outfh:
		outfh.write(('\t').join(['Gene', 'numberofposfactors']) + '\n')
		for gene in posfactors:
			pfs = []
			for tx in posfactors[gene]:
				pfs.append(posfactors[gene][tx])
			pfs = list(set(pfs))
			outfh.write(('\t').join([gene, str(len(pfs))]) + '\n')


	return posfactors


#Make a fasta file containing the "terminal fragments" of all transcripts
def makeTFfasta(gff, genomefasta, lasttwoexons):
	TFs = {} #{transcriptid: [chrm, strand, [next to last exon1start, next to last exon1stop], [last exon2start, last exon2stop]]}
	
	#Make gff database
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	print 'Indexing genome sequence...'
	if os.path.basename(genomefasta).endswith('.gz'):
		seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	else:
		seq_dict = SeqIO.to_dict(SeqIO.parse(genomefasta, 'fasta'))
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	#Get last two exons
	genecounter = 0
	for gene in genes:
		genecounter +=1
		if genecounter % 5000 == 0:
			print 'Gene {0}...'.format(genecounter)
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
				outfh.write('>' + TF + '\n' + str(seq) + '\n')

def runSalmon(threads, reads1, reads2, samplename):
	#paired end
	if reads2:
		print 'Running salmon for {0}...'.format(samplename)
		command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--seqBias', '--gcBias', '-1', reads1, '-2', reads2, '-o', samplename, '--index', 'txfasta.idx']

	#Single end
	elif not reads2:
		fldMean = '250' #fragment length distribution mean
		fldSD = '20' #fragment length distribution standard deviation	
		command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--fldMean', fldMean, '--fldSD', fldSD, '--seqBias', '-r', reads1, '-o', samplename, '--index', 'txfasta.idx']

	subprocess.call(command)

def calculatepsi(positionfactors, salmondir):
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
			transcriptid = str(line[0]).split('.')[0]
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
		if totaltpm >= 5:
			psi = scaledtpm / float(totaltpm)
			psis[gene] = float(format(psi, '.3f'))
		else:
			#psis[gene] = 'NA'
			psis[gene] = np.nan

	return psis

def getdpsis_python(psifile, samp_conds_file):
	#Given a table of psis, calculate dpsis and LME-based p values
	deltapsidict = OrderedDict() #{genename : pvalue} ordered so it's easy to match it up with df
	pvaluedict = OrderedDict() #{genename : pvalue} ordered so it's easy to match it up with q values

	psidf = pd.read_csv(psifile, sep = '\t', header = 0, index_col = False)

	cond1samps = []
	cond2samps = []
	with open(samp_conds_file, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			cond1samps.append(line[0])
			cond2samps.append(line[1])

	print 'Condition 1 samples: ' + (', ').join(cond1samps)
	print 'Condition 2 samples: ' + (', ').join(cond2samps)

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
			print 'Calculating pvalue for gene {0}...'.format(genecounter)
		
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

		#delta psi is difference between mean psi of two conditions (cond2 - cond1)
		cond1meanpsi = float(format(np.mean(rowdf.query('cond1 == 1').value.dropna()), '.3f'))
		cond2meanpsi = float(format(np.mean(rowdf.query('cond2 == 1').value.dropna()), '.3f'))
		deltapsi = cond2meanpsi - cond1meanpsi
		deltapsi = float(format(deltapsi, '.3f'))
		deltapsidict[row['Gene']] = deltapsi

		#Get LME pvalue
		#Lots of warnings about convergence, etc. Suppress them.
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
				if math.isnan(p):
					p = 1
			#These exceptions are needed to catch cases where either all psi values are nan (valueerror) or all psi values for one condition are nan (linalgerror)
			except (ValueError, np.linalg.LinAlgError):
				p = 1

			pvaluedict[row['Gene']] = float('{:.2e}'.format(p))

	#Correct pvalues using BH method
	pvalues = list(pvaluedict.values())
	fdrs = list(multipletests(pvalues, method = 'fdr_bh')[1])
	fdrs = [float('{:.2e}'.format(fdr)) for fdr in fdrs]

	#Add deltapsis, pvalues, and FDRs to df
	deltapsis = list(deltapsidict.values())
	psidf = psidf.assign(deltapsi = deltapsis)
	psidf = psidf.assign(pval = pvalues)
	psidf = psidf.assign(FDR = fdrs)

	#Write df
	fn = os.path.abspath(psifile) + '.python.pval'
	psidf.to_csv(fn, sep = '\t', index = False)



def getdpsis(psifile):
	#Given a table of psis, calculate dpsis and LME-based p values
	pvaluedict = OrderedDict() #{genename : pvalue} ordered so it's easy to match it up with q values

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

	#Store relationships of conditions and the samples in that condition
	#It's important that this dictionary be ordered because we are going to be iterating through it
	samp_conds = OrderedDict({'cond1' : ['ER1', 'ER2', 'ER3'],
		'cond2' : ['Cytosol1', 'Cytosol2', 'Cytosol3']})

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
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

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
					txexons += range(exon.start, exon.end + 1)
			elif tx.strand == '-':
				for exon in db.children(tx, featuretype = 'exon', order_by = 'start', reverse = True):
					txexons += range(exon.start, exon.end + 1)

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
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

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
	parser.add_argument('--mode', type = str, choices = ['makeTFfasta', 'runSalmon', 'calculatepsi', 'LME'])
	parser.add_argument('--gff', type = str, help = 'GFF of transcript annotation. Needed for makeTFfasta and calculatepsi.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format. Needed for makeTFfasta.')
	parser.add_argument('--lasttwoexons', action = 'store_true', help = 'Used for makeTFfasta. Do you want to only use the last two exons?')
	parser.add_argument('--txfasta', type = str, help = 'Fasta file of sequences to quantitate with salmon. Often the output of makeTFfasta mode.')
	parser.add_argument('--reads1', type = str, help = 'Comma separated list of forward read fastq files. Needed for runSalmon.')
	parser.add_argument('--reads2', type = str, help = 'Comma separated list of reverse read fastq files. Needed for runSalmon.')
	parser.add_argument('--samplename', type = str, help = 'Comma separated list of samplenames.  Needed for runSalmon.')
	parser.add_argument('--threads', type = str, help = 'Number of threads to use.  Needed for runSalmon.')
	parser.add_argument('--salmondir', type = str, help = 'Salmon output directory. Needed for calculatepsi.')
	parser.add_argument('--psifile', type = str, help = 'Psi value table. Needed for LME.')
	parser.add_argument('--sampconds', type = str, help = 'File containing sample names split by condition. Two column, tab delimited text file. Condition 1 samples in first column, condition2 samples in second column.')
	args = parser.parse_args()

	if args.mode == 'makeTFfasta':
		makeTFfasta(args.gff, args.genomefasta, args.lasttwoexons)

	elif args.mode == 'runSalmon':
		forreads = args.reads1.split(',')
		revreads = args.reads2.split(',')
		samplenames = args.samplename.split(',')
		if len(forreads) != len(revreads) != len(samplenames):
			print 'Need the same number of forward reads, reverse reads, and sample names!'
			sys.exit()

		command = ['salmon', 'index', '-t', args.txfasta, '-i', 'txfasta.idx', '--type', 'quasi', '-k', '31']
		print 'Indexing transcripts...'
		subprocess.call(command)
		print 'Done indexing!'

		for i in range(len(forreads)):
			freads = forreads[i]
			rreads = revreads[i]
			samplename = samplenames[i]
			runSalmon(args.threads, freads, rreads, samplename)

	elif args.mode == 'calculatepsi':
		print 'Calculating position factors for every transcript...'
		positionfactors = getpositionfactors(args.gff, 25)
		print 'Done with position factors!'
		salmondirs = [os.path.join(os.path.abspath(args.salmondir), d) for d in os.listdir(args.salmondir) if os.path.isdir(os.path.join(os.path.abspath(args.salmondir), d)) and d != 'txfasta.idx']
		psidfs = []
		for sd in salmondirs:
			samplename = os.path.basename(sd)
			
			#For ENCODE
			#samplename = os.path.basename(sd).split('_')
			#samplename = ('_rep').join([samplename[0], samplename[1]])
			
			psis = calculatepsi(positionfactors, sd)
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
		finalpsidf.to_csv('LABRATpsis.3end.txt', sep = '\t', index = False)

	elif args.mode == 'LME':
		getdpsis_python(args.psifile, args.sampconds)


		