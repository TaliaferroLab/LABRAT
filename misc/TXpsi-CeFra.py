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


#For every gene, the goal is to define what the overall usage of proximal/distal polyA sites is.  This is done by defining 
#a "psi" value for every gene.  For each transcript in a gene, the "terminal fragment" (TF) is the last two
#exons of the transcript concatenated together.  For each transcript, a position factor (PF) is calculated as (m / n + 1; will range 
#between 0 and 1).
#The tpm for its TF is multiplied by its PF.  This is then done for all TFs and their scaled values are summed and divided
#by the unscaled sum of tpms for all TFs to give psi, which will also range between 0 (exculsive usage of the proximal polyA site)
#and 1 (exclusive usage of the distal polyA site).  

#Does a gene pass filters?
def genefilters(gene, db):
	return True #not using filter at all
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
	#Not used in Drosophila
	'''
	if 'protein_coding' in transcript.attributes['transcript_type']:
		proteincoding = True
	'''
	

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



#Make a fasta file containing the "terminal fragments" of all transcripts
def makeTFfasta(gff, genomefasta):
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

			TFs[txname] = [transcript.chrom, transcript.strand, exons[-2], exons[-1]]

	#Get sequences of TFs
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

def runSalmon(transcriptfasta, threads, reads1, reads2, samplename):
	#transcriptfasta should probably be 'TFseqs.fasta' produced by makeTFfasta
	#index transcripts if index does not already exist
	
	if os.path.exists('/vol3/home/taliaferro/Annotations/dm6/dm6.entiretranscript.idx/') == False:
		command = ['salmon', 'index', '-t', transcriptfasta, '-i', 'transcripts.idx', '--type', 'quasi', '-k', '31']
		print 'Indexing transcripts...'
		subprocess.call(command)
		print 'Done indexing!'

	#paired end
	if reads2:
		print 'Running salmon for {0}...'.format(samplename)
		command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--seqBias', '--gcBias', '-1', reads1, '-2', reads2, '-o', samplename, '--index', '/vol3/home/taliaferro/Annotations/dm6/dm6.entiretranscript.idx/']
		#command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--seqBias', '--gcBias', '-1', reads1, '-2', reads2, '-o', samplename, '--index', 'transcripts.dm6.idx']
	
	#Single end
	elif not reads2:
		fldMean = '250' #fragment length distribution mean
		fldSD = '20' #fragment length distribution standard deviation	
		command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--fldMean', fldMean, '--fldSD', fldSD, '--seqBias', '-r', reads1, '-o', samplename, '--index', 'transcripts.idx']

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
		if totaltpm > 5:
			psi = scaledtpm / float(totaltpm)
			psis[gene] = psi
		else:
			psis[gene] = 'NA'

	return psis

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



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--mode', type = str, choices = ['makeTFfasta', 'runSalmon', 'calculatepsi', 'LME'])
	parser.add_argument('--gff', type = str, help = 'GFF of transcript annotation. Needed for makeTFfasta and calculatepsi.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format. Needed for makeTFfasta.')
	parser.add_argument('--reads1', type = str, help = 'Comma separated list of forward read fastq files. Needed for runSalmon.')
	parser.add_argument('--reads2', type = str, help = 'Comma separated list of reverse read fastq files. Needed for runSalmon.')
	parser.add_argument('--samplename', type = str, help = 'Comma separated list of samplenames.  Needed for runSalmon.')
	parser.add_argument('--threads', type = str, help = 'Number of threads to use.  Needed for runSalmon.')
	parser.add_argument('--salmondir', type = str, help = 'Salmon output directory. Needed for calculatepsi.')
	parser.add_argument('--psifile', type = str, help = 'Psi value table. Needed for LME.')
	args = parser.parse_args()

	if args.mode == 'makeTFfasta':
		makeTFfasta(args.gff, args.genomefasta)


	elif args.mode == 'runSalmon':
		readdir = '/vol3/home/taliaferro/data/CeFra/RawReads'
		celllines = ['K562']
		fracs = ['nucleus', 'total']
		methods = ['polyA', 'ribodep']
		reps = ['Rep1', 'Rep2']
		
		forwardreads = [] # list of all forward read files
		outputnames = []
		for cellline in celllines:
			for frac in fracs:
				for method in methods:
					for rep in reps:
						filename = [cellline, frac, method, rep, '1', 'fastq', 'gz']
						forwardreads.append(('.').join(filename))
						outputnames.append(('_').join([cellline, frac, method, rep]))

		reversereads = [] # list of all forward read files
		for cellline in celllines:
			for frac in fracs:
				for method in methods:
					for rep in reps:
						filename = [cellline, frac, method, rep, '2', 'fastq', 'gz']
						reversereads.append(('.').join(filename))

		forwardreads = [os.path.join(readdir, forread) for forread in forwardreads]
		reversereads = [os.path.join(readdir, revread) for revread in reversereads]

		for i in range(len(outputnames)):
			forreads = forwardreads[i]
			revreads = reversereads[i]
			outputname = outputnames[i]
			runSalmon('TFseqs.fasta', args.threads, forreads, revreads, outputname)



	elif args.mode == 'calculatepsi':
		print 'Calculating position factors for every transcript...'
		positionfactors = getpositionfactors(args.gff, 25)
		genecount = len(positionfactors)
		txcount = 0
		for gene in positionfactors:
			for transcript in positionfactors[gene]:
				txcount +=1
		print genecount, txcount
		print 'Done with position factors!'
		salmondirs = [os.path.join(os.path.abspath(args.salmondir), d) for d in os.listdir(args.salmondir) if os.path.isdir(os.path.join(os.path.abspath(args.salmondir), d))]
		psidfs = []
		for sd in salmondirs:
			#samplename = os.path.basename(sd)
			samplename = ('_').join(os.path.basename(sd).split('_')[1:])
			psis = calculatepsi(positionfactors, sd)
			psidf = pd.DataFrame.from_dict(psis, orient = 'index')
			psidf.reset_index(level = 0, inplace = True)
			psidf.columns = ['Gene', samplename]
			psidfs.append(psidf)
		bigpsidf = reduce(lambda x, y: pd.merge(x, y, on = 'Gene'), psidfs)
		bigpsidf.to_csv('DrosophilaPSIvalues.txt', sep = '\t', index = False)

	elif args.mode == 'LME':
		getdpsis(args.psifile)


		