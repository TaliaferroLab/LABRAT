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
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import math



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
#Also, this function is actually getting you m values, not position factors, but it's what we want here

def getpositionfactors(gff, lengthfilter):
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
			txname = str(transcript.id).replace('transcript:', '').split('.')[0]
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

			#For last two exons
			TFs[txname] = [transcript.chrom, transcript.strand, exons[-2], exons[-1]]

			#For all exons
			#TFs[txname] = [transcript.chrom, transcript.strand]
			#for exon in exons:
				#TFs[txname].append(exon)

	
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
	'''

	#Get sequences of all exons
	with open('TFseqs.fasta', 'w') as outfh:
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
	'''

def fastqtrimmer(threeprimetrim, forreads, revreads):
	#Maybe you want to trim the reads from the 3' end before giving them to salmon.
	#Trim <threeprimetrim> nt from the 3' end of the read
	threeprimetrim = int(threeprimetrim)

	counter = 0
	foutfilename = 'tempf.fastq'
	routfilename = 'tempr.fastq'
	with gzip.open(forreads, 'rb') as forinfh, gzip.open(revreads, 'rb') as revinfh, open(foutfilename, 'w') as foroutfh, open(routfilename, 'w') as revoutfh:
		try:
			for title, seq, qual in FastqGeneralIterator(forinfh):
				counter +=1
				if counter % 10000000 == 0:
					print 'On read {0} of {1}.'.format(counter, forreads)
				foroutfh.write('@{0}\n{1}\n+\n{2}\n'.format(title, seq[:threeprimetrim * -1], qual[:threeprimetrim* -1]))
		except ValueError:
			pass

		try:
			counter = 0
			for title, seq, qual in FastqGeneralIterator(revinfh):
				counter +=1
				if counter % 10000000 == 0:
					print 'On read {0} of {1}.'.format(counter, revreads)
				revoutfh.write('@{0}\n{1}\n+\n{2}\n'.format(title, seq[:threeprimetrim * -1], qual[:threeprimetrim* -1]))

		except ValueError:
			pass

	print 'Done trimming {0} and {1}.'.format(forreads, revreads)

def runSalmon(transcriptfasta, threads, reads1, reads2, samplename):
	#transcriptfasta should probably be 'TFseqs.fasta' produced by makeTFfasta
	#index transcripts if index does not already exist
	if os.path.exists('/vol3/home/taliaferro/data/T1DNeuron/LABRAT/hg38.tfseqs.idx') == False:
		command = ['salmon', 'index', '-t', transcriptfasta, '-i', '/vol3/home/taliaferro/data/T1DNeuron/LABRAT/hg38.tfseqs.idx', '--type', 'quasi', '-k', '31']
		print 'Indexing transcripts...'
		subprocess.call(command)
		print 'Done indexing!'

	#paired end
	if reads2:
		print 'Running salmon for {0}...'.format(samplename)
		command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--seqBias', '--gcBias', '-1', reads1, '-2', reads2, '-o', samplename, '--index', '/vol3/home/taliaferro/data/T1DNeuron/LABRAT/hg38.tfseqs.idx']
		#command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--seqBias', '--gcBias', '-1', reads1, '-2', reads2, '-o', samplename, '--index', '/vol3/home/taliaferro/Annotations/dm6/dm6.entiretranscript.idx/']

	#Single end
	elif not reads2:
		fldMean = '250' #fragment length distribution mean
		fldSD = '20' #fragment length distribution standard deviation	
		command = ['salmon', 'quant', '--libType', 'A', '-p', threads, '--fldMean', fldMean, '--fldSD', fldSD, '--seqBias', '-r', reads1, '-o', samplename, '--index', '../transcripts.mm10.idx']

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
			transcriptid = str(line[0])
			tpm = float(line[3])
			txtpms[transcriptid] = tpm

	#Put the transcript tpms for every gene together
	for gene in positionfactors:
		genetpms[gene] = []
		posfactorgenetpms[gene] = []
		for transcript in positionfactors[gene]:
			txtpm = txtpms[transcript]
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

	#cond1 = soma, cond2 = neurite

	#samp_conds = OrderedDict({'cond1': ['GFPctrlSomaA', 'GFPctrlSomaB', 'GFPctrlSomaC', 'GFPhypoxiaSomaA', 'GFPhypoxiaSomaB', 'GFPhypoxiaSomaC',
		#'T1DctrlSomaA', 'T1DctrlSomaB', 'T1DctrlSomaC', 'T1DhypoxiaSomaA', 'T1DhypoxiaSomaB', 'T1DhypoxiaSomaC'], 
		#'cond2' : ['GFPctrlNeuriteA', 'GFPctrlNeuriteB', 'GFPctrlNeuriteC', 'GFPhypoxiaNeuriteA', 'GFPhypoxiaNeuriteB', 'GFPhypoxiaNeuriteC',
		#'T1DctrlNeuriteA', 'T1DctrlNeuriteB', 'T1DctrlNeuriteC', 'T1DhypoxiaNeuriteA', 'T1DhypoxiaNeuriteB', 'T1DhypoxiaNeuriteC']})

	#samp_conds = OrderedDict({'cond1': ['GFPctrlSomaA', 'GFPctrlSomaB', 'GFPctrlSomaC', 'GFPhypoxiaSomaA', 'GFPhypoxiaSomaB', 'GFPhypoxiaSomaC'], 
		#'cond2' : ['GFPctrlNeuriteA', 'GFPctrlNeuriteB', 'GFPctrlNeuriteC', 'GFPhypoxiaNeuriteA', 'GFPhypoxiaNeuriteB', 'GFPhypoxiaNeuriteC']})

	samp_conds = OrderedDict({'cond1': ['GFPctrlSomaA', 'GFPctrlSomaB', 'GFPctrlSomaC'], 
		'cond2' : ['GFPctrlNeuriteA', 'GFPctrlNeuriteB', 'GFPctrlNeuriteC']})

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
	parser.add_argument('--mode', type = str, choices = ['makeTFfasta', 'runSalmon', 'calculatepsi', 'LME'])
	parser.add_argument('--gff', type = str, help = 'GFF of transcript annotation. Needed for makeTFfasta and calculatepsi.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in fasta format. Needed for makeTFfasta.')
	parser.add_argument('--reads1', type = str, help = 'Comma separated list of forward read fastq files. Needed for runSalmon.')
	parser.add_argument('--reads2', type = str, help = 'Comma separated list of reverse read fastq files. Needed for runSalmon.')
	parser.add_argument('--samplename', type = str, help = 'Comma separated list of samplenames.  Needed for runSalmon.')
	parser.add_argument('--threeprimetrim', type = int, help = 'Number of bases to trim from 3\' end of reads.  Optional for runSalmon.')
	parser.add_argument('--threads', type = str, help = 'Number of threads to use.  Needed for runSalmon.')
	parser.add_argument('--salmondir', type = str, help = 'Salmon output directory. Needed for calculatepsi.')
	parser.add_argument('--psifile', type = str, help = 'Psi value table. Needed for LME.')
	args = parser.parse_args()

	if args.mode == 'makeTFfasta':
		makeTFfasta(args.gff, args.genomefasta)

	elif args.mode == 'runSalmon':
		forreads = args.reads1.split(',')
		revreads = args.reads2.split(',')
		samplenames = args.samplename.split(',')
		if len(forreads) != len(revreads) != len(samplenames):
			print 'Need the same number of forward reads, reverse reads, and sample names!'
			sys.exit()

		for i in range(len(forreads)):
			freads = forreads[i]
			rreads = revreads[i]
			samplename = samplenames[i]
			if args.threeprimetrim:
				print 'Trimming sample {0}...'.format(samplename)
				fastqtrimmer(args.threeprimetrim, freads, rreads)
			runSalmon('TFseqs.hg38.fasta', args.threads, 'tempf.fastq', 'tempr.fastq', samplename)

		if args.threeprimetrim:
			os.remove('tempf.fastq')
			os.remove('tempr.fastq')

	elif args.mode == 'calculatepsi':
		print 'Calculating position factors for every transcript...'
		positionfactors = getpositionfactors(args.gff, 20)
		print 'Done with position factors!'
		salmondirs = [os.path.join(os.path.abspath(args.salmondir), d) for d in os.listdir(args.salmondir) if os.path.isdir(os.path.join(os.path.abspath(args.salmondir), d))]
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
		bigpsidf.to_csv('ENCODEpsis.3end.txt', sep = '\t', index = False)

	elif args.mode == 'LME':
		getdpsis(args.psifile)


		