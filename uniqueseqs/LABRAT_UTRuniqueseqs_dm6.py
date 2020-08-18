#This script is designed to get the "unique" UTR sequences that define each position factor for each gene.
#This is specifically made for an Ensembl Drosophila annotation (dm6, BDGP6). The one we have in Ensembl 88.


import gffutils
import os
from collections import defaultdict
import sys
from itertools import groupby, chain
from operator import itemgetter
import pickle
import subprocess
import pybedtools
from Bio import SeqIO
import gzip
import argparse
import random

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

	if 'protein_coding' in gene.attributes['biotype']:
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
	if 'protein_coding' in transcript.attributes['biotype']:
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
		outfh.write(('\t').join(['Gene', 'numberofposfactors']) + '\n')
		for gene in posfactors: #{ENSMUSG : {ENSMUST : positionfactor}}
			pfs = []
			for tx in posfactors[gene]:
				pfs.append(posfactors[gene][tx])
			pfs = list(set(pfs))
			txids = {} #{positionfactor : [list of transcriptIDs]}
			for pf in sorted(pfs):
				txids[pf] = []
				for tx in posfactors[gene]:
					if posfactors[gene][tx] == pf:
						txids[float(pf)].append(tx)

			alltxids = []
			for pf in sorted(list(txids.keys())):
				alltxids.append((',').join(txids[pf]))
			outfh.write(('\t').join([gene, str(len(pfs)), ('\t').join(alltxids)]) + '\n')


	return posfactors



#Given a transcript, retrieve it's UTR coords
def getUTRcoords(txid, db):
	exoncoords = [] #[[exon1start, exon1stop], [exon2start, exon2stop]]
	CDScoords = []
	UTRcoords = [] #[UTRstart, UTRstop]

	transcript = db['transcript:' + txid]
	noUTRcount = 0

	for exon in db.children(transcript, featuretype = 'exon', order_by = 'start'):
		exoncoords.append([exon.start, exon.end])
	for CDSexon in db.children(transcript, featuretype = 'CDS', order_by = 'start'):
		CDScoords.append([CDSexon.start, CDSexon.end])

	#3' UTR start is directly after CDS end
	if transcript.strand == '+':
		CDSend = max(CDScoords, key = itemgetter(1))[1]
		#If the transcript ends right were the CDS ends, define the UTR as the 25 nt downstream of CDSend
		if CDSend == transcript.end:
			UTRcoords = [transcript.end + 1, transcript.end + 26]
			print('NO UTR')
			return UTRcoords
		else:
			UTR3start = CDSend + 1
			UTRcoords = [UTR3start, transcript.end]
	elif transcript.strand == '-':
		CDSend = min(CDScoords, key = itemgetter(0))[0]
		if CDSend == transcript.start: # If the transcript ends right where the CDS ends, define the UTR as the 25 nt upstream of CDSstart
			UTRcoords = [transcript.start - 25, transcript.start]
			print('NO UTR')
			return UTRcoords
		else:
			UTR3start = CDSend - 1
			UTRcoords = [transcript.start, UTR3start]

	###Check to see if the UTR is fully contained within the coordinates of one exon
	singleexonUTR = False
	for exoncoord in exoncoords:
		exonstart = exoncoord[0]
		exonend = exoncoord[1]
		if exonstart <= UTRcoords[0] and exonend >= UTRcoords[1]:
			singleexonUTR = True
			return UTRcoords

	if singleexonUTR == False:
		#Get all positions that are both exonic and in the 3' UTR
		overlappingbp = [] #sorted exonic positions in UTR
		UTR3range = range(UTRcoords[0], UTRcoords[1] + 1)
		for exoncoord in exoncoords:
			exonrange = range(exoncoord[0], exoncoord[1] + 1)
			overlap = set(UTR3range).intersection(exonrange)
			for nt in sorted(list(overlap)):
				overlappingbp.append(nt)

		#Now get breaks in consecutive exonic positions
		#http://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
		UTRexoncoords = [] #[[exon1start, exon1stop], [exon2start, exon2stop]]
		for k, g in groupby(enumerate(overlappingbp), lambda ix: ix[0] - ix[1]):
			exonbp = list(map(itemgetter(1), g))
			if len(exonbp) > 1:
				UTRexoncoords.append([exonbp[0], exonbp[-1]])

	return UTRexoncoords


#Get a UTR for each transcript that has a posfactor
def getallUTRcoords(gff, posfactors):
	#posfactors = {} #{ENSMUSG : {ENSMUST : positionfactor}}
	UTRcoords = {} #{ENSMUSG : {posfactor : {ENSMUST1 : [coords], ENSMUST2 : [coords]}}}
	print(len(posfactors))
	#Make gff database
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	print('Getting UTR coords for all transcripts...')

	genecounter = 0
	for gene in posfactors:
		genecounter +=1
		print(genecounter, gene, len(posfactors[gene]))
		#Only genes that have more than one polyA site
		if len(posfactors[gene]) == 1:
			continue
		if gene not in UTRcoords:
			UTRcoords[gene] = {}
		for transcript in posfactors[gene]:
			posfactor = posfactors[gene][transcript]
			txUTRcoords = getUTRcoords(transcript, db)
			#Sometimes a transcript does not have a UTR (the end of the CDS == the end of the tx) meaning txUTRcoords is None
			if txUTRcoords == None:
				print('WHAT')
				continue
			if posfactor not in UTRcoords[gene]:
				UTRcoords[gene][posfactor] = {}
			UTRcoords[gene][posfactor][transcript] = txUTRcoords

	return UTRcoords



def mergeUTRs(txids, db, listofUTRcoords):
	#Use bedtools merge to merge UTRs

	#Need to make a bed
	#Get chrms and strand from txids
	chrms = []
	strands = []
	for tx in txids:
		tx = 'transcript:' + tx
		chrm = db[tx].chrom
		strand = db[tx].strand
		chrms.append(chrm)
		strands.append(strand)

	#They better be all the same chromosome and strand
	chrms = list(set(chrms))
	strands = list(set(strands))
	if len(chrms) > 1 or len(strands) > 1:
		print('ERROR: Merge items from different chromosomes or strands!')
		sys.exit()

	#Need each feature (exon) to be a single line
	exons = []
	print('listofUTRcoords: ', listofUTRcoords)
	for utr in listofUTRcoords:
		print('UTR: ', utr)
		#If this is a single exon UTR, this will NOT be a nested list
		if any(isinstance(i, list) for i in utr) == False:
			exons.append([utr[0], utr[1]])
		#If this is not a single exon UTR, we need each exon individually
		elif any(isinstance(i, list) for i in utr) == True:
			for exon in utr:
				exons.append([exon[0], exon[1]])


	#Make bed file
	with open('temp.bed', 'w') as outfh:
		for exon in exons:
			outfh.write(('\t').join([chrms[0], str(exon[0]), str(exon[1]), '.', '1000', strands[0]]) + '\n')

	#Now sort bed by chrm and start position
	command = ['sort', '-k1,1', '-k2,2n', 'temp.bed']
	with open('temp.sorted.bed', 'w') as outfh:
		subprocess.call(command, stdout = outfh)

	#Merge the UTR exons
	bed = pybedtools.BedTool('temp.sorted.bed')
	mergedbed = bed.merge()
	bedlines = []
	for exon in mergedbed:
		bedlines.append([int(exon[1]), int(exon[2])])

	mergedexons = []
	if len(bedlines) == 1:
		mergedexons = [bedlines[0][0], bedlines[0][1]]
	elif len(bedlines) > 1:
		for exon in bedlines:
			mergedexons.append([int(exon[0]), int(exon[1])])

	os.remove('temp.bed')
	os.remove('temp.sorted.bed')
	return mergedexons



def mergeUTRcoords(gff, UTRcoords):
	#Go through UTRcoords (output of getallUTRcoords).
	#If there are any position factors that have multiple transcripts associated with them, use mergeUTRs
	#to merge their UTRs. Then replace their entry in UTRcoords with this merged UTR.
	#After this function, all position factors will have a single coordinate list associated with them.
	#UTRcoords = {} #{ENSMUSG : {posfactor : {ENSMUST1 : [coords], ENSMUST2 : [coords]}}}

	#Make gff database
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	print('Merging UTRs that have the same polyA site...')

	cleanedUTRcoords = {} #{ENSMUSG : {posfactor1 : [coords], posfactor2 : [coords]}}

	for gene in UTRcoords:
		#It's possible to have a gene that only has one posfactor that actually has utr coordinates associated with it
		#This is because sometimes a tx would not have a UTR, but still it ends up in posfactors dict
		#If a gene does not have at least two different posfactors with utr coordinates associated with them, skip this gene.
		if len(UTRcoords[gene]) <= 1:
			continue
		cleanedUTRcoords[gene] = {}
		for posfactor in UTRcoords[gene]:
			#If there is only one UTR associated with this posfactor, we don't need to merge
			if len(UTRcoords[gene][posfactor]) == 1:
				txid = list(UTRcoords[gene][posfactor].keys())[0]
				cleanedUTRcoords[gene][posfactor] = UTRcoords[gene][posfactor][txid]
			elif len(UTRcoords[gene][posfactor]) > 1:
				txids = []
				utrs = []
				print(gene)
				for transcript in UTRcoords[gene][posfactor]:
					txids.append(transcript)
					utrs.append(UTRcoords[gene][posfactor][transcript])

				mergedexons = mergeUTRs(txids, db, utrs)
				cleanedUTRcoords[gene][posfactor] = mergedexons

	return cleanedUTRcoords



def subtractUTRs(dictofUTRcoords, geneid, db):
	#Use bedtools subtract to get unique regions in each UTR
	uniqueutrs = {} #{utr number : [coords]}

	#Get chromosome and strand from geneid
	geneid = 'gene:' + geneid
	chrm = db[geneid].chrom
	strand = db[geneid].strand

	#Get number of UTRs for this gene
	numberofutrs = len(dictofUTRcoords)

	#In this dict, the keys are the order of the UTRs.
	#Write each UTR to a bed.
	print(dictofUTRcoords)
	for i in range(numberofutrs):
		pf = float(i / (numberofutrs - 1))
		#Write UTR to bed
		utr = dictofUTRcoords[pf]
		print(utr)

		#Need each feature (exon) to be a single line
		utrexons = []
		#If this is a single exon UTR, this will NOT be a nested list
		if any(isinstance(i, list) for i in utr) == False:
			utrexons.append([utr[0], utr[1]])
		#If this is not a single exon UTR, we need each exon individually
		elif any(isinstance(i, list) for i in utr) == True:
			for exon in utr:
				utrexons.append([exon[0], exon[1]])

		#Make bed file
		with open('temp{0}.bed'.format(i), 'w') as outfh:
			for exon in utrexons:
				#If this is a feature of length 1 (start = stop in gff), add one to the stop so that bedtools handles it correctly
				if exon[0] == exon[1]:
					outfh.write(('\t').join([chrm, str(exon[0]), str(exon[1] + 1), '.', '1000', strand]) + '\n')
				else:
					outfh.write(('\t').join([chrm, str(exon[0]), str(exon[1]), '.', '1000', strand]) + '\n')

		#Now sort bed by chrm and start position
		command = ['sort', '-k1,1', '-k2,2n', 'temp{0}.bed'.format(i)]
		with open('temp{0}.sorted.bed'.format(i), 'w') as outfh:
			subprocess.call(command, stdout = outfh)

	#Now subtract UTRs
	#If there are four: 1 = 1; 2 = 2 - 1; 3 = 3 - 2 - 1; 4 = 4 - 3 - 2 - 1
	for i in range(numberofutrs):
		exons = []
		#For the first UTR, there is no subtraction
		if i == 0:
			bed = pybedtools.BedTool('temp{0}.sorted.bed'.format(i))
			for exon in bed:
				exons.append([int(exon[1]), int(exon[2])])
			uniqueutrs[i] = exons

		else:
			#This is the bed we start with
			bed = pybedtools.BedTool('temp{0}.sorted.bed'.format(i))

			#Get all previous utr numbers
			previousutrnumbers = []
			for j in range(i):
				previousutrnumbers.append(j)
				#They need to be reversed because the subtraction goes like 4 = 4-3-2-1
				previousutrnumbers.reverse()

			for j in previousutrnumbers:
				subtractionbed = pybedtools.BedTool('temp{0}.sorted.bed'.format(j))
				newbed = bed.subtract(subtractionbed, s = True)
				bed = newbed

			exons = []
			for exon in bed:
				strand == exon[5]
				#Here we have to fix an off by one error that occurs because the subtraction is done on a
				#bed but we will be writing to gff
				if strand == '+':
					exons.append([int(exon[1]) + 1, int(exon[2])])
				elif strand == '-':
					exons.append([int(exon[1]), int(exon[2]) - 1])

			uniqueutrs[i] = exons

	for i in range(numberofutrs):
		os.remove('temp{0}.bed'.format(i))
		os.remove('temp{0}.sorted.bed'.format(i))

	print(uniqueutrs)

	#The start of UTR B should not be allowed to be before the end of UTR A. Fix it so this is so.
	uniqueutrs_filt = {} #{utrnumber : [coords]}
	previousend = 0
	for i in range(numberofutrs):
		#For the first UTR, just take the start coord
		coords = []
		for exon in uniqueutrs[i]:
			coords.append(exon[0])
			coords.append(exon[1])
		
		#For some reason, there can be some empty UTRs here.  Rare, but it happens.
		if not coords:
			continue

		if strand == '+':
			utrend = max(coords)
		elif strand == '-':
			utrend = min(coords)
		if i == 0:
			previousend = utrend
			uniqueutrs_filt[i] = uniqueutrs[i]
			continue
		elif i > 0:
			utrcoords = []
			if strand == '+':
				for exon in uniqueutrs[i]:
					#If both start and stop are greater than the previous end, we are good
					if exon[0] >= previousend and exon[1] >= previousend:
						utrcoords.append([exon[0], exon[1]])
					#If start is less than previousend, the new exon start where the last one ended
					elif exon[0] < previousend and exon[1] >= previousend:
						utrcoords.append([previousend + 1, exon[1]])
					#If both start and stop are less than the previous end, skip this exon
					elif exon[0] < previousend and exon[1] < previousend:
						continue

			elif strand == '-':
				for exon in uniqueutrs[i]:
					#If both start and stop are less than the previous end, we are good
					if exon[0] <= previousend and exon[1] <= previousend:
						utrcoords.append([exon[0], exon[1]])
					#If start is greater than previousend, the new exon start where the last one ended
					elif exon[0] <= previousend and exon[1] > previousend:
						utrcoords.append([exon[0], previousend - 1])
					#If both start and stop are greater than the previous end, skip this exon
					elif exon[0] > previousend and exon[1] > previousend:
						continue

			previousend = utrend
			uniqueutrs_filt[i] = utrcoords

	return uniqueutrs_filt

def getuniquecoords(gff, cleanedUTRcoords):
	#Make gff database
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	print('Subtracting UTRs to get unique sequences for {0} genes...'.format(len(cleanedUTRcoords)))

	uniqueutrcoords = {} #{ENSMUSG : {posfactor1 : [coords], posfactor2 : [coords]}}

	genecounter = 0
	for gene in cleanedUTRcoords:
		genecounter +=1
		if genecounter % 1000 == 0:
			print('Subtracting UTRs for gene {0} of {1}...'.format(genecounter, len(cleanedUTRcoords)))
		utrs = cleanedUTRcoords[gene]
		uniqueutrs = subtractUTRs(utrs, gene, db)
		uniqueutrcoords[gene] = uniqueutrs

	return uniqueutrcoords


def makegff(gff, uniqueutrcoords):
	#Given a dictionary of unique coords for each gene, make a gff
	#Make gff database
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	print('Making gff of unique UTR sequences for {0} genes...'.format(len(uniqueutrcoords)))

	with open('uniqueutrcoords.gff', 'w') as outfh:
		#Write line indicating what this gff was made from
		outfh.write('#gff_annotation = {0}'.format(os.path.abspath(gff)) + '\n')
		for gene in uniqueutrcoords:
			chrm = db['gene:' + gene].chrom
			strand = db['gene:' + gene].strand

			#Get min and max utr coord across all unique sequences
			coords = []
			for uniqueseq in uniqueutrcoords[gene]:
				for exon in uniqueutrcoords[gene][uniqueseq]:
					coords.append(exon[0])
					coords.append(exon[1])

			mincoord = min(coords)
			maxcoord = max(coords)

			#Write gene-level (utr) line
			outfh.write(('\t').join([chrm, 'UTRuniqueseqs', 'UTR', str(mincoord), str(maxcoord), '.', strand, 
				'.', 'ID={0}'.format(gene) + ';' + 'number_of_uniqueseqs={0}'.format(len(uniqueutrcoords[gene]))]) + '\n')

			#Write unique seqs line
			for uniqueseq in uniqueutrcoords[gene]:
				#Get min and max uniqueseq coords across all exons for this uniqueseq
				coords = []
				for exon in uniqueutrcoords[gene][uniqueseq]:
					coords.append(exon[0])
					coords.append(exon[1])
					mincoord = min(coords)
					maxcoord = max(coords)

				outfh.write(('\t').join([chrm, 'UTRuniqueseqs', 'uniqueUTR', str(mincoord), str(maxcoord), '.', strand,
					'.', 'ID={0}'.format(gene + '_' + 'uniqueUTR{0}'.format(str(uniqueseq))) + ';' + 'number_of_uniqueseqs={0}'.format(len(uniqueutrcoords[gene])) +
					';' + 'Parent={0}'.format(gene) + ';' + 'number_of_exons={0}'.format(len(uniqueutrcoords[gene][uniqueseq]))]) + '\n')

				#Write exon lines
				exoncounter = 0
				for exon in uniqueutrcoords[gene][uniqueseq]:
					exoncounter +=1
					outfh.write(('\t').join([chrm, 'UTRuniqueseqs', 'uniqueUTRexon', str(exon[0]), str(exon[1]), '.', strand,
						'.', 'ID={0}'.format(gene + '_' + 'uniqueUTR{0}'.format(str(uniqueseq)) + '_' + 'exon{0}'.format(exoncounter)) + ';' +
						'Parent={0}'.format(gene + '_' + 'uniqueUTR{0}'.format(str(uniqueseq))) + ';' + 'number_of_exons={0}'.format(len(uniqueutrcoords[gene][uniqueseq]))]) + '\n')


def getsequences(utrgff, genomefasta):
	#Take the coords of the unique seqs and get the sequences
	#Make gff database
	print('Indexing gff...')
	gff_fn = utrgff
	db_fn = os.path.abspath(utrgff) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	print('Indexing genome sequence...')
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta, 'rt'), 'fasta'))
	print('Done indexing!')

	utrs = db.features_of_type('UTR')
	print('Retrieving unique UTR seqeuences for {0} genes...'.format(sum(1 for _ in utrs)))

	#Remake generator
	utrs = db.features_of_type('UTR')
	
	with open('uniqueutrseqs.fa', 'w') as outfh:
		for utr in utrs:
			for uniqueUTR in db.children(utr, featuretype = 'uniqueUTR'):
				seq = ''
				for exon in db.children(uniqueUTR, featuretype = 'uniqueUTRexon', order_by = 'start'):
					if utr.strand == '+':
						exonseq = seq_dict[utr.chrom].seq[exon.start - 1 : exon.end].upper()
						seq += exonseq
					elif utr.strand == '-':
						exonseq = seq_dict[utr.chrom].seq[exon.start - 1 : exon.end].reverse_complement().upper()
						newseq = exonseq + seq
						seq = newseq

				outfh.write('>' + uniqueUTR.id + '\n' + str(seq) + '\n')

	os.remove(db_fn)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'Genome annotation in gff format.')
	parser.add_argument('--genomefasta', type = str, help = 'Genome sequence in gzipped fasta format.')
	args = parser.parse_args()

	posfactors = getpositionfactors(args.gff, 25) #lengthfactor of 25 so that transcript ends within 25 nt of each other are given the same m value
	#pickle.dump(posfactors, open('posfactors.pkl', 'wb'))
	#posfactors = pickle.load(open('posfactors.pkl', 'rb'))
	UTRcoords = getallUTRcoords(args.gff, posfactors)
	#pickle.dump(UTRcoords, open('utrcoords.pkl', 'wb'))
	#UTRcoords = pickle.load(open('utrcoords.pkl', 'rb'))
	cleanedUTRcoords = mergeUTRcoords(args.gff, UTRcoords)
	#pickle.dump(cleanedUTRcoords, open('cleanedUTRcoords.pkl', 'wb'))
	#cleanedUTRcoords = pickle.load(open('cleanedUTRcoords.pkl', 'rb'))
	uniqueutrcoords = getuniquecoords(args.gff, cleanedUTRcoords)
	makegff(args.gff, uniqueutrcoords)
	getsequences('uniqueutrcoords.gff', args.genomefasta)