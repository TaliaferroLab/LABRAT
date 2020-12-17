#Count aligned reads (as bam) at pA sites using bedtools
#Only reads that are 300 nt or less upstream from a pA site will be counted for that pA site

import gffutils
import os
import sys
import pysam
import pickle
import pandas as pd
from functools import reduce
import argparse

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

def getpositionfactorintervals(gff, posfactors):
    #The mean insert sizes for lexogen 3' end data is 200-300 nt. Looking at the Corley data
    #produced with this protocol on a genome browser, it looks like it is very often less than that.
    #So almost all of the reads associated with a particular polyA site should be within 300 nt of that site.

    #So for each transcript that is assigned to a position factor, define an interval as the last 300 nt of that
    #transcript.  Count how many reads lie in that window, and calculate psi values using those read counts.
    #posfactors = {} #{ENSMUSG : {ENSMUST : positionfactor}}
    posfactors_intervals = {} #{ENST : [chrm, windowstart, windowstop, strand]}

    print('Indexing gff...')
    gff_fn = gff
    db_fn = os.path.abspath(gff_fn) + '.db'
    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

    db = gffutils.FeatureDB(db_fn)
    print('Done indexing!')

    genecounter = 0
    for gene in posfactors:
        genecounter +=1
        if genecounter % 1000 == 0:
            print('Gene {0} of {1}...'.format(genecounter, len(posfactors)))
        for txid in posfactors[gene]:
            txexons = [] #list of exonic coordinates for this transcript
            tx = db[txid]
            for exon in db.children(tx, featuretype = 'exon', order_by = 'start'):
                txexons += list(range(exon.start, exon.end + 1))

            if tx.strand == '-':
                txexons = list(reversed(txexons))

            if len(txexons) > 300:
                windowstart = txexons[-300]
            else:
                windowstart = txexons[0]
            windowstop = txexons[-1]

            if tx.strand == '+':
                posfactors_intervals[txid] = [str(tx.chrom), int(windowstart), int(windowstop), str(tx.strand)]
            elif tx.strand == '-':
                posfactors_intervals[txid] = [str(tx.chrom), int(windowstop), int(windowstart), str(tx.strand)]

    return posfactors_intervals

def calculatepsi(posfactors, posfactors_intervals, bamfile):
    #Take a bam of aligned reads and calculate psi by counting the number of reads
    #in the predefined 300 nt windows for each transcript

    txcounts = {} #{transcriptid : readcounts} 
    genecounts = {} #{geneid : [transcript counts]} (unscaled)
    posfactorgenecounts = {} #{geneid : [transcriptcounts scaled by posfactor]} (scaled)
    psis = {} #{geneid : psi}

    bam = pysam.AlignmentFile(bamfile, 'rb')

    print('Calculating psi for {0}...'.format(os.path.basename(bamfile)))
    genecounter = 0
    for gene in posfactors:
        genecounter +=1
        if genecounter % 1000 == 0:
            print('Gene {0} of {1}...'.format(genecounter, len(posfactors)))
        genecounts[gene] = []
        posfactorgenecounts[gene] = []
        for tx in posfactors[gene]:
            posfactor = posfactors[gene][tx]

            window = posfactors_intervals[tx]
            windowchrm = window[0]
            windowstart = window[1]
            windowstop = window[2]
            strand = window[3]

            overlappingreads = bam.fetch(windowchrm, windowstart, windowstop)
            readsinwindow = []
            for read in overlappingreads:
                if read.is_reverse:
                    readstrand = '-'
                elif not read.is_reverse:
                    readstrand = '+'

                if strand == readstrand:
                    readsinwindow.append(read)

            counts = len(readsinwindow)
            scaledcounts = counts * posfactor
            genecounts[gene].append(counts)
            posfactorgenecounts[gene].append(scaledcounts)

    #Calculate psi values
    for gene in posfactorgenecounts:
        totalcounts = sum(genecounts[gene])
        scaledcounts = sum(posfactorgenecounts[gene])
        if totalcounts < 100:
            psi = 'NA'
        else:
            psi = scaledcounts / totalcounts
            psis[gene] = round(psi, 3)

    #turn into df
    df = pd.DataFrame.from_dict(psis, orient = 'index', columns = [os.path.basename(bamfile)])
    
    return df

def domanybams(bamdir, posfactors, posfactors_intervals):
	dfs = []
	for f in os.listdir(bamdir):
		fp = os.path.join(os.path.abspath(bamdir), f)
		if fp.endswith('.bam'):
			df = calculatepsi(posfactors, posfactors_intervals, fp)
			dfs.append(df)

	#merge dfs 
	mergeddf = reduce(lambda df1, df2: pd.merge(df1, df2, left_index = True, right_index = True, how = 'inner'), dfs)

	return mergeddf    



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', help = 'Genome annotation in gff format.')
	parser.add_argument('--bamdir', help = 'Directory of bam files to quantify.')
	parser.add_argument('--output', help = 'Output file for psi values.')
	args = parser.parse_args()

	posfactors = getpositionfactors(args.gff, 25)
	positionfactorintervals = getpositionfactorintervals(args.gff, posfactors)
	df = domanybams(args.bamdir, posfactors, positionfactorintervals)
	df.to_csv(args.output, sep = '\t', header = True, index = True, index_label = 'Gene')


            









