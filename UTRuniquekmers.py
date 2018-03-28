#Given a bunch of unique UTR sequences (output of UTRuniqueseqs_PF.py), we want to see what sequences
#are prevalent in downstream unique regions relative to upstream unique regions, especially for a defined
#set of genes.  These represent potential binding sites of RBPs that are acting specifically on isoforms
#that utilize downstream polyA sites.

#We can count the abundance of kmers in each unique region, ending up with a "psi" score for each kmer in
#each gene, similar to the psi values calculated in LABRAT.  There's a problem here, though.  AU-rich kmers
#are naturally going to have higher psi values because of the lower GC content further and further into the 
#3' UTR.  To account for this, we will (1) compare psi values for kmers in two different sets of genes, and
#(2) compare the psi value of a kmer to the psi values of 50 GC- & CpG-matched kmers to get a z-score for the kmer.

#python3
import itertools
import sys
from Bio import SeqIO


#Count kmers in a single sequence
def countkmers(seq, k):
	kmers = {} #{kmer : count}
	seq = seq.replace('T', 'U').upper()

	#Make all possible kmers
	k = int(k)
	bases = ['A', 'U', 'G', 'C']
	allkmers = [''.join(x) for x in itertools.product(bases, repeat = k)]

	#Add all kmers to kmers dictionary
	for kmer in allkmers:
		kmers[kmer] = 0

	currentkmer = ''
	previouskmer = ''
	currentkmercount = 0
	for i in range(len(seq) -k + 1):
		kmer = seq[i:i+k]
		#Count homopolymeric stretches as only one instance of kmer, but allow nonoverlapping instances of
		#kmer in homopolymers
		if ((kmer != currentkmer and kmer != previouskmer) or currentkmercount >= k) and 'N' not in kmer: 
			kmers[kmer] +=1
			#To avoid over-counting dinucleotide repeats (UCUCUCUCUCUCUCUC, for example), keep track of what the kmer 2 rounds ago was.
			#currentkmer is the kmer from last round, previous kmer is from 2 rounds ago.
			previouskmer = currentkmer
			currentkmer = kmer
			currentkmercount = 1
		else:
			currentkmercount +=1

	return kmers

#Count kmers in all seqs in the uniqueutrfasta
def countkmers_fasta(fasta, k):
	k = int(k)

	kmerdict = {} #{ENSGENE : {UTRnumber : {kmer : count}}}
	with open(fasta, 'r') as infh:
		for record in SeqIO.parse(infh, 'fasta'):
			genename = record.id.split('.')[0]
			UTRnumber = int(record.id.split('uniqueUTR')[1])
			kmers = countkmers(str(record.seq), k)
			if genename not in kmerdict:
				kmerdict[genename] = {}
			kmerdict[genename][UTRnumber] = kmers

	return kmerdict

#Calculate psi values for each kmer given a dictionary like kmerdict {ENSGENE: {UTRnumber : {kmers}}}
def calculatepsi(kmerdict, k):
	#kmerdict = {} #{ENSGENE : {UTRnumber : {kmer : count}}}
	psis = {} #{kmer : psi}
	utrcounts = {} #{genename : number of utrs}
	rawkmercounts = {} #{kmer : raw number of times its seen in any utr}
	weightedkmercounts = {} #{kmer : weighted number of times its seen in any utr}

	#Get the number of utrs for each gene
	for gene in kmerdict:
		utrcounts[gene] = len(kmerdict[gene])

	#Get the raw number of times and weighted number of times each kmer is seen
	for gene in kmerdict:
		for utrnumber in kmerdict[gene]:
			for kmer in kmerdict[gene][utrnumber]:
				rawcount = kmerdict[gene][utrnumber][kmer]
				weightingfactor = float(utrnumber / (utrcounts[gene] - 1)) #just like calculating psis in LABRAT
				weightedcounts = rawcount * weightingfactor
				
				if kmer not in rawkmercounts:
					rawkmercounts[kmer] = rawcount
				else:
					rawkmercounts[kmer] += rawcount

				if kmer not in weightedkmercounts:
					weightedkmercounts[kmer] = weightedcounts
				else:
					weightedkmercounts[kmer] += weightedcounts

	#psi for each kmer is sum of all weighted counts divided by sum of all raw counts
	#Make all possible kmers
	k = int(k)
	bases = ['A', 'U', 'G', 'C']
	allkmers = [''.join(x) for x in itertools.product(bases, repeat = k)]

	for kmer in allkmers:
		if kmer in rawkmercounts:
			psi = weightedkmercounts[kmer] / rawkmercounts[kmer]
		else:
			psi = 'NA'
		psis[kmer] = psi

	print(psis)
	print(min(psis.values()), max(psis.values()))

	

#print(countkmers('TACGACGTACTACTGACTGACTGCAGTACGTACGTACTGACACACACTGACGTACTGACTGACTACTCA', 5))
kmerdict = countkmers_fasta(sys.argv[1], 5)
calculatepsi(kmerdict, 5)