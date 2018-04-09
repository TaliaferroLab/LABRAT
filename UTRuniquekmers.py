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
import argparse


#Count kmers in a single sequence
def countkmers(seq, k):
	#Get the total number of considered kmers.  This is different than len(seq) - k + 1 because we aren't considering every kmer (homopolymers, for example)
	consideredkmers = 0

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
			consideredkmers += 1
		else:
			currentkmercount +=1

	return kmers, consideredkmers

#Count kmers in all seqs in the uniqueutrfasta
def countkmers_fasta(fasta, k):
	k = int(k)

	kmerdict = {} #{ENSGENE : {UTRnumber : {kmer : count}}}
	consideredkmersdict = {} #{ENSGENE : {UTRnumber : total number of kmers}}
	with open(fasta, 'r') as infh:
		for record in SeqIO.parse(infh, 'fasta'):
			genename = record.id.split('.')[0]
			UTRnumber = int(record.id.split('uniqueUTR')[1])
			kmers, consideredkmers = countkmers(str(record.seq), k)
			
			if genename not in kmerdict:
				kmerdict[genename] = {}
			kmerdict[genename][UTRnumber] = kmers

			if genename not in consideredkmersdict:
				consideredkmersdict[genename] = {}
			consideredkmersdict[genename][UTRnumber] = consideredkmers

	return kmerdict, consideredkmersdict

#Calculate psi values for each kmer given a dictionary like kmerdict {ENSGENE: {UTRnumber : {kmers}}}
def calculatepsi(kmerdict, consideredkmersdict, k):
	#kmerdict = {} #{ENSGENE : {UTRnumber : {kmer : count}}}
	psis = {} #{kmer : psi}
	utrcounts = {} #{genename : number of utrs}
	rawkmerdensities = {} #{kmer : [raw densities across all UTRs]}
	weightedkmerdensities = {} #{kmer : [weighted densities across all UTRs]}

	#Get the number of utrs for each gene
	for gene in kmerdict:
		utrcounts[gene] = len(kmerdict[gene])

	#Get the raw kmer density for each kmer/UTR and the weighted kmer density for each kmer/UTR
	for gene in kmerdict:
		for utrnumber in kmerdict[gene]:
			for kmer in kmerdict[gene][utrnumber]:
				rawcount = kmerdict[gene][utrnumber][kmer]
				weightingfactor = float(utrnumber / (utrcounts[gene] - 1)) #just like calculating psis in LABRAT
				totalkmers = consideredkmersdict[gene][utrnumber]
				rawdensity = rawcount / totalkmers
				weightedkmerdensity = rawdensity * weightingfactor
				
				if kmer not in rawkmerdensities:
					rawkmerdensities[kmer] = [rawdensity]
				else:
					rawkmerdensities[kmer].append(rawdensity)

				if kmer not in weightedkmerdensities:
					weightedkmerdensities[kmer] = [weightedkmerdensity]
				else:
					weightedkmerdensities[kmer].append(weightedkmerdensity)

	#psi for each kmer is sum of all weighted counts divided by sum of all raw counts
	#Make all possible kmers
	k = int(k)
	bases = ['A', 'U', 'G', 'C']
	allkmers = [''.join(x) for x in itertools.product(bases, repeat = k)]

	for kmer in allkmers:
		if kmer in rawkmerdensities:
			psi = sum(weightedkmerdensities[kmer]) / sum(rawkmerdensities[kmer])
		else:
			psi = 'NA'
		psis[kmer] = psi

	with open('/Users/mtaliaferro/Desktop/temp2.txt', 'w') as outfh:
		for kmer in psis:
			outfh.write(kmer + '\t' + str(psis[kmer]) + '\t' + str(kmer.count('G') + kmer.count('C')) + '\n')

	return psis


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', type = str, help = 'Unique UTR sequences in fasta format.  Output of UTRuniqueseqs_PF.py.')
	parser.add_argument('-k', type = int, help = 'Kmer length.')
	parser.add_argument('--controlmode', choices = ['twogenelists', 'kmerzscore'], help = 'How to control for GC content?')
	parser.add_argument('--genelist1', type = str, help = 'For use with controlmode twogenelists. List of genes in set 1.')
	parser.add_argument('--genelist2', type = str, help = 'For use with controlmode twogenelists. List of genes in set 2.')
	args = parser.parse_args()

	#If we aren't using a control
	if not args.controlmode:
		kmerdict, consideredkmersdict = countkmers_fasta(args.fasta, args.k)
		psis = calculatepsi(kmerdict, consideredkmersdict, args.k)
		#Write results
		fn = args.fasta + '.utrkmers'
		with open(fn, 'w') as outfh:
			for kmer in psis:
				outfh.write(kmer + '\t' + str(psis[kmer]) + '\n')

	#Genelists are list of ensembl IDs without dots
	if args.controlmode == 'twogenelists':
		genes1 = []
		genes2 = []
		with open(args.genelist1, 'r') as infh:
			for line in infh:
				line = line.strip()
				genes1.append(line)

		with open(args.geneslist2, 'r') as infh:
			for line in infh:
				line = line.strip()
				genes2.append(line)

		with open(args.fasta, 'r') as infh:
			records1 = []
			records2 = []
			for record in SeqIO.parse(infh, 'fasta'):
				genename = record.id.split('.')[0]
				if genename in genes1:
					records1.append(record)
				if genename in genes2:
					records2.append(record)

		#Make temporary fasta files for the two gene lists
		SeqIO.write(records1, 'temp1.fasta', 'fasta')
		SeqIO.write(records2, 'temp2.fasta', 'fasta')

		#Run functions for two fasta files
		kmerdict, consideredkmersdict = countkmers_fasta('temp1.fasta', args.k)
		psis1 = calculatepsi(kmerdict, consideredkmersdict, args.k)

		kmerdict, consideredkmersdict = countkmers_fasta('temp2.fasta', args.k)
		psis2 = calculatepsi(kmerdict, consideredkmersdict, args.k)

		#Calculate dpsis
		dpsis = {}
		for kmer in psis1:
			if psis1[kmer] != 'NA' and psis2[kmer] != 'NA':
				dpsi = psis2[kmer] - psis1[kmer]
				dpsis[kmer] = dpsi
			else:
				dpsis[kmer] = 'NA'

		fn = args.fasta + '.utrkmers'
		with open(fn, 'w') as outfh:
			outfh.write('kmer' + '\t' + 'genes1psi' + '\t' + 'genes2psi' + '\t' + 'dpsi' + '\n')
			for kmer in dpsis:
				outfh.write(kmer + '\t' + str(psis1[kmer]) + '\t' + str(psis2[kmer]) + '\t' + str(dpsis[kmer]) + '\n')



