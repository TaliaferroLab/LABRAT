#Given a bunch of unique UTR sequences (output of UTRuniqueseqs_PF.py), we want to see what sequences
#are prevalent in downstream unique regions relative to upstream unique regions, especially for a defined
#set of genes.  These represent potential binding sites of RBPs that are acting specifically on isoforms
#that utilize downstream polyA sites.

#We can count the abundance of kmers in each unique region, ending up with a "psi" score for each kmer in
#each gene, similar to the psi values calculated in LABRAT.  There's a problem here, though.  AU-rich kmers
#are naturally going to have higher psi values because of the lower GC content further and further into the 
#3' UTR.  To account for this, we will  compare the psi value of a kmer to the psi values of 50 GC- & 
#CpG-matched kmers to get a z-score for the kmer. We can then compare z-scores between two sets of genes,
#one of which is differentailly localized depending on APA site and one of which is not.

#python3
import itertools
import sys
from Bio import SeqIO
import argparse
import os
import numpy as np


#Count kmers in a single sequence
def countkmers(seq, k):
	#Get the total number of considered kmers.  This is different than len(seq) - k + 1 because we aren't considering every kmer (homopolymers, for example)
	consideredkmers = 0

	#If the length of the sequence is less than k, return None
	if len(seq) < k:
		return None, None

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
	seqcounter = 0

	kmerdict = {} #{ENSGENE : {UTRnumber : {kmer : count}}}
	consideredkmersdict = {} #{ENSGENE : {UTRnumber : total number of kmers}}
	with open(fasta, 'r') as infh:
		for record in SeqIO.parse(infh, 'fasta'):
			seqcounter +=1
			if seqcounter % 5000 == 0:
				print('Counting kmers in sequence {0}...'.format(seqcounter))
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

def cleankmerdicts(kmerdict, consideredkmersdict):
	#Sometimes a UTR is shorter than k, meaning that there are 0 considered kmers.
	#This is a problem because we later divide by the number of considered kmers, which would be 0
	#So we need to remove any UTR which has 'None' for kmers and/or consideredkmers
	#We then need to change the UTR number of any UTR after that UTR in that gene.  If UTR2 is short,
	#UTR3 becomes UTR2, UTR4 becomes UTR3, etc.
	cleankmerdict = {}
	cleanconsideredkmersdict = {}

	for gene in kmerdict:
		newgenedict_kmerdict = {} #{UTRnumber : {kmer : count}}}
		newgenedict_consideredkmersdict = {} #{UTRnumber : total number of kmers}}
		for UTRnumber in kmerdict[gene]:
			oldnumberofutrs = len(kmerdict[gene])
			#Remove this UTR if its value is None
			#k = UTRnumber, v = kmerdict
			newgenedict_kmerdict = {k : v for k, v in kmerdict[gene].items() if v}
			#k = UTRnumber, v = totalnumberofkmers
			newgenedict_consideredkmersdict = {k : v for k, v in consideredkmersdict[gene].items() if v}
			newnumberofutrs = len(newgenedict_kmerdict)
			if len(newgenedict_kmerdict) != len(newgenedict_consideredkmersdict):
				print('ERROR!!! ', gene)

		#If, after deleting all UTRs with None, there's only one left, delete the whole gene
		if newnumberofutrs == 1:
			continue
		#If there's at least two left, reorder them and then make new dictionaries based on the reordering (3 --> 2, 2 --> 1, etc.)
		elif newnumberofutrs >= 2:
			reorderedgenedict_kmer = {} #{UTRnumber : {kmer : count}}}
			reorderedgenedict_consideredkmers = {} # {UTRnumber : total number of kmers}}
			i = 0
			for kmerkey, consideredkmerkey in zip(sorted(newgenedict_kmerdict), sorted(newgenedict_consideredkmersdict)):
				reorderedgenedict_kmer[i] = newgenedict_kmerdict[kmerkey]
				reorderedgenedict_consideredkmers[i] = newgenedict_consideredkmersdict[consideredkmerkey]
				i +=1

			cleankmerdict[gene] = reorderedgenedict_kmer
			cleanconsideredkmersdict[gene] = reorderedgenedict_consideredkmers

	return cleankmerdict, cleanconsideredkmersdict



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
	genecounter = 0
	for gene in kmerdict:
		genecounter +=1
		if genecounter % 1000 == 0:
			print('Calculating weighted kmer densities for gene {0}...'.format(genecounter))
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

	return psis

def normalizepsis(psis):
	#We need to normalize psis based on their GC content.  This is done with a z-score normalization to all other kmers
	#with the same GC count

	#Get mean and sd for every GC count
	GCdict = {} #{GC count : [psis]}
	GCmeandict = {} #{GC count : meanpsi}
	GCsddict = {} #{GC count : psi sd}
	for kmer in psis:
		GCcount = kmer.count('G') + kmer.count('C')
		psi = psis[kmer]
		if GCcount not in GCdict:
			GCdict[GCcount] = [psi]
		else:
			GCdict[GCcount].append(psi)

	for GCcount in GCdict:
		psilist = GCdict[GCcount]
		psimean = np.mean(psilist)
		psisd = np.std(psilist)
		GCmeandict[GCcount] = psimean
		GCsddict[GCcount] = psisd

	#Calculate z-scores
	zs = {} #{kmer : zscore}
	for kmer in psis:
		psi = psis[kmer]
		GCcount = kmer.count('G') + kmer.count('C')
		psimean = GCmeandict[GCcount]
		psisd = GCsddict[GCcount]
		zscore = (psi - psimean) / psisd
		zs[kmer] = zscore

	return zs


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
		kmerdict, consideredkmersdict = cleankmerdicts(kmerdict, consideredkmersdict)
		psis = calculatepsi(kmerdict, consideredkmersdict, args.k)
		zs = normalizepsis(psis)
		#Write results
		fn = args.fasta + '.utrkmers'
		with open(fn, 'w') as outfh:
			for kmer in psis:
				outfh.write(kmer + '\t' + str(psis[kmer]) + '\t' + str(kmer.count('G') + kmer.count('C')) + '\t' + str(kmer.count('CG')) + '\t' + str(zs[kmer]) + '\n')

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

		os.remove('temp1.fasta')
		os.remove('temp2.fasta')



