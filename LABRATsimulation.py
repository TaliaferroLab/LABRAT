#Simulate an RNAseq dataset, then run LABRAT and other tools
#Two conditions, three reps each, we need psi values for genes in each replicate
#We are only going to consider genes that have 2 posfactors for simplicity

#For affected genes, there will be two ranges of psis for the two conditions
#One range will be randomly sampled from between 0.1 and 0.49.
#The other range will be randomly sampled from 0.5 to 0.9.
#Affected genes will have a further filter such that the minimum deltapsi for each condition is abs(0.1)
#Dpsi is condition2 - condition1

#If there are multiple transcripts associated with a posfactor, expression will be distributed among them randomly

#We want 1250 genes with positive dPSI, 1250 genes with negative dPSI, and 2500 genes with no change in psi

#usage
#python3 LABRATsimulation.py numberofposfactors.txt <alltranscriptcdnasinput.fa> <chosentxsoutput.fa> <txcountsoutput.txt> <genepsisout.txt> <txtpmsoutput.txt> <numberofreads>
#python3
import argparse
import random
import sys
import numpy as np
from Bio import SeqIO
from collections import defaultdict

def pickgenes(posfactorfile):
	#given a numberofposfactors file, get the genes that have 2 posfactors, and randomly pick 
	#250 to get positive dpsi, 250 to get negative dpsi, and 500 to get no change in dpsi
	pf2genes = [] #list of genes with 2 pfs
	with open(posfactorfile, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			if line[0] == 'Gene' or 'PAR' in line[0]: #PAR genes screw up parsing later
				continue
			gene = line[0].split('.')[0]
			pf = int(line[1])
			if pf == 2: #weird PAR genes, screws up parsing later
				pf2genes.append(gene)

	print('{0} genes have 2 posfactors.'.format(len(pf2genes)))

	#Now pick 250 for pos deltapsis
	posdpsigenes = random.sample(pf2genes, 1250)
	#Remove these genes from pf2genes
	pf2genes = [gene for gene in pf2genes if gene not in posdpsigenes]
	negdpsigenes = random.sample(pf2genes, 1250)
	pf2genes = [gene for gene in pf2genes if gene not in posdpsigenes]
	controldpsigenes = random.sample(pf2genes, 2500)

	print('Picked {0} positive deltapsi genes, {1} negative dpsi genes, and {2} control genes.'.format(len(posdpsigenes), len(negdpsigenes), len(controldpsigenes)))
	return posdpsigenes, negdpsigenes, controldpsigenes

def makepsis(sign):
	#Make dpsi values (3 per condition, 6 total) for a gene that should have a positive, negative, or control dpsi
	#Return [cond1psi, cond1psi, cond1psi, cond2psi, cond2psi, cond2psi]
	possiblelowpsis = list(np.linspace(0.1, 0.49, 40))
	possiblelowpsis = [round(psi, 2) for psi in possiblelowpsis]
	possiblehipsis = list(np.linspace(0.5, 0.9, 41))
	possiblehipsis = [round(psi, 2) for psi in possiblehipsis]

	if sign == 'positive' or sign == 'negative':
		minpsi = 0
		maxpsi = 0
		while (maxpsi - minpsi) < 0.1:
			lowpsis = random.choices(possiblelowpsis, k = 3)
			hipsis = random.choices(possiblehipsis, k = 3)
			minpsi = max(lowpsis)
			maxpsi = min(hipsis)

		if sign == 'positive':
			return lowpsis + hipsis
		elif sign == 'negative':
			return hipsis + lowpsis

	elif sign == 'control':
		mincond1psi = 0
		maxcond1psi = 1
		mincond2psi = 0
		maxcond2psi = 1
		while abs(maxcond1psi - mincond2psi > 0.25) or abs(maxcond2psi - mincond1psi > 0.25):
			possiblepsis = list(np.linspace(0.1, 0.9, 81))
			possiblepsis = [round(psi, 2) for psi in possiblepsis]
			psis = random.choices(possiblepsis, k = 6)
			mincond1psi = min(psis[0:3])
			maxcond1psi = max(psis[0:3])
			mincond2psi = min(psis[3:6])
			maxcond2psi = max(psis[3:6])

		return psis

def makegeneexpression(genes):
	#Given a list of genes (probably posdpsigenes, negdpsigenes, and controldpsigenes all summed together at once),
	#make fractional expression values for them (the fraction of the total pool that they take up, i.e. tpm / 1e6)
	#Return dictionary of {gene : expression}
	fracs = np.random.dirichlet(np.ones(len(genes)),size=1).tolist()[0]
	d = dict(zip(genes, fracs))

	return d

def makeupstreamanddownstreamtxs(posfactorfile, geneswecareabout):
	upstreamanddownstreamtxs = {} #{geneid : [[upstreamtxs], [downstreamtxs]]}
	with open(posfactorfile, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			if line[0] == 'Gene' or 'PAR' in line[0]:
				continue
			gene = line[0].split('.')[0]
			if gene in geneswecareabout:
				txs = line[2]
				upstreamtxs = txs.split('_')[0].split(',')
				downstreamtxs = txs.split('_')[1].split(',')
				upstreamtxs = [tx.split('.')[0] for tx in upstreamtxs]
				downstreamtxs = [tx.split('.')[0] for tx in downstreamtxs]
				upstreamanddownstreamtxs[gene] = [upstreamtxs, downstreamtxs]

	return upstreamanddownstreamtxs

def maketxexpression(upstreamtxs, downstreamtxs, psi, geneexp):
	txexpression = {} #{tx : expression}
	#For a single gene, give the expression of the whole gene and the relative expression of the
	#upstream txs (as a group) vs the downstream txs (as a group), assign expression levels to each tx
	upstreamgroupexpression = geneexp * (1 - psi)
	downstreamgroupexpression = geneexp * psi

	#Randomly assign expression to each tx in a group, such that their sum equals the desired expression of the group
	upstreamexps = np.random.dirichlet(np.ones(len(upstreamtxs)), size = 1).tolist()[0]
	#Scale so sum is upstreamgroupexpression
	upstreamexps = [exp * upstreamgroupexpression for exp in upstreamexps]

	downstreamexps = np.random.dirichlet(np.ones(len(downstreamtxs)), size = 1).tolist()[0]
	downstreamexps = [exp * downstreamgroupexpression for exp in downstreamexps]

	txs = upstreamtxs + downstreamtxs
	exps = upstreamexps + downstreamexps
	txexpression = dict(zip(txs, exps))

	return txexpression

def gettxlengths(txfasta):
	txlengths = {} #{tx : length}
	with open(txfasta, 'r') as infh:
		for record in SeqIO.parse(infh, 'fasta'):
			txname = str(record.id)
			txlength = len(record.seq)
			txlengths[txname] = txlength

	return txlengths

def scalefracsbylength(txexpression, txlengths):
	#tpm = counts / length....so counts = tpm * length, longer things are going to need more counts
	scaledtxexpression = {} #{tx : scaled expression (by length)}
	txs = []
	scaledexps = []
	for tx in txexpression:
		txs.append(tx)
		exp = txexpression[tx]
		txlength = txlengths[tx]
		#Going from counts to tpm by dividing counts by length is not exactly correct
		#What we want to do is divide counts by the *effective* length of the transcript
		#This is hard to estimate a priori, but empirically, this can be found by dividing the tx length by ~1.18
		scaledexp = (exp * txlength) / 1.18
		scaledexps.append(scaledexp)

	#OK now rescale these so that they sum to 1 (or pretty close)
	scaledexps_norm = [float(exp) / sum(scaledexps) for exp in scaledexps]
	scaledtxexpression = dict(zip(txs, scaledexps_norm))

	return scaledtxexpression

def maketxcounts(scaledtxexpression, totalreads):
	txcounts = {} #{tx : counts}
	for tx in scaledtxexpression:
		counts = int(round(scaledtxexpression[tx] * totalreads))
		txcounts[tx] = counts

	return txcounts


#pick genes
posdpsigenes, negdpsigenes, controldpsigenes = pickgenes(sys.argv[1])
allgenes = posdpsigenes + negdpsigenes + controldpsigenes

#assign psi values
psis = {} # {Gene : [cond1psi, cond1psi, cond1psi, cond2psi, cond2psi, cond2psi]}
for gene in posdpsigenes:
	p = makepsis('positive')
	psis[gene] = p
for gene in negdpsigenes:
	p = makepsis('negative')
	psis[gene] = p
for gene in controldpsigenes:
	p = makepsis('control')
	psis[gene] = p


#Make gene expression per gene, in "frac" units
geneexp_allgenes = makegeneexpression(allgenes)

#Get upstream and downstream txs
upstreamanddownstreamtxs = makeupstreamanddownstreamtxs(sys.argv[1], allgenes)


#Make tx expression, in "frac" units, which I think is essentially TPM / 1e6
#"Frac" of all transcripts in the experiment
#You have to do this once for each replicate (once for each psi in psis)
txreadcount = defaultdict(list) #{txid : [numreads in samp1.....numreads in samp6]}
txexpression_allsamples = defaultdict(list) #{tx : [sample1frac, sample2frac, ...]}
for i in list(range(0, 6)):
	print(i)
	txexpression = {} #{transcript : expression}
	for gene in allgenes:
		upstreamtxs = upstreamanddownstreamtxs[gene][0]
		downstreamtxs = upstreamanddownstreamtxs[gene][1]
		psi = psis[gene][i] 
		geneexp = geneexp_allgenes[gene]
		txexp = maketxexpression(upstreamtxs, downstreamtxs, psi, geneexp)
		for tx in txexp:
			txexpression_allsamples[tx].append(txexp[tx])
		txexpression.update(txexp)

	#For each tx, scale expression by txlength, then multiply scaled expression by total number of counts(reads) in the experiment (50 M)
	txlengths = gettxlengths(sys.argv[2])
	#Scale true frac by txlength in preparation for making readcounts
	scaledtxexpression = scalefracsbylength(txexpression, txlengths)
	#Scale true frac by txlength in preparation for making readcounts
	totalreads = int(sys.argv[7])
	txcounts = maketxcounts(scaledtxexpression, totalreads)
	for tx in txcounts:
		txreadcount[tx].append(txcounts[tx])

#output results
with open(sys.argv[2], 'r') as fainfh, open(sys.argv[3], 'w') as faoutfh, open(sys.argv[4], 'w') as expfh, open(sys.argv[5], 'w') as psisfh, open(sys.argv[6], 'w') as tpmfh:
	#make dictionary of all possible tx sequences
	alltxseqs = {} #{tx : seq}
	for record in SeqIO.parse(fainfh, 'fasta'):
		alltxseqs[str(record.id)] = str(record.seq)

	for tx in txreadcount:
		#write fasta sequence of tx
		faoutfh.write('>' + tx + '\n' + alltxseqs[tx] + '\n')
		#write readcounts
		readcounts = txreadcount[tx]
		readcounts = [str(r) for r in readcounts]
		expfh.write(('\t').join([tx] + readcounts) + '\n')

	#Write psi values
	psisfh.write(('\t').join(['Gene', 'Cond1Rep1psi', 'Cond1Rep2psi', 'Cond1Rep3psi', 'Cond2Rep1psi', 'Cond2Rep2psi', 'Cond2Rep3psi', 'genetype']) + '\n')
	for gene in psis:
		if gene in posdpsigenes:
			genetype = 'pos'
		elif gene in negdpsigenes:
			genetype = 'neg'
		else:
			genetype = 'ctrl'
		p = psis[gene]
		p = [str(psi) for psi in p]
		outlist = [gene] + p + [genetype]
		psisfh.write(('\t').join(outlist) + '\n')

	#Write expected tpms
	tpmfh.write(('\t').join(['transcript', 'Cond1Rep1tpm', 'Cond1Rep2tpm', 'Cond1Rep3tpm', 'Cond2Rep1tpm', 'Cond2Rep2tpm', 'Cond2Rep3tpm']) + '\n')
	for tx in txexpression_allsamples:
		tpms = [str(round(t * 1e6, 2)) for t in txexpression_allsamples[tx]]
		tpms = '\t'.join(tpms)
		tpmfh.write(tx + '\t' + tpms + '\n')

