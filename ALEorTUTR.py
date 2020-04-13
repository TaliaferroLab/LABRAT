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
import pickle
import operator

with open('/Users/mtaliaferro/Desktop/posfactors.hg38.pkl', 'rb') as infh:
	posfactors = pickle.load(infh)

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

	ales = 0
	tutrs = 0
	mixed = 0
	error = 0
	with open('/Users/mtaliaferro/Desktop/Threeprimeendtypes.hg38.txt', 'w') as outfh:
		for gene in genetypes:
			outfh.write(gene + '\t' + genetypes[gene] + '\n')
			if genetypes[gene] == 'ALE':
				ales +=1
			elif genetypes[gene] == 'TUTR':
				tutrs +=1
			elif genetypes[gene] == 'mixed':
				mixed +=1
			elif genetypes[gene] == 'ERROR':
				error +=1

	print len(genetypes), ales, tutrs, mixed, error

	return genetypes

exoniccoords = getexoniccoords(posfactors, sys.argv[1])
classifygenes(exoniccoords, sys.argv[1])
