import gffutils
import os
import sys
from Bio import SeqIO
import gzip

#We want to know if the genes that use alternative polyadenylation to go to the ER are more or less likely to have an ER signal sequence
#than those that do not.  We first need to get the amino acid sequences for all transcripts for all genes.  We will then take those sequences
#and look for signal sequences using SignalP 4.1.  This information will then be brought back and those genes that have transcripts with 
#position factors (i.e. those that make it into the alt polyA analysis), will be queried.

def writeallAA(gff, genomefasta):
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

	
	aaseqs = {} #{ensg_enst : aaseq}
	genes = db.features_of_type('gene')

	genecounter = 0
	for gene in genes:
		genecounter +=1
		if genecounter % 5000 == 0:
			print 'Gene {0}...'.format(genecounter)
		
		#Only protein-coding genes
		if 'protein_coding' not in gene.attributes['gene_type']:
			continue

		for transcript in db.children(gene, featuretype = 'transcript'):
			#Only protein-coding transcripts
			if 'protein_coding' not in transcript.attributes['transcript_type']:
				continue
			#Stitch together CDS pieces
			cdsseq = ''

			if transcript.strand == '+':
				for cds in db.children(transcript, featuretype = 'CDS', order_by = 'start'):
					try:
						seq = seq_dict[transcript.chrom].seq[cds.start - 1 : cds.end]
					except KeyError: #This chromosome isn't in the fasta
						print 'Chromosome {0} is not in the genome sequence fasta.'.format(transcript.chrom)
						continue
					cdsseq += seq

			elif transcript.strand == '-':
				for cds in db.children(transcript, featuretype = 'CDS', order_by = 'start', reverse = True):
					try:
						seq = seq_dict[transcript.chrom].seq[cds.start - 1 : cds.end].reverse_complement()
					except KeyError: #This chromosome isn't in the fasta
						print 'Chromosome {0} is not in the genome sequence fasta.'.format(transcript.chrom)
						continue
					cdsseq += seq

			if cdsseq:
				aaseq = str(cdsseq.translate())
				aaseqs[gene.id + '_' + transcript.id] = aaseq

	fcounter = 1
	txcounter = 0
	totaltxcounter = 0
	longtxcounter = 0
	outfh = open('AAseqs_{0}.fa'.format(fcounter), 'w')
	for tx in aaseqs:
		txcounter +=1
		totaltxcounter +=1
		if txcounter == 10000:
			txcounter = 1
			fcounter += 1
			outfh.close()
			outfh = open('AAseqs_{0}.fa'.format(fcounter), 'w')
		if len(aaseqs[tx]) < 8000: 
			outfh.write('>' + tx + '\n' + aaseqs[tx] + '\n')
		else: 
			longtxcounter +=1

	outfh.close()

	print totaltxcounter, longtxcounter


writeallAA(sys.argv[1], sys.argv[2])


