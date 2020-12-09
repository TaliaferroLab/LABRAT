#Application of LABRAT to alevin-quantified scRNAseq data

from struct import Struct
import numpy as np
import pandas as pd
import gzip
import sys
import os
from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
import gffutils
import pickle
import argparse


def read_quants_bin(base_dir, clipped=False, mtype="data"):
    '''
    Read the quants Sparse Binary output of Alevin and generates a dataframe
    Parameters
    ----------
    base_location: string
        Path to the folder containing the output of the alevin run
    clipped: bool (default False)
        Clip off all zero rows and columns
    mtype: "[data(default), tier, var, mean]"
        Alevin's matrix type to load into memory

    This function is heavily based on a function in https://github.com/k3yavi/vpolo/blob/master/vpolo/alevin/parser.py
    '''
    if not os.path.isdir(base_dir):
        print("{} is not a directory".format( base_dir ))
        sys.exit(1)

    base_location = os.path.join(base_dir, "alevin")
    
    #print(base_location)
    if not os.path.exists(base_location):
        print("{} directory doesn't exist".format( base_location ))
        sys.exit(1)

    data_type = "f"
    if mtype == "data":
        quant_file = os.path.join(base_location, "quants_mat.gz")
    elif mtype == "tier":
        data_type = "B"
        quant_file = os.path.join(base_location, "quants_tier_mat.gz")
    elif mtype == "mean":
        quant_file = os.path.join(base_location, "quants_mean_mat.gz")
    elif mtype == "var":
        quant_file = os.path.join(base_location, "quants_var_mat.gz")
    else:
        print("wrong mtype:".format( mtype ))
        sys.exit(1)

    if not os.path.exists(quant_file):
        print("quant file {} doesn't exist".format( quant_file ))
        sys.exit(1)

    if mtype in ["mean", "var"]:
        cb_file = os.path.join(base_location, "quants_boot_rows.txt")
    else:
        cb_file = os.path.join(base_location, "quants_mat_rows.txt")

    if not os.path.exists(cb_file):
        print("quant file's index: {} doesn't exist".format( cb_file ))
        sys.exit(1)

    gene_file = os.path.join(base_location, "quants_mat_cols.txt")
    
    if not os.path.exists(gene_file):
        print("quant file's header: {} doesn't exist".format(gene_file))
        sys.exit(1)

    cb_names = pd.read_csv(cb_file, header=None)[0].values
    gene_names = pd.read_csv(gene_file, header=None)[0].values
    num_genes = len(gene_names)
    num_cbs = len(cb_names)
    num_entries = int(np.ceil(num_genes/8))


    with gzip.open( quant_file ) as f:
        line_count = 0
        tot_umi_count = 0
        umi_matrix = []

        header_struct = Struct( "B" * num_entries)
        while True:
            line_count += 1
            if line_count%100 == 0:
                print ("\r Done reading " + str(line_count) + " cells", end= "")
                sys.stdout.flush()
            try:
                num_exp_genes = 0
                exp_counts = header_struct.unpack_from( f.read(header_struct.size) )
                for exp_count in exp_counts:
                    num_exp_genes += bin(exp_count).count("1")

                data_struct = Struct( data_type * num_exp_genes)
                sparse_cell_counts_vec = list(data_struct.unpack_from( f.read(data_struct.size) ))[::-1]
                cell_umi_counts = sum(sparse_cell_counts_vec)

            except:
                print ("\nRead total " + str(line_count-1) + " cells")
                print ("Found total " + str(round(tot_umi_count)) + " reads")
                break

            if cell_umi_counts > 0.0:
                tot_umi_count += cell_umi_counts

                cell_counts_vec = []
                for exp_count in exp_counts:
                    for bit in format(exp_count, '08b'):
                        if len(cell_counts_vec) >= num_genes:
                            break

                        if bit == '0':
                            cell_counts_vec.append(0.0)
                        else:
                            abund = sparse_cell_counts_vec.pop()
                            cell_counts_vec.append(abund)

                if len(sparse_cell_counts_vec) > 0:
                    print("Failure in consumption of data")
                    print("left with {} entry(ies)".format(len(sparse_cell_counts_vec)))
                umi_matrix.append( cell_counts_vec )
            else:
                print("Found a CB with no read count, something is wrong")
                sys.exit(1)


    alv = pd.DataFrame(umi_matrix)
    alv.columns = gene_names
    alv.index = cb_names
    if clipped:
        alv = alv.loc[:, (alv != 0).any(axis=0)]

    barcodes = list(alv.index.values)
    barcodes = [os.path.basename(os.path.abspath(base_dir)) + '_' + b for b in barcodes]
    alv = alv.set_index([barcodes])

    #for testing, can use alv.head() here to only return first n cells
	return alv

def createbigdf(alevindir):
    #Given a directory (alevindir) that contains alevin output subdirectories
    #concatentate all alevin outputs into one big df

    #directory structure:
    # alevindir
    #     |
    #     |---sample1
    #     |      |
    #     |      |----alevin
    #     |              |----quants_mat.gz, etc.
    #     |---sample2
    #     |---sample3
    #     etc.

    #index of big df will be sampleid + '_' + cell barcode
    # e.g. sample1_CAGTCTCGTTGCTCAG

    alevindir = os.path.abspath(alevindir)
    dfs = [] #list of sample dfs
    samples = os.listdir(alevindir)
    for samp in samples:
        samppath = os.path.join(alevindir, samp)
        if not os.path.isdir(samppath):
            continue
        print('Reading data in {0}...'.format(samppath))
        p = os.path.join(alevindir, samp)
        df = read_quants_bin(p, False, 'data')
        dfs.append(df)

    bigdf = pd.concat(dfs)
    bigdf = bigdf.reset_index()
    bigdf = bigdf.rename(columns = {'index': 'cellid'})

    return bigdf

def addclusters(bigdf, clusters):
	#Add labels of cluster ID to df of counts
	clusterdf = pd.read_csv(clusters, sep = '\t', header = 0, index_col = None)
	clusterdf = clusterdf.rename(columns = {'sample' : 'cellid'})
	newbigdf = pd.merge(clusterdf, bigdf, how = 'inner', on = 'cellid')
	
	return newbigdf

def collapseclusters(bigdf):
	#Sum transcript counts across every cell in a cluster
	conditions = list(set(bigdf['condition'].tolist()))
	allcountsums = [] #list of countsums dictionaries, in order of conditions
	for condition in conditions:
		countsums = {} #{transcriptid : countsum}
		countsums['condition'] = condition
		df = bigdf.loc[bigdf['condition'] == condition]
		colnames = list(df)
		#select only transcriptID columns
		colnames = [c for c in colnames if c != 'cellid' and c != 'condition']
		for c in colnames:
			counts = df[c].tolist()
			sumcounts = sum(counts)
			countsums[c] = [sumcounts]
		allcountsums.append(countsums)

	#Turn the countsum dictionaries into dfs
	dfs = []
	for countsum in allcountsums:
		df = pd.DataFrame.from_dict(countsum)
		dfs.append(df)

	#Concatenate dfs
	countsumdf = pd.concat(dfs)

	return countsumdf

def calculatepsi_fromcollapsedclusters(countsumdf, posfactors):
	#Calculate psi values for every gene in posfactors
	#This function takes a df that has one row per cluster/condition
	#that was produced by collapseclusters().

	conditions = countsumdf['condition'].tolist()
	allpsidicts = [] #list of all psidicts

	for condition in conditions:
		df = countsumdf.loc[countsumdf['condition'] == condition]
		psidict = {} #{gene : psi}
		psidict['condition'] = condition
		for gene in posfactors:
			geneid = gene.split('.')[0]
			rawcounts = []
			scaledcounts = []
			for tx in posfactors[gene]:
				txid = tx.split('.')[0]
				m = posfactors[gene][tx]
				counts = df[txid].tolist()[0]
				scaled_counts = counts * m
				rawcounts.append(counts)
				scaledcounts.append(scaled_counts)

			if sum(rawcounts) <= 0: #add count filter here
				psi = 'NA'
			else:
				psi = sum(scaledcounts) / sum(rawcounts)
			
			psidict[geneid] = [psi]
		
		allpsidicts.append(psidict)

	#Turn the psidicts into dfs
	dfs = []
	for psidict in allpsidicts:
		df = pd.DataFrame.from_dict(psidict)
		dfs.append(df)

	#Concatenate dfs
	psidf = pd.concat(dfs)

	return psidf

def calculatedpsi_fromcollapsedclusters(psidf, conditionA, conditionB):
	#Given a dataframe of psi values (from calculatepsi_fromcollapsedclusters()),
	#calculate a dpsi value for each gene

	#dpsi is defined as conditionB - conditionA
	
	dpsis = {} #{gene : [condApsi, condBpsi, dpsi]}
	colnames = list(psidf.columns.values)
	colnames = [c for c in colnames if c != 'condition']

	for gene in colnames:
		condApsi = psidf.loc[psidf['condition'] == conditionA].at[0, gene]
		condBpsi = psidf.loc[psidf['condition'] == conditionB].at[0, gene]
		if condApsi == 'NA' or condBpsi == 'NA':
			dpsi = 'NA'
		else:
			condApsi = round(condApsi, 3)
			condBpsi = round(condBpsi, 3)
			dpsi = round(condBpsi - condApsi, 3)
		dpsis[gene] = [condApsi, condBpsi, dpsi]

	#Turn into df
	df = pd.DataFrame.from_dict(dpsis, orient = 'index')
	df = df.reset_index()
	df.columns = ['gene', '{0}_psi'.format(conditionA), '{0}_psi'.format(conditionB), 'deltapsi']
	
	return df

def dotest_bootstrapclusters(bigdf, posfactors, fractosample, numberofsamples, conditionA, conditionB):
	#Given a dataframe of counts per cell (usually bigdf, produced by addclusters())
	#for each condition (cluster), randomly select x% of the cells. Then collapse the clusters
	#and calculate psivalues using collapseclusters() and calculatepsi_fromcollapsedclusters()
	#Record psi for each gene in each cluster.
	#Repeat n times
	#Compare psi values across clusters to derive p value

	condAdf = bigdf.loc[bigdf['condition'] == conditionA]
	condBdf = bigdf.loc[bigdf['condition'] == conditionB]

	condApsis = {} #{gene : [psi values from samples]}
	condBpsis = {} #{gene : [psi values from samples]}

	for i in range(numberofsamples):
		print('Sample ' + str(i + 1))
		#sample rows (cells)
		condAsamp = condAdf.sample(frac = fractosample, replace = False, axis = 'index')
		condBsamp = condBdf.sample(frac = fractosample, replace = False, axis = 'index')
		#Put the condA and condB back together again
		df = pd.concat([condAsamp, condBsamp], axis = 0, ignore_index = True)
		collapseddf = collapseclusters(df)
		psidf = calculatepsi_fromcollapsedclusters(collapseddf, posfactors)
		
		#Record psi values for each gene
		colnames = list(psidf.columns.values)
		genes = [c for c in colnames if c != 'condition']
		for gene in genes:
			condApsi = psidf.loc[psidf['condition'] == conditionA].at[0, gene]
			condBpsi = psidf.loc[psidf['condition'] == conditionB].at[0, gene]
			if gene not in condApsis:
				condApsis[gene] = [condApsi]
			else:
				condApsis[gene].append(condApsi)
			if gene not in condBpsis:
				condBpsis[gene] = [condBpsi]
			else:
				condBpsis[gene].append(condBpsi)

	#Test if sampled psis are different across conditions
	#wilcoxon rank sum test
	pvals = {} #{gene : pval}
	for gene in genes:
		Apsis = condApsis[gene]
		Bpsis = condBpsis[gene]
		Apsis = [psi for psi in Apsis if psi != 'NA']
		Bpsis = [psi for psi in Bpsis if psi != 'NA']

		#if either is just NAs, p = NA
		if not Apsis or not Bpsis:
			p = 'NA'
		else:
			p = ranksums(Apsis, Bpsis)[1]
			p = round(float('{:.2e}'.format(p)), 3)			
		pvals[gene] = p

	
	#Correct pvalues using BH method, but only using pvalues that are not NA
	genes = list(pvals.keys())
	pvalues = list(pvals.values())
	pvaluesformultitest = [pval for pval in pvalues if str(pval) != 'NA']
	fdrs = list(multipletests(pvaluesformultitest, method = 'fdr_bh')[1])
	fdrs = [float('{:.2e}'.format(fdr)) for fdr in fdrs]

	#Now we need to incorporate the places where p = NA into the list of FDRs (also as NA)
	fdrswithnas = []
	fdrindex = 0
	for pvalue in pvalues:
		if str(pvalue) != 'NA':
			fdrswithnas.append(fdrs[fdrindex])
			fdrindex +=1
		elif str(pvalue) == 'NA':
			fdrswithnas.append('NA')
	fdrswithnas = [[fdr] for fdr in fdrswithnas] #fdr must be in list for conversion to df

	correctedpvals = dict(zip(genes, fdrswithnas))

	#make df of pval and FDR
	pvalues = [[p] for p in pvalues] #p must be in list for conversion to df
	pvals = dict(zip(genes, pvalues))
	rawpdf = pd.DataFrame.from_dict(pvals, orient = 'index')
	rawpdf.columns = ['raw_pval']
	fdrdf = pd.DataFrame.from_dict(correctedpvals, orient = 'index')
	fdrdf.columns = ['fdr']

	#Merge rawpdf and fdrdf
	df = rawpdf.merge(fdrdf, how = 'inner', left_index = True, right_index = True)
	df = df.reset_index()
	df.columns = ['gene', 'raw_pval', 'FDR']

	return df	





def calculatepsi_cellbycell(bigdf, posfactors):
	#Calculate psi values for every gene in every cell
	#In practice, this is not going to be super useful because many
	#gene level counts are going to be 0, giving a psi of NA
	#This function takes a df of transcript counts per cell (produced by createbigdf() and addclusters())
	
	#Make a small df of just cellid and condition
	conditiondf = bigdf[['cellid', 'condition']]
	
	allpsidicts = {} #{cellid : {gene : psi}}
	cellcount = 0

	for idx, row in bigdf.iterrows():
		cellcount +=1
		if cellcount % 10 == 0:
			print('Calculating psi for cell {0}'.format(cellcount))
		cellid = row['cellid']
		allpsidicts[cellid] = {}
		for gene in posfactors:
			geneid = gene.split('.')[0]
			rawcounts = []
			scaledcounts = []
			for tx in posfactors[gene]:
				txid = tx.split('.')[0]
				m = posfactors[gene][tx]
				counts = row[txid]
				scaled_counts = counts * m
				rawcounts.append(counts)
				scaledcounts.append(scaled_counts)

			if sum(rawcounts) <= 0: #add counter filter here
				psi = 'NA'
			else:
				psi = sum(scaledcounts) / sum(rawcounts)
				psi = round(psi, 3)

			allpsidicts[cellid][geneid] = psi

	#Turn nested dict into df
	df = pd.concat({k: pd.DataFrame.from_dict(v, 'index') for k, v in allpsidicts.items()}, axis = 1)
	df = df.transpose()
	#drop first index and get second (cellid) as column name
	df.reset_index(level = 1, inplace = True, drop = True)
	df.index.name = 'cellid'
	df.reset_index(level = 0, inplace = True)

	#Merge condition column back in
	df = pd.merge(conditiondf, df, how = 'inner', on = 'cellid')

	return df

def calculatedpsi_cellbycell(psidf, conditionA, conditionB):
	#Given a dataframe a cell-level psi values (from calculatepsi_cellbycell())
	#calculate deltapsi (conditionB - conditionA)
	#This is defined as the difference in median psi values between conditions
	#Also report the number of cells in each condition that had psi values and went
	#into the calculation.

	dpsis = {} #{gene : [dpsi, number of cells with psi for condA, number of cells with psi for condB]}

	colnames = list(psidf.columns.values)
	genes = [c for c in colnames if c != 'cellid' and c != 'condition']
	condAdf = psidf.loc[psidf['condition'] == conditionA]
	condBdf = psidf.loc[psidf['condition'] == conditionB]

	for gene in genes:
		condApsis = condAdf[gene].tolist()
		condBpsis = condBdf[gene].tolist()
		condApsis = [psi for psi in condApsis if psi != 'NA']
		condBpsis = [psi for psi in condBpsis if psi != 'NA']
		condAcells = len(condApsis)
		condBcells = len(condBpsis)

		#if either all condA psis are NA or all condB psis are NA, delta psi = NA
		if condBpsis and not condApsis:
			condAmed = 'NA'
			condBmed = round(np.median(condBpsis), 3)
			dpsi = 'NA'
		elif condApsis and not condBpsis:
			condAmed = round(np.median(condApsis), 3)
			condBmed = 'NA'
			dpsi = 'NA'
		elif not condApsis and not condBpsis:
			condAmed = 'NA'
			condBmed = 'NA'
			dpsi = 'NA'
		elif condApsis and condBpsis:
			condAmed = round(np.median(condApsis), 3)
			condBmed = round(np.median(condBpsis), 3)
			dpsi = round(condBmed - condAmed, 3)

		dpsis[gene] = [condAmed, condBmed, dpsi, condAcells, condBcells]

	#Turn into df
	df = pd.DataFrame.from_dict(dpsis, orient = 'index')
	df = df.reset_index()
	df.columns = ['gene', '{0}_psi'.format(conditionA), '{0}_psi'.format(conditionB), 'deltapsi', '{0}_cells'.format(conditionA), '{0}_cells'.format(conditionB)]
	
	return df

def dotest_cellbycell(psidf):
	#Given a matrix of psi values for every gene in every cell (psidf)
	#in which the cells have been split into 2 conditions (column condition)
	#perform a statistical test	to ask whether the psi values for a gene are 
	#different across conditons. This is a bit tricky because most psi values
	#will be NA, due to the fact that most genes have 0 counts in any given cell

	pvals = {} #{gene : pval}

	conditions = psidf['condition'].tolist()
	conditions = list(set(conditions))
	conditionA = conditions[0]
	conditionB = conditions[1]
	condAdf = psidf.loc[psidf['condition'] == conditionA]
	condBdf = psidf.loc[psidf['condition'] == conditionB]

	colnames = list(psidf.columns.values)
	genes = [c for c in colnames if c != 'cellid' and c != 'condition']

	for gene in genes:
		condApsis = condAdf[gene].tolist()
		condApsis = [p for p in condApsis if p != 'NA']
		condBpsis = condBdf[gene].tolist()
		condBpsis = [p for p in condBpsis if p != 'NA']

		#wilcoxon rank-sum test
		#If either condApsis or condBpsis is empty, p = NA
		if not condApsis or not condBpsis:
			p = 'NA'
		else:
			p = ranksums(condApsis, condBpsis)[1]
			p = round(float('{:.2e}'.format(p)), 3)
		pvals[gene] = p

	#Correct pvalues using BH method, but only using pvalues that are not NA
	genes = list(pvals.keys())
	pvalues = list(pvals.values())
	pvaluesformultitest = [pval for pval in pvalues if str(pval) != 'NA']
	fdrs = list(multipletests(pvaluesformultitest, method = 'fdr_bh')[1])
	fdrs = [float('{:.2e}'.format(fdr)) for fdr in fdrs]

	#Now we need to incorporate the places where p = NA into the list of FDRs (also as NA)
	fdrswithnas = []
	fdrindex = 0
	for pvalue in pvalues:
		if str(pvalue) != 'NA':
			fdrswithnas.append(fdrs[fdrindex])
			fdrindex +=1
		elif str(pvalue) == 'NA':
			fdrswithnas.append(np.nan)
	fdrswithnas = [[fdr] for fdr in fdrswithnas] #fdr must be in list for conversion to df

	correctedpvals = dict(zip(genes, fdrswithnas))

	#Make df
	pvalues = [[p] for p in pvalues] #p must be in list for conversion to df
	pvals = dict(zip(genes, pvalues))
	rawpdf = pd.DataFrame.from_dict(pvals, orient = 'index')
	rawpdf.columns = ['raw_pval']
	fdrdf = pd.DataFrame.from_dict(correctedpvals, orient = 'index')
	fdrdf.columns = ['fdr']

	#Merge rawpdf and fdrdf
	df = rawpdf.merge(fdrdf, how = 'inner', left_index = True, right_index = True)
	df = df.reset_index()
	df.columns = ['gene', 'raw_pval', 'FDR']

	return df



		

	


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

    


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--mode', type = str, choices = ['cellbycell', 'subsampleClusters'], help = 'How to perform tests? Either compare psi values of individual cells or subsample cells from clusters.')
	parser.add_argument('--gff', type = str, help = 'GFF of transcript annotation.')
	parser.add_argument('--alevindir', type = str, help = 'Directory containing subdirectories of alevin output.')
	parser.add_argument('--conditions', type = str, help = 'Two column, tab-delimited file with column names \'sample\' and \'condition\'. First column contains cell ids and second column contains cell condition or cluster.')
	parser.add_argument('--conditionA', type = str, help = 'Must be found in the second column of the conditions file. Delta psi is calculated as conditionB - conditionA.')
	parser.add_argument('--conditionB', type = str, help = 'Must be found in the second column of the conditions file. Delta psi is calculated as conditionB - conditionA.')
	args = parser.parse_args()

	#bigdf is a df (obv) of counts with transcripts as columns and cell ids as rows
	print('Collecting count data...')
	bigdf = createbigdf(args.alevindir)
		
	print('Adding condition/cluster info...')
	bigdf = addclusters(bigdf, args.conditions)
	
	print('Assigning transcripts to polyA sites...')
	pklfile = os.path.abspath(args.gff) + '.posfactors.pkl'
	#Check to see if there is a pickled file a position factors
	if os.path.exists(pklfile):
		with open(pklfile, 'rb') as infh:
			posfactors = pickle.load(infh)
	else:
		posfactors = getpositionfactors(args.gff, 25)
		with open(pklfile, 'wb') as outfh:
			pickle.dump(posfactors, outfh)

	if args.mode == 'cellbycell':
		print('Calculating psi values for each cell...')
		psidf = calculatepsi_cellbycell(bigdf, posfactors)
		print('Writing table of psi values...')
		psidf.to_csv('psis.cellbycell.txt.gz', sep = '\t', header = True, index = False, compression = 'gzip')
		print('Calculating delta psi values...')
		deltapsidf = calculatedpsi_cellbycell(psidf, args.conditionA, args.conditionB)
		print('Performing statistical tests...')
		fdrdf = dotest_cellbycell(psidf)
		#Merge deltapsidf and fdrdf
		df = deltapsidf.merge(fdrdf, how = 'inner', on = 'gene', left_index = False, right_index = False)
		#Write df
		df.to_csv('results.cellbycell.txt', sep = '\t', header = True, index = False)
		print('Done!')

	elif args.mode == 'subsampleClusters':
		print('Collapsing counts across all cells within a cluster...')
		countsumdf = collapseclusters(bigdf)
		print('Calculating psi values...')
		psidf = calculatepsi_fromcollapsedclusters(countsumdf, posfactors)
		print('Calculating delta psi values...')
		deltapsidf = calculatedpsi_fromcollapsedclusters(psidf, args.conditionA, args.conditionB)
		print('Performing statistical tests...')
		fdrdf = dotest_bootstrapclusters(bigdf, posfactors, 0.4, 5, args.conditionA, args.conditionB)
		#Merge deltapsidf and fdrdf
		df = deltapsidf.merge(fdrdf, how = 'inner', on = 'gene', left_index = False, right_index = False)
		#Write df
		df.to_csv('results.subsampleclusters.txt', sep = '\t', header = True, index = False)
		print('Done!')
