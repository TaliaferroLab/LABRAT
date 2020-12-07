#Application of LABRAT to alevin-quantified scRNAseq data

from struct import Struct
import numpy as np
import pandas as pd
import gzip
import sys
import os
#import sce
from scipy.io import mmread
from scipy.sparse import csr_matrix
import gffutils
import pickle


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
                print ("Found total " + str(tot_umi_count) + " reads")
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

    return alv.head(n = 100) #only 100 cells for now

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

    

#bigdf = createbigdf(sys.argv[1])
#bigdf.to_csv('bigdf.txt', sep = '\t', header = True, index = True)
bigdf = pd.read_csv('bigdf.txt.gz', sep = '\t', header = 0, index_col = 0, compression = 'gzip')
#posfactors = getpositionfactors(sys.argv[2], 25)
#with open('posfactors.pkl', 'wb') as outfh:
#    pickle.dump(posfactors, outfh)
with open('posfactors.pkl', 'rb') as infh:
   posfactors = pickle.load(infh)

bigdf = addclusters(bigdf, sys.argv[2])
#countsumdf = collapseclusters(bigdf)
#psidf = calculatepsi_fromcollapsedclusters(countsumdf, posfactors)

psidf = calculatepsi_cellbycell(bigdf, posfactors)