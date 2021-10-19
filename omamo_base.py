#!/usr/bin/python3

from sys import argv

import pyoma.browser.db as db

import pandas as pd
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import statistics as stats

import matplotlib
matplotlib.use('Agg')

import pickle


def inf_cont(address):

	'''Create a dictionary {GO term ID: information content}'''

#Create a data frame from the text file, remove the third column
	df = pd.read_csv(address, sep = '\t')
	del df['freq']

#Convert the dataframe to an array
	df = df.to_numpy()

#Create a dictionary where {GO term ID: informaiton content}
	dictionary = {}
	for i in df:
    		dictionary[i[0]] = i[1]
	
	return dictionary



def pairwise_orthologs(species, model_organism):
    
    '''Method that finds pairwise orthologs between two species
    Output: NumPy array, 0th column is the gene of the species of interest and
    1st column is the gene of a potential model organism'''
    
    COLUMN_COUNTER = 2
    ROW_COUNTER = 0
    
    zeroth_column_list = []
    first_column_list = []
    
    model_genome = [i[0] for i in oma_db.all_proteins_of_genome(model_organism)]
    species_genome = [i[0] for i in oma_db.all_proteins_of_genome(species)]
    
    for gene in species_genome:
        genes_orthologs = [pair[1] for pair in oma_db.get_vpairs(gene)]
        
        for ortholog in genes_orthologs:
            if ortholog in model_genome:
                ROW_COUNTER += 1
                zeroth_column_list.append(gene)
                first_column_list.append(ortholog)
                
    output = np.empty((ROW_COUNTER, COLUMN_COUNTER), dtype = int)
    
    output[:,0] = zeroth_column_list
    output[:,1] = first_column_list
    
    return(output)


def GO_overlap(array):
    
    '''This function determines common GO terms.
    Input: array, a numpy array of pairwise orthologs
    Output is a list of lists:
    [[pair, [bio process GO overlap], perc_similarity],...]'''
    
    output = []
    
    for pair in array:
        
        annos_gene1 = oma_db.get_gene_ontology_annotations(pair[0])
        annos_gene2 = oma_db.get_gene_ontology_annotations(pair[1])
        
        term1 = set(oma_db.gene_ontology.term_by_id(annos_gene1[i]['TermNr']) for i in range(len(annos_gene1)))
        term2 = set(oma_db.gene_ontology.term_by_id(annos_gene2[i]['TermNr']) for i in range(len(annos_gene2)))

	#Consider parent terms
        
        term1_plus_parents = set()
        for i in term1:
            term1_plus_parents.update(go.get_superterms_incl_queryterm(i))
        term1_plus_parents_ids = set(t.id for t in term1_plus_parents)
        
        term2_plus_parents = set()
        for i in term2:
            term2_plus_parents.update(go.get_superterms_incl_queryterm(i))
        term2_plus_parents_ids = set(t.id for t in term2_plus_parents)
        
        overlap = [int(('{0}'.format(a)).strip('GO:')) for a in term1_plus_parents_ids.intersection(term2_plus_parents_ids)]
        overlap_defined = []
        
        #Only consider GO terms that are present in the go_annotations file and with information content >= 5
        for i in overlap:
            if i in list(inf_content.keys()):
                if inf_content[i] >=5:
                    overlap_defined.append(i)
                
        union = [int(('{0}'.format(a)).strip('GO:')) for a in term1_plus_parents_ids.union(term2_plus_parents_ids)]
        union_defined = []
        
        for i in union:
            if i in list(inf_content.keys()):
                if inf_content[i] >=5:
                    union_defined.append(i)
        
        #The output only includes BP overlapping GO terms
        bio_process_overlap = []
        for t in overlap_defined:
            if oma_db.gene_ontology.term_by_id(t).aspect == 1:
                bio_process_overlap.append(t)
        
        if len(bio_process_overlap) != 0:
            small_list = [pair]
            small_list.append(bio_process_overlap)
            
            intrsctn_inf_cont = sum(list(inf_content[i] for i in overlap_defined))
            union_inf_cont = sum(list(inf_content[e] for e in union_defined))
                
            try:
                similarity = intrsctn_inf_cont/union_inf_cont
                small_list.append(similarity)
                output.append(small_list) 
            except ZeroDivisionError:
                pass
        
    return output


def high_similarity_fnctn_ogs(GO_overlap_output, threshold):
    
    '''This function finds orthologous pairs where similarity > threshold
    Input: GO_overlap array & lower thresshold
    Output: list of lists [[pair, [bio process GO overlap], perc_similarity],...] '''
    
    output = []
    
    for i in GO_overlap_output:
         if i[2] >= threshold:
                output.append(i)
    
    return output


def frequency_GO_filter(high_sim_fn_output, lower_limit, upper_limit):
    
    '''Gives a list of overlapping GO terms that occur between lower_ and upper_limit times + count
    Input: high functional similarity array & lower and upper threshold count values
    Output: [(GO_ID1, count1), (GO_ID2, count2),...]'''
    
    large_list = []
    output = []
    
    for i in high_sim_fn_output:
        large_list += i[1]
        
    large_set = set(large_list)
    
    for i in large_set:
        count = large_list.count(i)
        if count >= lower_limit and count <= upper_limit:
            output.append((i, count))
            
    return output


def gene_name_mapper(entry_nr):
    
    '''Map OMA ID to UniProt Gene Name'''
	
    xref_mapper = oma_db.id_mapper['XRef']

    x = xref_mapper.map_entry_nr(entry_nr)
    
    for i in x:
        if i['source'] == 'UniProtKB/SwissProt':
            return(i['xref'])
    for i in x:
        if i['source'] == 'UniProtKB/TrEMBL':
            return(i['xref'])
    for i in x:
        if i['source'] == 'Gene Name':
            return(i['xref'])
    for i in x:
        if i['source'] == 'EnsemblGenomes':
            return(i['xref'])

    return(str(entry_nr))


def table2(species, frequency_GO_filter_output,high_similarity_fnctn_ogs_output, threshold):
    
    '''Creates pandas dataframe [GO ID, [human genes], [model genes], no. of orthologs, av. functional similarity & st. dev., sum of similarity values]'''
    
    output = [] 
    
    for i in frequency_GO_filter_output:
        small_list = []
        human_genes = []
        sp_genes = []
        try:
            if inf_content[i[0]] >= threshold:
                small_list.append(i[0])
                small_list.append(species)
                s = []
                for e in high_similarity_fnctn_ogs_output:
                    if i[0] in e[1]:
                        human_genes.append(e[0][0])
                        sp_genes.append(e[0][1])
                        s.append(e[2])
                av_sim = sum(s)/len(human_genes)
                if len(s) > 1:
                    st_dev = round(stats.stdev(s),4)
                else:
                    st_dev = '0.0000'
                av_and_st_dev = str(round(av_sim,4)) + ' ± ' + str(st_dev)
                h_genes_uniprot = []
                for t in human_genes:
                    h_genes_uniprot.append(gene_name_mapper(t))
                small_list.append(','.join(h_genes_uniprot))
                sp_genes_uniprot = []
                for t in sp_genes:
                    sp_genes_uniprot.append(gene_name_mapper(t))
                small_list.append(','.join(sp_genes_uniprot))
                small_list.append(len(human_genes))
                small_list.append(av_and_st_dev)
                small_list.append(round(sum(s), 2))
                output.append(small_list)
        except KeyError:
            pass
    
    df = pd.DataFrame(output, columns = ['GO ID', 'Species', 'Human Genes', 'Species Genes', 'No. of OGs','Average func. similarity ± st. dev', 'Score'])
                
    return df


def table1(species, frequency_GO_filter_output,high_similarity_fnctn_ogs_output, threshold):
    
    '''Create table ['GO ID', 'Species', 'Human Gene', 'Species Gene', 'Functional Similarity']'''
    output = []
    
    for i in high_similarity_fnctn_ogs_output:
        small_list = []
        for e in frequency_GO_filter_output:
            try:
                if inf_content[e[0]] >= 5:
                    if e[0] in i[1]:
                        small_list.append(species)
                        small_list.append(gene_name_mapper(i[0][0]))
                        small_list.append(gene_name_mapper(i[0][1]))
                        small_list.append(round(i[2],4))
                        output.append(small_list)
                        break
            except KeyError:
                pass
        
    
    df = pd.DataFrame(output, columns = ['Species', 'Human Gene', 'Species Gene', 'Func. Similarity' ])
                    
    return df
     

if __name__ == "__main__":

	'''ARG1 - OmaServer Address, ARG2 - go_positive_annos file address, ARG3 - species code'''

	oma_address = argv[1]
	oma_db = db.Database(oma_address)
	go = oma_db.gene_ontology
	
	list_speices = [z.uniprot_species_code for z in oma_db.tax.genomes.values()]
	print("There are",len(list_speices),"species in the OMA database.")

	go_file_address = argv[2]
	inf_content = inf_cont(go_file_address)

	species = argv[3]

	print('Finding orthologs in', species)

	og_pairs = pairwise_orthologs('HUMAN', species)
	print('There are {} pairwise orthologs in {}.\n'.format(len(og_pairs), species))
	
	file = 'OGs_{}.pkl'.format(species)
	open_file = open(file, 'wb')
	pickle.dump(og_pairs, open_file)
	open_file.close()
	print('ORTHOLOGOUS PAIRS FOR {} HAVE BEEN SAVED AS A PICKLE FILE'.format(species))

	overlap = GO_overlap(og_pairs)
	similar_ogs = high_similarity_fnctn_ogs(overlap, 0.05)
	high_fr_go = frequency_GO_filter(similar_ogs,1,5000)

	output1 = table1(species, high_fr_go,similar_ogs,5)
	filename1 = species + '1.csv'
	output1.to_csv(filename1,sep='\t')
	print('Table 1 for {} has been saved.'.format(species))
	output2 = table2(species, high_fr_go,similar_ogs,5)
	filename2 = species + '2.csv'
	output2.to_csv(filename2,sep='\t')
	print('Table 2 for {} has been saved.\n'.format(species))
	print('DATA FOR {} HAS BEEN COLLECTED.\n'.format(species))
	

	print('===DATA FOR SPECIES {} HAS BEEN COLLECTED==='.format(species))

	
