import warnings

import urllib.request
import time
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool

import pandas as pd
import matplotlib.pyplot as plt
import numpy

import GEOparse

def remove_repeated_var_inds(adata_tmp):
    i_dels = [True]* adata_tmp.n_vars
    adata_tmp.var['i'] = list(range(0, adata_tmp.n_vars))
    reps  = adata_tmp.var.index.value_counts()
    rep_names = list(reps[reps > 1].index)

    for rep in rep_names:
        i_del_ = adata_tmp.var.loc[ rep, 'i']
        for i_del in i_del_:
            i_dels[i_del] = False
    
    adata_tmp.var.drop(columns=['i'], inplace=True)
    return adata_tmp[:, i_dels]


def join_map_mart(adata_tmp, annot, gene_annot='external_gene_name', how='left'):
    '''
    Takes the index of the first arg AnnData and reduces the df
    in the right to those indexes
    '''
    maps_gene_names = pd.merge(
        pd.DataFrame({gene_annot:adata_tmp.var.index}),
        annot,
        how=how, on=gene_annot, suffixes=('_xxx','_yyy'))
    maps_gene_names.rename(columns={'external_gene_name':'gene_symbol'}, inplace=True)
    maps_gene_names.drop_duplicates(subset=['gene_symbol'], inplace=True)
    maps_gene_names.set_index('gene_symbol', inplace=True)
    return maps_gene_names


def search_str_index(s_list:list, search_str:str) -> list:
    """
    Get the indexes for in a list that contain a str
    
    Parameters:
        s_list (list): the list to search in
        search_str (str): the str to search for, it only must be contained

    Returns:
        matched_indexes (list): A list of the indexes where it was found
    """
    
    matched_indexes = []
    i = 0
    length = len(s_list)

    while i < length:
        if type(s_list[i]) is not str:
            i += 1
            continue
        if search_str in s_list[i]:
            matched_indexes.append(i)
        i += 1
        
    return matched_indexes
    


def get_geo_exprs(gse_str='', data_dir='/root/datos/maestria/netopaas/lung_scRNA'):
    """
    Builds a metada matrix from geo SOFT file for scRNA-seq data that
    doesn't have the expression matrix attached. A csv is saved in
    f'{data_dir}/{gse_str}/{gse_str}_metadata.csv'
    
    Parameters:
        gse_str (str): The string of the GSE to be gotten
        data_dir (str): Where the CSV is to be saved

    Returns:
        metadata (dict): A dict with first level the features, second the sample
    
    """
    gse = GEOparse.get_GEO(geo=gse_str, destdir=f"{data_dir}/{gse_str}/", silent=True)

    # Get expression data and metadata matrices
    exprs = []
    gsmNames = []
    metadata = {}
    sup_dict = {}
    for gsm_name, gsm in gse.gsms.items():
        if gsm.metadata['type'][0]=='SRA':
             # Expression data
            # print(gsm.__dict__)
            if len(gsm.table)>0:
                # TODO will there really be a table here anytime?
                # if so run code here
                pass
            else:
                # Get the supplementary files with their type because no table is attached
                sup_file_url = gsm.metadata['supplementary_file_1'][0]
                
                # TODO it is in this array but no standard index, search for it later
                l1 = gsm.metadata['data_processing']
                s = 'upplementary'
                matched_indexes = search_str_index(l1, s)
                sup_file_type = l1[matched_indexes[0]].split(':')[1]
                
                l1 = gsm.metadata['data_processing']
                s = 'Genome_build'
                matched_indexes = search_str_index(l1, s)
                genome_build = l1[matched_indexes[0]].split(':')[1]
                
                if not 'sup_file_url' in sup_dict:
                    sup_dict['sup_file_url'] = {}
                sup_dict['sup_file_url'][gsm_name] = sup_file_url
                
                if not 'sup_file_type' in sup_dict:
                    sup_dict['sup_file_type'] = {}
                sup_dict['sup_file_type'][gsm_name] = sup_file_type
                
                if not 'genome_build' in sup_dict:
                    sup_dict['genome_build'] = {}
                sup_dict['genome_build'][gsm_name] = genome_build
                
                
                # print('No expression table, saving supplementary file url'
                #              f'{sup_file_url} with type: {sup_file_type}')
            if hasattr(gsm, 'SRA'):
                warnings.warn("There is an SRArun access table, consider using "
                              "your snakemake workflow to parallely download them")
                
    # Metadata
            for key,value in gsm.metadata.items():
                if (key=='characteristics_ch1' or key=='characteristics_ch2') and (len([i for i in value if i!=''])>1 or value[0].find(': ')!=-1):
                    tmpVal = 0
                    for tmp in value:
                        splitUp = [i.strip() for i in tmp.split(':')]
                        if len(splitUp)==2:
                            if not splitUp[0] in metadata:
                                metadata[splitUp[0]] = {}
                            metadata[splitUp[0]][gsm_name] = splitUp[1]
                        else:
                            if not key in metadata:
                                metadata[key] = {}
                            metadata[key][gsm_name] = splitUp[0]
                else:
                    if not key in metadata:
                        metadata[key] = {}
                    if len(value)==1:
                        metadata[key][gsm_name] = ' '.join([j.replace(',',' ') for j in value])
            ftp_name = sup_dict['sup_file_url'][gsm_name].split('/')[-1]
            
            key = 'sup_type'
            if not key in metadata:
                metadata[key] = {}
            metadata[key][gsm_name] = sup_dict['sup_file_type'][gsm_name]
            
            key = 'local_path'
            if not key in metadata:
                metadata[key] = {}
            metadata[key][gsm_name] = f'{data_dir}/{gse_str}/{ftp_name}'
            
            key = 'genome_build'
            if not key in metadata:
                metadata[key] = {}
            metadata[key][gsm_name] = sup_dict['genome_build'][gsm_name]
            
            metadata['local_path'][gsm_name] = f'{data_dir}/{gse_str}/{ftp_name}'
    pd.DataFrame(metadata).to_csv(f'{data_dir}/{gse_str}/{gse_str}_metadata.csv')

    return metadata

def download_url(args, gse='', data_dir='/root/datos/maestria/netopaas/lung_scRNA'):
    t0 = time.time()
    #Extract url and file path from args
    url = args[0]
    path = args[1]
    try:
        # For getting ftp that does not need 
        gsm_path, response = urllib.request.urlretrieve(url,
                                      path)
        return(url, time.time() - t0)
    except Exception as e:
        print('Exception in download_url():', e)


def download_parallel(args, gse=''):
    cpus = int(cpu_count()/3)
    print("CPUS: ", cpus)
    
    # Fix the download func for a specific GSE so that we dnot have o make a list of the same GSE
    download_url_gse = lambda args : download_url(args, gse=gse)
    with ThreadPool(cpus -1 ) as pool:
        for result in pool.imap_unordered(download_url_gse, args):
            print('url:', result[0], 'time (s):', result[1])
