import os
import anndata
import pandas as pd
import scanpy as sc
import numpy as np
from .SCSA import Process,Annotator


class pySCSA(object):

    def __init__(self,adata:anndata.AnnData,
                foldchange:float=1.5,pvalue:float=0.05,
                output1:str='temp/scsa1.csv',
                output2:str='temp/scsa2.txt',
                model_path:str='',
                outfmt:str='txt',Gensymbol:bool=True,
                species:str='Human',weight:int=100,tissue:str='All',target:str='cellmarker',
                celltype:str='normal',norefdb:bool=False,cellrange:str=None,
                noprint:bool=True,list_tissue:bool=False) -> None:

        r"""Initialize the pySCSA class

        Arguments:
            adata: AnnData object of scRNA-seq after preprocessing
            foldchange: Fold change threshold for marker filtering. (2.0)
            pvalue: P-value threshold for marker filtering. (0.05)
            output: Output file for marker annotation.(temp/rna_anno.txt)
            model_path: Path to the Database for annotation. If not provided, the model will be downloaded from the internet.
            outfmt: Output format for marker annotation. (txt)
            Gensymbol: Using gene symbol ID instead of ensembl ID in input file for calculation.
            species: Species for annotation. Only used for cellmarker database. ('Human',['Mouse'])
            weight: Weight threshold for marker filtering from cellranger v1.0 results. (100)
            tissue: Tissue for annotation. you can use `get_model_tissue` to see the available tissues. ('All')
            target: Target to annotation class in Database. (cellmarker,[cancersea,panglaodb])
            celltype: Cell type for annotation. (normal,[cancer])
            norefdb: Only using user-defined marker database for annotation.
            noprint: Do not print any detail results.
            list_tissue: List all available tissues in the database.
            cellrange: Cell sub_type for annotation. (if you input T cell, it will only provide T helper cell, T cytotoxic cell, T regulatory cell, etc.)
        
        """

        #create temp directory
        try:
            if not os.path.isdir('temp'):
                print("...Creating directory {}".format('temp'))
                os.makedirs('temp', exist_ok=True)
        except OSError as e:
            print("...Unable to create directory {}. Reason {}".format('temp',e))

        self.adata=adata
        self.foldchange=foldchange
        self.pvalue=pvalue
        self.output1=output1
        self.output2=output2
        self.outfmt=outfmt
        self.Gensymbol=Gensymbol
        self.species=species
        self.weight=weight
        self.tissue=tissue
        self.celltype=celltype
        self.norefdb=norefdb
        self.noprint=noprint
        self.list_tissue=list_tissue
        self.target=target
        self.cellrange=cellrange
        if model_path =='':
            self.model_path=data_downloader(url='https://figshare.com/ndownloader/files/41369037',
                                            path='temp/pySCSA_2023_v2_plus.db',title='whole')
        else:
            self.model_path=model_path

    def cell_anno(self,clustertype:str='leiden',
                  cluster:str='all',rank_rep=False)->pd.DataFrame:
        r"""Annotate cell type for each cluster.
        
        Arguments:
            clustertype: Clustering name used in scanpy. (leiden)
            cluster: Only deal with one cluster of marker genes. (all,[1],[1,2,3],[...])
        """

        data_preprocess(self.adata,clustertype=clustertype,path=self.output1,rank_rep=rank_rep)

        # print('...Auto annotate cell')
        
        p = Process()
        p.run_cmd_p(foldchange=self.foldchange,
                    weight=self.weight,
                    pvalue=self.pvalue,
                    tissue=self.tissue,
                    species=self.species,
                    target=self.target,
                    norefdb=self.norefdb,
                    MarkerDB=None,
                    db=self.model_path,
                    noprint=self.noprint,
                    input=self.output1,
                    output=self.output2,
                    source="scanpy",
                    cluster=cluster,
                    fc=self.foldchange,
                    outfmt=self.outfmt,
                    celltype=self.celltype,
                    Gensymbol=self.Gensymbol,
                    list_tissue=self.list_tissue,
                    cellrange=self.cellrange)
        

        result=pd.read_csv(self.output2, sep='\t')
        self.result=result
        return result


def data_preprocess(adata, path, clustertype='leiden',
                    rank_rep=False):
    r"""data preprocess for SCSA
    
    Parameters
    ----------
    - adata: `AnnData`
        AnnData object
    - path: `str`   
        the save path of datasets

    Returns
    -------
    - adata: `AnnData`
        AnnData object
    """
    dirname, _ = os.path.split(path)
    try:
        if not os.path.isdir(dirname):
            print("......Creating directory {}".format(dirname))
            os.makedirs(dirname, exist_ok=True)
    except OSError as e:
        print("......Unable to create directory {}. Reason {}".format(dirname,e))

    sc.settings.verbosity = 2  # reduce the verbosity
    if rank_rep==False and 'rank_genes_groups' not in adata.uns.keys():
        sc.tl.rank_genes_groups(adata, clustertype, method='wilcoxon')
    elif rank_rep==True:
        sc.tl.rank_genes_groups(adata, clustertype, method='wilcoxon')
    else:
        pass
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    dat = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
    dat.to_csv(path)
    return dat

def fullAnno(anno):
    pd_li = []
    for i in set(anno['Cluster']):
        anno1 = anno.loc[anno['Cluster']==i].iloc[:2].reset_index(drop=True)
        pd1 = pd.DataFrame(
        {
            "Cluster": [anno1.loc[0, 'Cluster']],
            "CellType1": [anno1.loc[0, 'Cell Type']],
            "Z-score1": [anno1.loc[0, 'Z-score']],
            "CellType2": [anno1.loc[1, 'Cell Type']],
            "Z-score2": [anno1.loc[1, 'Z-score']],
            "Quality": "good" if anno1.loc[0, 'Z-score']>anno1.loc[1, 'Z-score']*2 else "bad"        
        })
        pd_li.append(pd1)
    full_anno = pd.concat(pd_li, axis=0, ignore_index=True)   
    return(full_anno)  

def simpleAnno(anno):
    pd_li = []
    for i in set(anno['Cluster']):
        anno1 = anno.loc[anno['Cluster']==i].iloc[:2].reset_index(drop=True)
        pd1 = pd.DataFrame(
        {
            "Cluster": [anno1.loc[0, 'Cluster']],
            "CellType": [anno1.loc[0, 'Cell Type']]
        })
        pd_li.append(pd1)
    full_anno = pd.concat(pd_li, axis=0, ignore_index=True)   
    return(full_anno)  

def scsaResults(adata, anno, clustertype):
    obs_df = adata.obs.copy()
    obs_df['cellbarcode_scsa'] = obs_df.index
    anno['Cluster'] = anno['Cluster'].astype("str").astype("category")
    obs_full = pd.merge(obs_df, anno, left_on= clustertype, right_on='Cluster', sort=False)
    obs_full.index = obs_full['cellbarcode_scsa']
    obs_full = obs_full.iloc[:,obs_df.shape[1]:]
    obs_full = obs_full.loc[obs_df.index.tolist()]
    return(obs_full)

def run_SCSA(adata, celltype, model_path,
             clustertype, rank_rep, full_anno=True,
             ):
    scsa = pySCSA(adata=adata,
                      foldchange=1.5,
                      pvalue=0.01,
                      celltype=celltype,
                      target='cellmarker',
                      tissue='All',
                      model_path=model_path                    
    )
    all_result = scsa.cell_anno(clustertype=clustertype,
               cluster='all',rank_rep=rank_rep)
    if full_anno:
        anno = fullAnno(all_result)
    else:
        anno = simpleAnno(all_result)
    result = scsaResults(adata, anno, clustertype)
    return result