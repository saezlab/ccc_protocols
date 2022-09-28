import os
import subprocess

import pandas as pd
import scanpy as sc

class CovidBalf():
    """Load the BALF covid data set from Liao 2020 (https://doi.org/10.1038/s41591-020-0901-9)"""
    def __init__(self, data_path: str = './'):
        """Init method

        Parameters
        ----------
        self.data_path : str
            Location to store balf files locally formatted as "full/path/to/balf_data_dir/"
        """
        self.data_path = data_path

    def download_data(self):
        """Downloads the metadata and expression data."""
        # download the metadata
        metadata_link = 'https://raw.githubusercontent.com/zhangzlab/covid_balf/master/all.cell.annotation.meta.txt'
        cmd = 'wget ' + metadata_link + ' -O ' + os.path.join(self.data_path, 'metadata.txt')
        subprocess.run(cmd, shell=True, 
                       stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        # download the expression data
        sample_links = [
            'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4339nnn/GSM4339769/suppl/GSM4339769%5FC141%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5',
            'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4339nnn/GSM4339770/suppl/GSM4339770%5FC142%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5',
            'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4339nnn/GSM4339771/suppl/GSM4339771%5FC143%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5', 
            'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4339nnn/GSM4339772/suppl/GSM4339772%5FC144%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5', 
            'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4339nnn/GSM4339773/suppl/GSM4339773%5FC145%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5',
            'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4339nnn/GSM4339774/suppl/GSM4339774%5FC146%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5', 
            'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4475nnn/GSM4475048/suppl/GSM4475048%5FC51%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5', 
            'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4475nnn/GSM4475049/suppl/GSM4475049%5FC52%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5', 
            'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4475nnn/GSM4475050/suppl/GSM4475050%5FC100%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5', 
            'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4475nnn/GSM4475051/suppl/GSM4475051%5FC148%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5', 
            'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4475nnn/GSM4475052/suppl/GSM4475052%5FC149%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5',
            'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4475nnn/GSM4475053/suppl/GSM4475053%5FC152%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5']

        for sl in sample_links:
            cmd = 'wget ' + sl + ' -P ' + os.path.join(self.data_path)
            subprocess.run(cmd, shell=True, 
               stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    def format_data(self):
        """Format the metadata and expression data into Python objects.

        Returns
        -------
        md: pd.DataFrame
            the formatted cell metadata
        balf_samples : Dict[str, AnnData]
            keys are the balf_samples, values are the raw UMI counts stored in an AnnData object
        """
        
        # format the metadata
        md = pd.read_csv(os.path.join(self.data_path, 'metadata.txt'), sep='\t')
        md.rename(columns = {'sample': 'Sample_ID', 'group': 'Context', 'celltype': 'cell_type', 'ID': 'cell_barcode'}, 
                        inplace = True)
        md['Context'] = md.Context.map({'HC': 'Healthy_Control', 'M': 'Moderate_Covid', 'S': 'Severe_Covid'}).tolist()
        md['Context'] = pd.Categorical(md.Context, categories = ['Healthy_Control', 'Moderate_Covid', 'Severe_Covid'])
        md = md.sort_values(by = ['Context', 'Sample_ID'], ascending = True)

        md = md[md.Sample_ID.isin([sample_id for sample_id in md.Sample_ID.unique() if sample_id.startswith('C')])] 

        md.set_index('cell_barcode', inplace = True)

        # format the expression data
        balf_samples = dict()
        for filename in os.listdir(self.data_path):
            if filename.endswith('.h5'):
                sample = filename.split('_')[1]

                # subset and format metadata
                md_sample = md[md.Sample_ID == sample]
                md_sample.index = [cell_barcode.split('_')[0] + '-1' for cell_barcode in md_sample.index]

                adata = sc.read_10x_h5(self.data_path + filename)
                adata = adata[md_sample.index, ] # only include cells present in the metadata
                adata.obs = md_sample[['Sample_ID', 'Context', 'cell_type']]
                adata.var_names_make_unique()
                balf_samples[sample] = adata
        balf_samples = {sample_id: balf_samples[sample_id] for sample_id in md.Sample_ID.unique()} # reorder
        
        for sample, adata in balf_samples.items():
            sc.pp.filter_genes(adata, min_cells=3)
        
        return md, balf_samples
