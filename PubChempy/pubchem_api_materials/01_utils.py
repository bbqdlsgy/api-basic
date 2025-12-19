from collections import defaultdict
import pandas as pd 
import os
import subprocess
from tqdm import tqdm
import re

import requests
import urllib
import urllib.request
from urllib.error import HTTPError, URLError

from rdkit import Chem
from rdkit.Chem import PandasTools

curr_dir = os.path.dirname(os.path.abspath(__file__))

def get_digits(str1):
    c = ""
    for i in str1:
        if i.isdigit():
            c += i
    return c

def get_3d_conformer (file: str):
    """
        get 3d conformer from pubchem using pubchempy 
    Args:
        file (str): file (.csv | .tsv | .xlsx)
    """

    file_path = os.path.join(curr_dir, file)
    
    #-- output directory
    output = os.path.join('/'.join(file_path.split('/')[:-1]), "Patent_US10968172_3D_PubChem_jw")
    os.makedirs(output)
        
    if os.path.isfile(file_path):
        if file_path.split('.')[-1] == 'csv':
            df = pd.read_csv(file_path)

        if file_path.split('.')[-1] == 'tsv':
            df = pd.read_csv(file_path, sep='\t')
                   
        if file_path.split('.')[-1] == 'xlsx':
            df = pd.read_excel(file_path)
        
    for i in df.index:
        
        ID = df.loc[i]['ID']
        cid = df.loc[i]['CID']
        
        smi = df.loc[i]['SMILES']
        
        try:
            url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smi}/SDF?record_type=3d&response_type=save&response_basename=Conformer3D_COMPOUND_CID_{cid}'
            response = urllib.request.urlopen(url)
            savename = f"{output}/{ID}_{cid}_jw.sdf"
            
            # url이 가리키는 주소에 접근해서 해당 자원을 로컬 컴퓨터에 저장하기
            urllib.request.urlretrieve(url, savename)
        except HTTPError as herr:
            print(f'{ID} ' + 'Reason: ', herr.reason)
            url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d&response_type=save&response_basename=Conformer3D_COMPOUND_CID_{cid}'
            response = urllib.request.urlopen(url)
            savename = f"{output}/{ID}_{cid}_jw.sdf"
            
            # url이 가리키는 주소에 접근해서 해당 자원을 로컬 컴퓨터에 저장하기
            urllib.request.urlretrieve(url, savename)

def get_3d_conformer_all (file: str):

    file_path = os.path.join(curr_dir, file)

    #-- parent directory
    parent_dir = os.path.join('/'.join(file_path.split('/')[:-1]), "Patent_US10968172_3D_PubChem_All_jw")
    os.makedirs(parent_dir)
    
    if os.path.isfile(file_path):
        if file_path.split('.')[-1] == 'csv':
            df = pd.read_csv(file_path)

        if file_path.split('.')[-1] == 'tsv':
            df = pd.read_csv(file_path, sep='\t')
                   
        if file_path.split('.')[-1] == 'xlsx':
            df = pd.read_excel(file_path)
        
    for i in df.index:
        
        ID = df.loc[i]['ID']
        cid = df.loc[i]['CID']
        smi = df.loc[i]['SMILES']

        #-- output directory
        output = os.path.join(parent_dir, ID)
        os.makedirs(output)
        
        try:
            id_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/conformers/TXT?record_type=3d"
            response = urllib.request.urlopen(id_url)
            saveid = f"{output}/{ID}_3d_id.txt"
            urllib.request.urlretrieve(id_url, saveid)
            
            response = requests.get(id_url)
            for line in response.text.splitlines():
                print(f"{line}_{cid}")
                id_cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/conformers/{line}/SDF?response_type=save&response_basename=Conformer3D_COMPOUND_CID_{cid}_{line}"
                response = urllib.request.urlopen(id_cid_url)
                saveid = f"{output}/{ID}_{cid}_{line}.sdf"
                urllib.request.urlretrieve(id_cid_url, saveid)
        
        except Exception as err:
            print(f'{saveid} Reason: {err}')

def single_3d_conformer_all (cid):

    #-- output directory
    output = os.path.join(curr_dir, str(cid))
    os.makedirs(output)
    
    try:
        id_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/conformers/TXT?record_type=3d"
        response = urllib.request.urlopen(id_url)
        saveid = f"{output}/CID{cid}_3d_id.txt"
        urllib.request.urlretrieve(id_url, saveid)
        
        response = requests.get(id_url)
        for line in response.text.splitlines():
            print(f"{line}_{cid}")
            id_cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/conformers/{line}/SDF?response_type=save&response_basename=Conformer3D_COMPOUND_CID_{cid}_{line}"
            response = urllib.request.urlopen(id_cid_url)
            saveid = f"{output}/CID{cid}_{line}.sdf"
            urllib.request.urlretrieve(id_cid_url, saveid)
    
    except Exception as err:
        print(f'{saveid} Reason: {err}')
                    
if __name__ == "__main__":
    
    # get_3d_conformer(file='Patent_US10968172_cid_contained.csv')
    # get_3d_conformer_all(file='Patent_US10968172_cid_contained.csv')
    single_3d_conformer_all(cid=171916831)
    # pass