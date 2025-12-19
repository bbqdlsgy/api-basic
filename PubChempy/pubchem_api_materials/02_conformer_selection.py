from __future__ import annotations
import os
from natsort import natsorted #-- sorts an array by using a "natural order" algorithm

import pandas as pd 

from rdkit import Chem
from rdkit.Chem import PandasTools

curr_dir = os.path.dirname(os.path.abspath(__file__))

#-- Select conformer from downloaded 3D conformer (if use get_3d_conformer function)
def extract_conformer_features (conformer_dir: str, output: bool = False):

    parent_dir = os.path.join(curr_dir, conformer_dir)
    
    PUBCHEM_MMFF94_ENERGY: dict = dict()
    PUBCHEM_PHARMACOPHORE_FEATURES: dict = dict()
    
    for _file in natsorted(os.listdir(parent_dir)):
  
        child_file = os.path.join(parent_dir, _file)
        
        if child_file.split('.')[-1] == 'sdf':    
            
            #-- PUBCHEM_MMFF94_ENERGY dictionary key materials
            _fn = child_file.split('.')[0].split('/')[-1]
            
            sdf: pd.DataFrame = PandasTools.LoadSDF(child_file)
            # print(sdf['PUBCHEM_MMFF94_ENERGY'])
            # print(sdf.columns)
            
            #-- PUBCHEM_MMFF94_ENERGY
            if _fn not in PUBCHEM_MMFF94_ENERGY.keys():                
                PUBCHEM_MMFF94_ENERGY[_fn]=sdf.loc[0]['PUBCHEM_MMFF94_ENERGY']
                PUBCHEM_PHARMACOPHORE_FEATURES[_fn]=sdf.loc[0]['PUBCHEM_PHARMACOPHORE_FEATURES']

        
    if PUBCHEM_MMFF94_ENERGY.keys() == PUBCHEM_PHARMACOPHORE_FEATURES.keys():
        #-- 데이터프레임에 넣을 데이터 준비 및 만들기        
        data = {
            'Conformer_ID': [k for k in PUBCHEM_MMFF94_ENERGY.keys()],
            'MMFF94_Energy': [float(v) for k,v in PUBCHEM_MMFF94_ENERGY.items()],
            'Ph4_Features_cnt': [v.split('\n')[0] for k,v in PUBCHEM_PHARMACOPHORE_FEATURES.items()],
            'Ph4_Features': [v.split('\n')[1:] for k,v in PUBCHEM_PHARMACOPHORE_FEATURES.items()]
        }
        df = pd.DataFrame(data)
        
        # #-- 원래 index 저장
        _df_index = df.index
        print(_df_index)
        
        #-- 원래 index 삭제 
        final_df = df.reset_index(drop=True)
        #-- 저장했던 원래 index를 가장 앞에 삽입
        final_df.insert(0, '1st_index', _df_index)
        
        #-- MMFF94_Energy 를 오름차순 정렬
        final_df = final_df.sort_values(['MMFF94_Energy'], ascending=True)
        
        #-- 데이터프레임 저장
        if output:
            final_df.to_csv(f"{parent_dir}/final_total.csv", index=False)
        
        return final_df
    
    else:
        
        return None
    
#-- Select conformer from downloaded 3D conformers (if use get_3d_conformer_all function)
def extract_conformer_all_features (conformer_all_dir: str, output: bool = False):
    """
    Description: Select conformer from downloaded 3D conformers (if use get_3d_conformer_all function)
    """
    parent_dir = os.path.join(curr_dir, conformer_all_dir)
    
    PUBCHEM_MMFF94_ENERGY: dict = dict()
    PUBCHEM_PHARMACOPHORE_FEATURES: dict = dict()
    
    for _dir in natsorted(os.listdir(parent_dir)):
        
        child_dir = os.path.join(parent_dir, _dir)
        
        for _file in natsorted(os.listdir(child_dir)):
            
            child_file = os.path.join(child_dir, _file)
            
            if child_file.split('.')[-1] == 'sdf':    
                
                #-- PUBCHEM_MMFF94_ENERGY dictionary key materials
                _fn = child_file.split('.')[0].split('/')[-1]
                
                sdf = PandasTools.LoadSDF(child_file)
                # print(sdf['PUBCHEM_MMFF94_ENERGY'])
                # print(sdf.columns)
                
                #-- PUBCHEM_MMFF94_ENERGY
                if _fn not in PUBCHEM_MMFF94_ENERGY.keys():                
                    PUBCHEM_MMFF94_ENERGY[_fn]=sdf.loc[0]['PUBCHEM_MMFF94_ENERGY']
                    PUBCHEM_PHARMACOPHORE_FEATURES[_fn]=sdf.loc[0]['PUBCHEM_PHARMACOPHORE_FEATURES']

        
    if PUBCHEM_MMFF94_ENERGY.keys() == PUBCHEM_PHARMACOPHORE_FEATURES.keys():
        #-- 데이터프레임에 넣을 데이터 준비 및 만들기        
        data = {
            'Conformer_ID': [k for k in PUBCHEM_MMFF94_ENERGY.keys()],
            'MMFF94_Energy': [float(v) for k,v in PUBCHEM_MMFF94_ENERGY.items()],
            'Ph4_Features_cnt': [v.split('\n')[0] for k,v in PUBCHEM_PHARMACOPHORE_FEATURES.items()],
            'Ph4_Features': [v.split('\n')[1:] for k,v in PUBCHEM_PHARMACOPHORE_FEATURES.items()]
        }
        df = pd.DataFrame(data)
        
        # #-- 원래 index 저장
        _df_index = df.index
        print(_df_index)
        
        #-- 원래 index 삭제 
        final_df = df.reset_index(drop=True)
        #-- 저장했던 원래 index를 가장 앞에 삽입
        final_df.insert(0, '1st_index', _df_index)
        
        #-- MMFF94_Energy 를 오름차순 정렬
        final_df = final_df.sort_values(['MMFF94_Energy'], ascending=True)
        
        #-- 데이터프레임 저장
        if output:
            final_df.to_csv(f"{parent_dir}/final_total.csv", index=False)
        
        return final_df
    
    else:
        
        return None
                
            
        
    # sdf = PandasTools.LoadSDF()

if __name__ == "__main__":
    
    final_df = extract_conformer_features(conformer_dir = '171916831', output=True)
    print(final_df)

    # final_df = extract_conformer_all_features(conformer_all_dir = 'Patent_US10968172_3D_PubChem_All_jw', output=True)
    # print(final_df)
