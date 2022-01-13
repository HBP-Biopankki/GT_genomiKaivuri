# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 13:31:41 2022

@author: HUS45090218
"""
import numpy as np
import pandas as pd
from time import sleep
import gzip
import os, sys
from io import StringIO


###Input metodi TBA
user_input = '''
'''

user_input_DF = pd.read_csv(StringIO(user_input), sep=',')

user_kromosomit = list(set(user_input_DF['CR'].to_list()))
user_paikat     = list(set(user_input_DF['pos'].to_list()))

user_kromosomit = [bytes(str(i),'utf-8') for i in user_kromosomit]


###  Config

VCF_file_type = None
vcf_path = r'O:\Raw_genotypes\test'
temp_path = os.path.join(vcf_path, 'tmp_folder')


### tmp_folderi
if os.path.exists(temp_path): 
    print('Temp folder on olemassa: ', temp_path)
    # exit()
else: 
    os.mkdir(temp_path)
    
    

### Tiedostot
all_files = pd.DataFrame(os.listdir(vcf_path), columns=['file_nimi'])
if all_files['file_nimi'].str.contains('vcf.gz').any(): 
    VCF_file_type = 'vcf.gz'
    all_files = all_files[all_files['file_nimi'].str.contains('vcf.gz')]
elif all_files['file_nimi'].str.contains('vcf').any(): 
    VCF_file_type = 'vcf'
    all_files = all_files[all_files['file_nimi'].str.contains('vcf')]
else: 
    print('VCF tieostojen tunnistuksessa virhe.')
    sleep(5)
    exit() #





###Suurempi .vcf.gz prosssointi


for index, read_file in enumerate(all_files['file_nimi']):     
    print('Working on file ', read_file)
    output = gzip.open(os.path.join(temp_path, 'small_'+read_file) , 'wb') 

    with gzip.open(os.path.join(vcf_path, read_file), 'rb') as f:
        for idx, line in enumerate(f):
            if line[:1] == b'#': 
                output.write(line)
            else:
                tmp_palat = line.split(b'\t')
                if tmp_palat[0] in(user_kromosomit):
                    if min([abs(i - int(tmp_palat[1].decode()) ) for i in user_paikat]) <= 10:
                        # print('match')
                        output.write(line)
    output.close()
                    
###Suurempi .vcf prosssointi    TBA
#Tarvitaanko tätä oikeasti? 
    
                
###Yhdistetaan kaikkki filet. 
def get_header(tiedosto):
    vcf_header = []
    with gzip.open(os.path.join(temp_path, tiedosto), "rt") as file_open:
          for line in file_open:
            if line.startswith("#CHROM"):
                  vcf_header = [x for x in line.split('\t')]
                  break
    file_open.close()
    return vcf_header


yhteiset = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
master_tulos = pd.DataFrame(columns= yhteiset )
small_files = pd.DataFrame(os.listdir(temp_path), columns=['file_nimi'])


for file in small_files['file_nimi']: 
    print(file)
    header = get_header(file)
    print(header)
    try: #quick and dirty... 
        tmp_data = pd.read_csv(os.path.join(temp_path, file), comment='#', sep='\t', compression='gzip', header=None, names=header)
        master_tulos = pd.merge(master_tulos, tmp_data, how = 'outer', on = yhteiset ) #DEV_FLAG: Vaikuttaa toimivan, tuplaTarkista
    except: 
        pass

    

###Siivotaan 0/0 tai NaN pois... 
poistoon = []
for column in master_tulos.columns: 
    if column not in yhteiset:
        # print(master_tulos[column].unique())
        if all([i in [np.nan,'./.' , '0/0']  for i in master_tulos[column].unique().tolist()]): 
            poistoon.append(column)

master_tulos = master_tulos.drop(columns = poistoon)
    
###Talletetaan 





########
def dev_stuff(): 

    def get_metaFiles() -> pd.DataFrame:
        '''hakee  '''
        chip_mk2 = pd.read_csv(r'C:\Downloads\Axiom_stuff\Axiom_Finngen2.na36.r4.a2.annot.csv', sep=';')
        chip_mk1 = pd.read_csv(r'C:\Downloads\Axiom_stuff/finngen1_variants_rsid.txt', sep='\t')
        return chip_mk2, chip_mk1
    
    def foo(): 
        kohde = 'C:\Downloads\Axiom_stuff\Axiom_Finngen2.na36.r4.a2.annot.csv'
        my_data = pd.read_csv(kohde,skiprows = 20, sep=',' )
        my_data.columns

    
    def crete_config()-> None: 
        pass
    
    def main_options(): #?
        pass



    
    

if __name__ == '__main__': 
    main()

