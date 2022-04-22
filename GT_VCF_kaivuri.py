# -*- coding: utf-8 -*-
"""
@author: HUS45090218
"""
import numpy as np
import pandas as pd
from time import sleep
import gzip
import os
from io import StringIO
from datetime import datetime

#Kansio missa raaka-data on. 
vcf_path = r'O:\Raw_genotypes\DF4'

def main(): 
    
    ###Input metodi TBA
    user_input = '''\
CR,pos
13,32330992
13,32331004
13,32332764
13,32338202
13,32338236
13,32338581
13,32339658
13,32339924
13,32339973
13,32340622
13,32340630
13,32356550
13,32363529
13,32363535
13,32370400
13,32371097
13,32379913
13,32380005
16,23634954
17,43045700
17,43045767
17,43047643
17,43057062
17,43057078
17,43057120
17,43067646
17,43067684
17,43071225
17,43074350
17,43090946
17,43091034
17,43091496
17,43091772
17,43091905
17,43092046
17,43092142
17,43092302
17,43092386
17,43092440
17,43092608
17,43092845
17,43092932
17,43093056
17,43093726
17,43093844
17,43093919
17,43094023
17,43094439
17,43094603
17,43094794
17,43095848
17,43095919'''
    
    user_input_DF = pd.read_csv(StringIO(user_input), sep=',')
    user_input_DF
    
    user_input_DF.columns
    user_kromosomit = list(set(user_input_DF['CR'].to_list())) #Dev flag: pitaa hakea tupleina 
    user_paikat     = list(set(user_input_DF['pos'].to_list()))
    
    user_kromosomit = [bytes(str(i),'utf-8') for i in user_kromosomit]
    
    
    ###  Config
    VCF_file_type = None
    # vcf_path = vcf_path # koodin alussa. 
    temp_path = os.path.join(vcf_path, 'tmp_folder')
    
    ### tmp_folderi
    if os.path.exists(temp_path): 
        print('Temp folder on olemassa: ', temp_path)
        # os._exit(1)
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
                        
    ###Suurempi .vcf prosssointi    TBA , ehka
    ###
    ###
    ###
                    
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
        # print(header)
        try: #quick and dirty... 
            tmp_data = pd.read_csv(os.path.join(temp_path, file), comment='#', sep='\t', compression='gzip', header=None, names=header)
            
            # Etsitaan yhteiset columnit
            common = list(set(master_tulos.columns.to_list()).intersection(header) )
            master_tulos = pd.merge(master_tulos, tmp_data, how = 'outer', on = common) 
        except Exception: 
            print(Exception)
            pass
    
    
    ###Siivotaan 0/0 tai NaN pois... 
    poistoArvot = [np.nan, '0/0'] #
    poistoon = []
    for column in master_tulos.columns: 
        if column not in yhteiset:
            if all([i in poistoArvot  for i in master_tulos[column].unique().tolist()]): 
                poistoon.append(column)

    
    master_tulos = master_tulos.drop(columns = poistoon)
        
    ###Talletetaan 
    def talletus():
        master_tulos.to_csv(os.path.join(temp_path, 'master.csv'), index=False)
    def load_prior(): 
        master_tulos = pd.read_csv(os.path.join(temp_path, 'master.csv'))
        return master_tulos 
    
    
    ###frekvenssit    
    
    tmp = master_tulos.T
    tmp.columns = tmp.iloc[2] +'^'+ tmp.iloc[6] 
    tmp= tmp.iloc[9:]

    freqs = pd.DataFrame()
    for col in tmp.columns: 
        freqs = pd.merge(freqs, tmp[col].value_counts().to_frame(), how='outer', left_index = True, right_index=True)

    freqs = freqs.T
    foo = freqs.reset_index()
    bar = foo['index'].str.split('^', expand=True).rename(columns={0: 'Ax_code', 1:'filter'})
    out = pd.merge(bar,foo , how='inner', left_index=True, right_index=True).sort_values(['Ax_code','filter'])
    freqs.to_csv(os.path.join(temp_path, 'frekvenssit.csv'))
    
    out.to_csv(os.path.join(temp_path, 'frekvenssit.csv'), index=False)


    ### Hienompi haku. 
    
    samalla_koodilla = pd.merge(user_input_DF,master_tulos,  how='inner', right_on =['#CHROM','POS' ], left_on = ['CR', 'pos'] ).drop_duplicates().sort_values(['CR', 'pos'])
    
    #Puuttuvat
    outer_join = user_input_DF.merge(samalla_koodilla, how = 'outer', indicator = True)
    puuttuvat = outer_join[~(outer_join._merge == 'both')].drop('_merge', axis = 1)[['CR', 'pos']]
    slave = master_tulos[~master_tulos['ID'].isin(samalla_koodilla['ID'])]
    slave['POS_plus_1'] = slave['POS'] + 1
    
    plus_1 = pd.merge(puuttuvat,slave,  how='inner',  left_on = ['CR', 'pos'], right_on =['#CHROM','POS_plus_1' ], ).drop_duplicates().sort_values(['CR', 'pos'])

    #Kirjoitetaan CSV datana: 
    valmis = samalla_koodilla.append(plus_1.drop(columns=['POS_plus_1'] )) 
    valmis.to_csv(os.path.join(temp_path, 'Finaalit.csv'), index=False)

    ###Kirjoitetaan VCF datana: 
    vcf_header = get_vcf_header()

    valmis.columns = valmis.columns.str.replace('\n', '')
    valmis = valmis.drop(columns=['CR', 'pos'])
    with open(os.path.join(temp_path, 'test_VCF.vcf'), 'w') as vcf:
        vcf.write(vcf_header)
    valmis.to_csv(os.path.join(temp_path, 'Finaalit.vcf'), sep="\t", header=True, mode='a', index=False)



def get_vcf_header() -> str:
    return f"""##fileformat=VCFv4.2
##fileDate={datetime.today().strftime('%Y%m%d')}
##source=HBP_manual_splicer_by_SJM
##reference=Hg38
##phasing=none
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FILTER=<ID=PASS,Description="All filters passed">
##source=apt-format-result:1.20.6(2.8.6)
##calls-file=/nfs/services2/Projects/PRO100104_FinnGen_SAX/DA/b05/cluster_b05_V1/genotype-inliers/AxiomGT1.calls.txt
##annotation-file=/nfs/falcon/lib/release/Axiom_FinnGen1.r1/analysis/Axiom_FinnGen1.na36.r1.a1.annot.db
##notes="Reporting for indels: For insertions, the reported chromosome position will be the location of first base that was inserted. For deletions, it will be the location of the preceding base even for a deletion within a nucleotide repeat (e.g., TTT[-/T]TTTTT)."
##export-vcf-file=/nfs/services2/Projects/PRO100104_FinnGen_SAX/DA/b05/cluster_b05_V1/genotype-inliers/AxiomGT1.calls.vcf
##FILTER=<ID=axiom_bp,Description="Did not pass Axiom Best Practices">
##contig=<ID=1,species="Homo sapiens">
##contig=<ID=2,species="Homo sapiens">
##contig=<ID=3,species="Homo sapiens">
##contig=<ID=4,species="Homo sapiens">
##contig=<ID=5,species="Homo sapiens">
##contig=<ID=6,species="Homo sapiens">
##contig=<ID=7,species="Homo sapiens">
##contig=<ID=8,species="Homo sapiens">
##contig=<ID=9,species="Homo sapiens">
##contig=<ID=10,species="Homo sapiens">
##contig=<ID=11,species="Homo sapiens">
##contig=<ID=12,species="Homo sapiens">
##contig=<ID=13,species="Homo sapiens">
##contig=<ID=14,species="Homo sapiens">
##contig=<ID=15,species="Homo sapiens">
##contig=<ID=16,species="Homo sapiens">
##contig=<ID=17,species="Homo sapiens">
##contig=<ID=18,species="Homo sapiens">
##contig=<ID=19,species="Homo sapiens">
##contig=<ID=20,species="Homo sapiens">
##contig=<ID=21,species="Homo sapiens">
##contig=<ID=22,species="Homo sapiens">
##contig=<ID=14_GL000009v2_random,species="Homo sapiens">
##contig=<ID=17_KI270909v1_alt,species="Homo sapiens">
##contig=<ID=19_KI270938v1_alt,species="Homo sapiens">
##contig=<ID=22_KI270879v1_alt,species="Homo sapiens">
##contig=<ID=4_GL000008v2_random,species="Homo sapiens">
##contig=<ID=7_KI270803v1_alt,species="Homo sapiens">
##contig=<ID=8_KI270821v1_alt,species="Homo sapiens">
##contig=<ID=MT,species="Homo sapiens">
##contig=<ID=X,species="Homo sapiens">
##contig=<ID=Y,species="Homo sapiens">
"""

########################################################################
########################################################################
def dev_stuff(): 

    def get_metaFiles() -> pd.DataFrame:
        '''Haetaan metadatasta kohdat ja sen perusteella. '''
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
    