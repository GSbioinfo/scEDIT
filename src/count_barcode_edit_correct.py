import pandas as pd
from pyfaidx import Fasta
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import FastaIO
import subprocess, sys, os, re
import shutil
import sqlite3
import numpy as np
from scipy.misc import derivative
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import CubicSpline
import math
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

if len(sys.argv)<3:
    print("$ipython3 count_barcode.py intermit_csv_dir subject_samplID")
    sys.exit(0)

run_dir=sys.argv[1]
subSampID = sys.argv[2]
#BaF3_NSG_sc/SRR12717666/result_out/
print("Calculating edit count and final barcode counting")
loc_count_file = run_dir
filter_out_dir= run_dir.replace("/final_csv/","/final_count/")
filter_out_dir= filter_out_dir.replace("/intermit_csv/","/final_count/")
if not os.path.exists(filter_out_dir):
            os.makedirs(filter_out_dir)
            #os.makedirs(runinput_dir+'/tmp_fastqs/')
barfile_list= os.listdir(loc_count_file)
editfile_list = os.listdir(loc_count_file)
barfile_list = [barfi for barfi in barfile_list if barfi.find('.tsv')>0]
barfile_list.sort()
bartable_db = sqlite3.connect(filter_out_dir+'bartable_80_rank_db')
bls=bartable_db.cursor()
bls.execute('''DROP TABLE IF EXISTS barlist''')
editfile_list = [editfil for editfil in editfile_list if re.match('EditDetails_*.*.txt', editfil)]
editfile_list.sort()
index_start = 1 
j=0
dfe=pd.read_table(loc_count_file+editfile_list[0],sep='\t')

#print(len(dfe.columns.to_list()))
column_names2 = dfe.columns.to_list()
for fil in editfile_list[1:]:
    df1 = pd.read_table(loc_count_file+fil,sep='\t',header=0, index_col=False)
    dfe = pd.concat([dfe,df1], axis=0)
#df=df.sort_values(by=['RowSum'],ascending=False).reset_index()
#df=df.groupby('CellID').sum().sort_values(by=['RowSum'],ascending=False).reset_index()
dfe['CB_AMP']= dfe['CellBarcode'].astype(str)+'!' + dfe['AmpID'].astype('str')
dfg = dfe[['AmpEditCount', 'CorrectedEditCount','CB_AMP']].groupby('CB_AMP').sum()
dfg = dfg.loc[dfg['AmpEditCount']>1]
dfg_str = dfe[['EditStr|Loc,','CB_AMP']].groupby('CB_AMP')['EditStr|Loc,'].apply('$'.join).reset_index()
dfg= pd.merge(dfg,dfg_str, on='CB_AMP')
dfg['AdjeditCount'] = dfg['CorrectedEditCount']
dfg['Locus1'] = 0
dfg['Locus2'] = 0
dfe =[]
for index, row in dfg.iterrows():   
    listsplit = []
    for xlis in str(row['EditStr|Loc,']).split('$'):
        for mlix in xlis.split(','):
            if len(mlix)>0:
                listsplit.append(mlix)
        
    cigar_list = [edtxt.split('|')[0]+'!'+edtxt.split('|')[1]+'!'+edtxt.split('|')[2] for edtxt in listsplit]
    cig_count = [int(edtxt.split('|')[3]) for edtxt in listsplit]
    data = {'CIG': cigar_list, 'CNT': cig_count}
    cigcnt_df = pd.DataFrame(data).groupby('CIG').sum()
    cigcnt_df = cigcnt_df.reset_index()
    cigcnt_df = cigcnt_df.loc[cigcnt_df['CNT']>1]
    if cigcnt_df.empty:
        continue
    
    if len(cigcnt_df.index)>1:
        cigcnt_df.sort_values(by="CNT",ascending=False,inplace=True)
        dfg['AdjeditCount'].iat[index] = cigcnt_df['CNT'].sum()
        dfg['Locus1'].iat[index] = int(cigcnt_df['CNT'].iat[0])
        dfg['Locus2'].iat[index] = int(cigcnt_df['CNT'].iat[1])
    else:
        dfg['AdjeditCount'].iat[index] = int(cigcnt_df['CNT'].iat[0])
        dfg['Locus1'].iat[index] = int(cigcnt_df['CNT'].iat[0])

dfg[['CellID','AmpIDEdite']]=dfg['CB_AMP'].str.split("!",n=1,expand=True)
dfg[['AmpID','Edite']]=dfg['AmpIDEdite'].str.split("!",n=1,expand=True)
dfg['NewCount'] = dfg['AdjeditCount'] - dfg['CorrectedEditCount']
dfg[['CellID','AmpIDEdite','AmpEditCount','AdjeditCount','Locus1','Locus2','EditStr|Loc,']].to_csv(filter_out_dir+'Loci_Edits_'+subSampID+'_data.txt',sep='\t',index=False)
AdjustedEdit_df = dfg.pivot(index='CellID',columns='AmpIDEdite',values='NewCount').fillna(0).reset_index()
AdjustedAmp_df = (-1*dfg.pivot(index='CellID',columns='AmpID',values='NewCount')+0).fillna(0).reset_index()
    
Adjusted_df = pd.merge(AdjustedAmp_df,AdjustedEdit_df, on='CellID')
Adjcolname=Adjusted_df.columns.to_list()

Adjusted_df = Adjusted_df[['CellID']+sorted(Adjcolname[1:])]
#Adjusted_df['RowSum'] = 0
Adjcolname=Adjusted_df.columns.to_list()

index_start = 1 
j=0
df=pd.read_table(loc_count_file+barfile_list[0],sep='\t')

print(len(df.columns.to_list()))
column_names2 = df.columns.to_list()
for fil in barfile_list[1:]:
    df1 = pd.read_table(loc_count_file+fil,sep='\t',header=0, index_col=False)
    df=pd.concat([df,df1], axis=0)
#df=df.sort_values(by=['RowSum'],ascending=False).reset_index()
temp_df = df[Adjcolname].copy(deep=True)

temp_df = pd.concat([temp_df,Adjusted_df], axis=0)
temp_df=temp_df.groupby('CellID').sum().reset_index()

refgeneAMPs=df.columns.to_list()
refgeneAMPs = [rfgen for rfgen in refgeneAMPs if rfgen.find("!HostGene")>0]
refGen_df = df[['CellID']+sorted(refgeneAMPs)+['RowSum']].copy(deep=True)
refGen_df=refGen_df.groupby('CellID').sum().reset_index()
df= pd.merge(temp_df,refGen_df, on="CellID")
df=df.sort_values(by=['RowSum'],ascending=False)

df['RowSumRank_max']=df['RowSum'].rank(method = 'max', ascending=False)
df['RowSumRank_min']=df['RowSum'].rank(method = 'min', ascending=False)
df['RowSumRank']=df['RowSum'].rank(ascending=False)
df['Rank_dff']=df['RowSumRank'].diff()
df['RowSum_CSpline']=df['RowSum'].interpolate(method='cubicspline')
df['RowSum_diff']=df['RowSum_CSpline'].diff()
df['dydx']=df['RowSum_diff']/df['Rank_dff']
df['ddydx']=df['RowSum_diff'].diff()/df['Rank_dff']
df.to_sql('barlist', bartable_db, if_exists='append')
df_orig = df.copy(deep=True)
#df=[]
df.columns

amp_filter = pd.DataFrame(df[df.columns[[':'in x for x in df.columns]]].sum(axis=0))


#read_cutoff = len(amp_filter)*5

amp_filter.columns = ['SUM']
pd.DataFrame(amp_filter['SUM'].describe()).transpose().to_csv(filter_out_dir+'Total_ampCount_'+subSampID+'.txt',sep='\t',index=True)

select_amps=amp_filter.loc[amp_filter['SUM']>amp_filter['SUM'].quantile([0.05]).iloc[0]]
select_amps.to_csv(filter_out_dir+'FinalCutOffAmpLsit_'+subSampID+'.txt',sep='\t',index=False)
read_cutoff = len(select_amps.index.to_list())*5

df= df[['RowSum', 'RowSumRank_max','RowSumRank_min','RowSumRank','Rank_dff','RowSum_CSpline','RowSum_diff','dydx','ddydx']]


ddypos=df.loc[df['ddydx']<0]


# Log10 transform RowSum and 
df['LogRowSum']= df['RowSum'].apply(math.log10)

df['LogRowSum_diff']=df['LogRowSum'].diff()

df['LogRowSum_dd']=df['LogRowSum_diff'].diff()
df_sumrank = df[['RowSumRank','RowSum']].copy(deep=True)
plt.loglog(df_sumrank['RowSumRank'], df_sumrank['RowSum'], label='Smoothed Data',linewidth=4, color='black')
plt.xlabel("Rank by read count")
plt.ylabel("Read count per cell")
plt.legend(bbox_to_anchor=(1.55, 1.0))
plt.savefig(filter_out_dir+subSampID+'_cell_Precall_plot.png')

df = df[df['RowSum']>read_cutoff]

df['RowSumRank'].astype('int32').max()*0.1
slide_window = int(df['RowSumRank'].astype('int32').max()*0.1)#500
slider = int(df['RowSumRank'].astype('int32').max()*0.1)

cell_call_df = pd.DataFrame(columns =df.columns)
while slider <= df['RowSumRank'].astype('int32').max():
    if slider == slide_window:
        slider_shiter= 0
    else:
        slider_shiter = int(df['RowSumRank'].astype('int32').max()*0.01) 

    df_slice = df.loc[(df['RowSumRank']>(slider-slider_shiter)) & (df['RowSumRank']<(slider+slide_window))]
    slider=slider+slide_window
    ddline_rankdf=df_slice.sort_values(by='LogRowSum_diff')[0:3]
    df_select=df_slice.loc[(df_slice['RowSum']<=ddline_rankdf['RowSum'].max())&(df_slice['RowSum']>=ddline_rankdf['RowSum'].min())]
    df_select=df_select.loc[(df_select['LogRowSum_diff']<0)&(df_select['LogRowSum_dd']>0)]
    cell_call_df=pd.concat([cell_call_df,df_select])
if cell_call_df.empty:
    print('Auntomated cell calling failed due to nature of the data.\n')
    print('Need to perfrom cell calling manually\n')
    print('Writing full data in '+filter_out_dir+' \n')
    df_orig[['CellID'] + amp_filter.index.to_list()+['RowSum']].to_csv(filter_out_dir+'Orignal_'+subSampID+'_data.tsv.gz',sep='\t',index=False,compression='gzip')
    sys.exit()
cell_call_df=cell_call_df.loc[cell_call_df['RowSumRank']<(df['RowSumRank'].astype('int32').max()-slide_window)]
cell_call_df=cell_call_df.sort_values(by=['LogRowSum_diff','RowSum'],ascending=[True,False])
ycut=cell_call_df['RowSum'].iloc[0]
xcut=cell_call_df['RowSumRank'].iloc[0]

plt.loglog(df_sumrank['RowSumRank'], df_sumrank['RowSum'], label='Smoothed Data',linewidth=4)
plt.loglog(xcut, ycut,marker='.', ls='none', ms=10, color='red')
plt.loglog([1,xcut],[ycut,ycut],color='green')
plt.loglog([xcut,xcut],[1,ycut],color="green")

#plt.yscale("log")
#plt.xscale("log")
plt.xlabel("Rank by read count")
plt.ylabel("Read count per cell")
plt.legend(bbox_to_anchor=(1.55, 1.0))
plt.savefig(filter_out_dir+subSampID+'_cell_call_plot.png')
#plt.show()


called_df = df_orig.loc[df_orig['RowSum']>ycut]
called_df=called_df[['CellID'] + amp_filter.index.to_list()+['RowSum']].copy(deep=True)

called_df.to_csv(filter_out_dir+'filtered_'+subSampID+'_data.tsv',sep='\t',index=False)
called_df['CellID'].to_csv(filter_out_dir+'Filter_barcode_'+subSampID+'_list.txt',sep='\t',index=False,header=False)
df_orig[['CellID'] + amp_filter.index.to_list()+['RowSum']].to_csv(filter_out_dir+'Orignal_'+subSampID+'_data.tsv.gz',sep='\t',index=False,compression='gzip')

bartable_db.close()
