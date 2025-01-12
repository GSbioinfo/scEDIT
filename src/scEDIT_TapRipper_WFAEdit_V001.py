import pandas as pd
from pyfaidx import Fasta
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import FastaIO
import subprocess, sys, os 
import shutil
import re

print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))
if len(sys.argv)<3:
    print("$ipython3 scEDITE_TapRipper_Wrapper_V001.py Sample_sheet.xlsx run_dir")
    sys.exit(0)
samplesheetname=sys.argv[1]
run_dir=sys.argv[2]
print(samplesheetname)

script_dir = os.path.dirname(sys.argv[0]) #os.path.realpath(__file__)
print("directoy of the python script:")
print(script_dir)
run_import_df = pd.read_excel(samplesheetname)
if 'sgRNA_file' not in run_import_df.columns:
    print('Missing sgRNA file location')
    sys.exit(0)
print(run_import_df)

print("*******************************************************************")
print("Warning: scEDIT Version 001 runs multiple cppinputfile.\n")
print("Warning: Running the script will delete and create new dir name "+run_dir+"\n")
print("Make sure that your output directory is not same as previous run.\n")
print("Do you want to continue y/n:  ")
#deci_in=input()
deci_in='y'
#**************Creating run directory for the analysis************
if deci_in.lower()!='y':
    print("Aborting run")
    sys.exit(0)

if run_import_df.empty:
    print("Error:Empty sample information sheet.\n")
    sys.exit(0)

if run_dir[-1:]=="/":
    run_dir=run_dir[:-1]
cwd=os.getcwd()
#run_dir=cwd+'/'+run_dir
if os.path.exists(run_dir):
    shutil.rmtree(run_dir)
    os.makedirs(run_dir)
else:
    os.makedirs(run_dir)
dirname = run_dir#os.path.dirname(run_dir)
#****************************************************************

subjectID_list = run_import_df["SubjectID"].unique()
cmd_list=[]
for subID in subjectID_list:
    
    sample_df = run_import_df.loc[(run_import_df["SubjectID"] == subID)]
    for sample in sample_df["SampleName"].to_list():
        sg_minIden = 60
        ampful_minIden =80
        ampPrim_minIden = 90
        bar_minIdentity = 77
        
        file_list=[]
        file_location = sample_df.loc[(sample_df["SampleName"]==sample)].to_dict(orient="list")
        print(file_location["FastqFile_location_R1_fastqgz"][0])
        print(file_location["Custom_Panel"][0])
        #***********make output directories, trimmed seq output dir for each animal ***************
        outdir=run_dir+'/'+subID+'/'+sample
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        result_dir=outdir+'/result_out'
        if not os.path.exists(result_dir):
            os.makedirs(result_dir)
        
        #*************Create input directory to store vl trimmer input files for each animal********************** 
        runinput_dir=outdir+'/runinput_dir' #--'alt'
        if not os.path.exists(runinput_dir):
            os.makedirs(runinput_dir)
            os.makedirs(runinput_dir+'/tmp_fastqs/')
        temfastq= runinput_dir+'/tmp_fastqs/'

        #*************Create sub output dictory for each animal ********************** 
        myDirectoryList=['fastqs/','dump_fastqs/','intermit_csv/','final_count/', 'CellFastqs/']
        for di in myDirectoryList:
            cmd1='mkdir -p %s' %result_dir+'/'+di
            xxx=subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell =True)
            xxx.communicate()
    
        fastqR1=file_location["FastqFile_location_R1_fastqgz"][0].split('/')[-1]
        fastqR2=file_location["FastqFile_location_R2_fastqgz"][0].split('/')[-1]
        seqstat='seqkit stat '+file_location["FastqFile_location_R1_fastqgz"][0]
        sqkt1=subprocess.Popen(seqstat,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        #sqkt1.communicate()
        index_start=1
        fastqStats=pd.DataFrame()
        for df in pd.read_table(sqkt1.stdout,sep='\s+', header=0,iterator=True,index_col=False):
            fastqStats=pd.concat([fastqStats,df],ignore_index=True)
        #fastqStat=pd.DataFrame(sqkt1.stdout)
        sqkt1.communicate()
        numbread= int(str(fastqStats['num_seqs'][0]).replace(',',''))
        if numbread > 5.5e7:
            minNumRead = 1e7
        else:
            minNumRead = 5e6

        if numbread <minNumRead:
            nsplit = 1
        else:
            nsplit = int(numbread/minNumRead)
        print(nsplit) 
        #parallel 40834
        #
        barfile=file_location["Barcode_File_csv"][0]
        barcode = pd.read_csv(barfile,sep=',',header=None)
        barcode.columns = ["BarSeq","BarID"]

        record_bar=[]
        for bars in barcode.itertuples():
            rec1=SeqRecord(Seq(bars[1]),id="BC"+str(bars[2])+"|"+str(bar_minIdentity),description="",name="",dbxrefs=[])
            record_bar.append(rec1)
 
        genome_file=file_location["Genome_file_fa"][0]
        genome_index=Fasta(genome_file,sequence_always_upper=True,one_based_attributes=False)
        ampliconFile=file_location["TargetAmplicon_File_xlsx"][0]
        amps_df=pd.read_excel(ampliconFile)
        amp_panel = file_location["Custom_Panel"][0]
        amps_df=amps_df.loc[(amps_df['Custom_Panel']==amp_panel)]
        constfast_file = file_location['ConstSeq_File_fasta'][0]
        file_list.append('ConstSeq_file:'+constfast_file)
        
        fprimer=[]
        rprimer=[]
        fulseq =[]
        amp_fileHeader = amps_df.columns

        for amps in amps_df.itertuples():
            if 'chr' in str(amps[1]):
                chromo = str(amps[1])
            else:
                chromo = "chr"+str(amps[1])
            tstart = amps[2]-1
            tend = amps[3] 
            if 'Primer_Ident' in amp_fileHeader:
                seqid = amps[4]+"_"+amps[5]+":"+str(tstart)+"-"+str(tend)+"!"+amps[8]+"|"+str(amps[9])
            else:
                seqid = amps[4]+"_"+amps[5]+":"+str(tstart)+"-"+str(tend)+"!"+amps[8]+"|"+str(ampPrim_minIden)
            seqidtot = amps[4]+"_"+amps[5]+":"+str(tstart)+"-"+str(tend)+"!"+amps[8]+"|"+str(ampful_minIden)
            seqid = seqid.replace(" ","_")
            seqidtot = seqidtot.replace(" ","_")
            seqtot= str(genome_index[chromo][tstart:tend])
            if ('Primer_fwd_seq' in amp_fileHeader) and ('Primer_rev_seq' in amp_fileHeader):
                fdseq = str(amps[6])
                revseq = str(amps[7])    
            else:
                fdseq = str(genome_index[chromo][tstart:tstart+25])
                revseq = str(-genome_index[chromo][tend-50:tend])
                
            fulrec=SeqRecord(Seq(seqtot),id=seqidtot,description="",name=re.sub("_.*","",str(amps[4])),dbxrefs=[])
            fulseq.append(fulrec)
            
            forrec=SeqRecord(Seq(fdseq),id=seqid,description="",name="",dbxrefs=[])
            fprimer.append(forrec)
            
            revrec=SeqRecord(Seq(revseq),id=seqid,description="",name="",dbxrefs=[])
            rprimer.append(revrec)
        
        sgrna_df=pd.read_excel(file_location['sgRNA_file'][0])
        windowseq = []
        sgrnaseq =[]
        nonsg_seq = []
        sgrna_fileHeader = sgrna_df.columns

        for sgrna in sgrna_df.itertuples():
            sg_seq = sgrna.guide_off_target_Seq
            win_seq = sgrna.Lookup_window_sequence 
            #tstart = sgrna[2]-1
            #tend = sgrna[3]-1
            
            for flseq in fulseq:
                seqidwin = 'NOTFOUND'
                seqid = 'NOTFOUND'
                #if('on-target' in sgrnaseq)
                if win_seq in flseq.seq:
                    winstart = flseq.seq.find(win_seq)
                    winend = winstart+len(win_seq)
                    windirect = 'For'
                    seqidwin = flseq.name+":"+windirect+":"+str(winstart)+"-"+str(winend)+":"+str(len(flseq.seq))+'$'+re.sub('\|.*','',flseq.id)+'|'+str(sg_minIden)
                if Seq(win_seq).reverse_complement() in flseq.seq:
                    win_seq = str(Seq(win_seq).reverse_complement())
                    winstart = flseq.seq.find(win_seq)
                    winend = winstart+len(win_seq)
                    windirect = 'Rev'
                    seqidwin = flseq.name+":"+windirect+":"+str(winstart)+"-"+str(winend)+":"+str(len(flseq.seq))+'$'+re.sub('\|.*','',flseq.id)+'|'+str(sg_minIden)
                if sg_seq in flseq.seq:
                    sgstart = flseq.seq.find(sg_seq)
                    sgend = sgstart+len(sg_seq)
                    sgdirect = 'For'
                    seqid = flseq.name+":"+sgdirect+":"+str(sgstart)+"-"+str(sgend)+":"+str(len(flseq.seq))+'$'+re.sub('\|.*','',flseq.id)+'|'+str(sg_minIden)  
                if Seq(sg_seq).reverse_complement() in flseq.seq:
                    sgstart = flseq.seq.find(Seq(sg_seq).reverse_complement())
                    sgend = sgstart+len(sg_seq)
                    sgdirect = 'Rev'
                    seqid = flseq.name+":"+sgdirect+":"+str(sgstart)+"-"+str(sgend)+":"+str(len(flseq.seq))+'$'+re.sub('\|.*','',flseq.id)+'|'+str(sg_minIden)  
                if seqidwin != 'NOTFOUND':
                    fulrec=SeqRecord(Seq(win_seq),id=seqidwin,description="",name=sgrna[5],dbxrefs=[])
                    if Seq(sg_seq) in Seq(win_seq):
                        sgStart = win_seq.find(sg_seq)
                        if sgStart>10:
                            ng_seq = win_seq[:sgStart]
                        else:
                            ng_seq = win_seq[sgStart+len(sg_seq):]
                    if Seq(sg_seq).reverse_complement() in Seq(win_seq):
                        sgStart = Seq(win_seq).find(Seq(sg_seq).reverse_complement())
                        if sgStart>10:
                            ng_seq = win_seq[:sgStart]
                        else:
                            ng_seq = win_seq[sgStart+len(sg_seq):]  
                    partWinSeq = SeqRecord(Seq(ng_seq),id=re.sub('\|[0-9][0-9]','|80',seqidwin),description="",name=sgrna[5],dbxrefs=[])
                    windowseq.append(fulrec)
                    nonsg_seq.append(partWinSeq)
                if seqid != 'NOTFOUND':
                    sgrna_rec=SeqRecord(Seq(sg_seq),id=seqid,description="",name=sgrna[5],dbxrefs=[])
                    sgrnaseq.append(sgrna_rec)

        win_fullseq_file = runinput_dir+'/'+'WIN_fullseq_'+amp_panel+'.fasta'
        handle_fulseq_file=open(win_fullseq_file, 'w')
        fasta_out = FastaIO.FastaWriter(handle_fulseq_file, wrap=None)
        fasta_out.write_file(windowseq)
        handle_fulseq_file.close()
        file_list.append('WINFullSeq:'+win_fullseq_file)
        
        win_partseq_file = runinput_dir+'/'+'WIN_partseq_'+amp_panel+'.fasta'
        handle_partseq_file=open(win_partseq_file, 'w')
        fasta_out = FastaIO.FastaWriter(handle_partseq_file, wrap=None)
        fasta_out.write_file(nonsg_seq)
        handle_partseq_file.close()
        file_list.append('WINPartSeq:'+win_partseq_file)

        sgrna_seq_file = runinput_dir+'/'+'sgRNA_seq_'+amp_panel+'.fasta'
        handle_fulseq_file=open(sgrna_seq_file, 'w')
        fasta_out = FastaIO.FastaWriter(handle_fulseq_file, wrap=None)
        fasta_out.write_file(sgrnaseq)
        handle_fulseq_file.close()
        file_list.append('sgRNASeq:'+sgrna_seq_file)
        amp_fullseq_file = runinput_dir+'/'+'AMP_fullseq_'+amp_panel+'.fasta'
        handle_fulseq_file=open(amp_fullseq_file, 'w')
        fasta_out = FastaIO.FastaWriter(handle_fulseq_file, wrap=None)
        fasta_out.write_file(fulseq)
        handle_fulseq_file.close()
        file_list.append('FullAMPSeq:'+amp_fullseq_file)
        cpp_sgrna = ''
        if "sgRNA_file" in file_location.keys():
            cpp_sgrna = '-g'
            sgrna_file = file_location["sgRNA_file"][0]
            
            gene_sgran_df = pd.read_excel(sgrna_file)
            sgrna_seqs = []
            for sgrnait in gene_sgran_df.itertuples():
                seqid = str(sgrnait[0])+'_sgRNA'
                grna_seq = str(sgrnait[1])
                grna_rec = SeqRecord(Seq(grna_seq),id=seqid,description="",name="",dbxrefs=[])
                sgrna_seqs.append(grna_rec)
            grna_fullseq_file = runinput_dir+'/'+'sgRNA_fullseq_'+amp_panel+'.fasta'
            handle_fulseq_file=open(grna_fullseq_file, 'w')
            fasta_out = FastaIO.FastaWriter(handle_fulseq_file, wrap=None)
            fasta_out.write_file(sgrna_seqs)
            handle_fulseq_file.close()
            #file_list.append('FullsgRNASeq:'+grna_fullseq_file)
        else:
            print("Error:Empty sample information sheet.\n")
            sys.exit(0)
        amp_forprime_file = runinput_dir+'/'+'AMP_Forward_primers_'+amp_panel+'.fasta'
        handle_forprimer_file=open(amp_forprime_file, 'w')
        fasta_out = FastaIO.FastaWriter(handle_forprimer_file, wrap=None)
        fasta_out.write_file(fprimer)
        handle_forprimer_file.close()
        file_list.append('ForwardPrimers:'+amp_forprime_file)
        
        amp_revprime_file = runinput_dir+'/'+'AMP_Revers_primers_'+amp_panel+'.fasta'
        handle_revprimer_file=open(amp_revprime_file, 'w')
        fasta_out = FastaIO.FastaWriter(handle_revprimer_file, wrap=None)
        fasta_out.write_file(rprimer)
        handle_revprimer_file.close()
        file_list.append('ReversePrimers:'+amp_revprime_file)
        handle_record_bar=open(runinput_dir+'/'+'CellBarcode.fasta', 'w')
        fasta_out = FastaIO.FastaWriter(handle_record_bar, wrap=None)
        fasta_out.write_file(record_bar)
        handle_record_bar.close()
        file_list.append('Barcode_file:'+runinput_dir+'/'+'CellBarcode.fasta')
        
        if nsplit < 2:
            #file_list=[]
            tapripper_infile=[]
            file_list.append('R1:'+file_location["FastqFile_location_R1_fastqgz"][0])
            file_list.append('R2:'+file_location["FastqFile_location_R2_fastqgz"][0])
            output_forward_fastq=result_dir+'/'+'fastqs/'+fastqR1.replace('.fastq.gz','_R1_output.fastq')
            output_reverse_fastq=result_dir+'/'+'fastqs/'+fastqR2.replace('.fastq.gz','_R2_output.fastq')
            file_list.append("outR1:"+output_forward_fastq)
            file_list.append("outR2:"+output_reverse_fastq)
            dump_forward_fastq = result_dir+'/'+'dump_fastqs/'+fastqR1.replace('.fastq.gz','_R1_dump.fastq')
            dump_reverse_fastq = result_dir+'/'+'dump_fastqs/'+fastqR2.replace('.fastq.gz','_R2_dump.fastq')
            file_list.append('dumpR1:'+ dump_forward_fastq)
            file_list.append('dumpR2:'+ dump_reverse_fastq)
            file_list.append('ResultTSV:'+result_dir+'/'+'intermit_csv/'+subID+'_'+sample+'_Result_out.tsv')
            maininput_file = runinput_dir+'/'+'cppinput_filelist.txt'
            tapripper_infile.append(maininput_file)
            cppinputfile=open(maininput_file,'w')
            cppinputfile.write('\n'.join(file_list))
            cppinputfile.flush() 
            cppinputfile.close()
        else:
            #cmd = seqkit split2 -1 reads_1.fq.gz -2 reads_2.fq.gz -p 2 -O out -f
            #file_list=[]
            #temfastq=result_dir+'/tmp_fastqs/'
            cmd_split= 'seqkit split2 -1 '+ file_location["FastqFile_location_R1_fastqgz"][0]+' -2 '+file_location["FastqFile_location_R2_fastqgz"][0]+' -p '+str(nsplit)+' -e ".gz" -f '+'-O '+temfastq
            seqsplit=subprocess.Popen(cmd_split,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

            seqsplit.communicate()
            splitFR1=[]
            splitFR2=[]
            tapripper_infile=[]
            for sp in range(1, nsplit+1):
                par_filelist=[]
                splitFR1.append(temfastq+fastqR1.replace('.fastq.gz','')+'.part_'+"{:03d}".format(sp)+'.fastq.gz')
                splitFR2.append(temfastq+fastqR2.replace('.fastq.gz','')+'.part_'+"{:03d}".format(sp)+'.fastq.gz')
                s_splitFR1=temfastq+fastqR1.replace('.fastq.gz','')+'.part_'+"{:03d}".format(sp)+'.fastq.gz'#'|'.join(splitFR1)
                s_splitFR2=temfastq+fastqR2.replace('.fastq.gz','')+'.part_'+"{:03d}".format(sp)+'.fastq.gz'#'|'.join(splitFR2)
            #splitFR2=[]
            #splitFR1=[]
                par_filelist.append('R1:'+s_splitFR1)
                par_filelist.append('R2:'+s_splitFR2)
                output_forward_fastq=result_dir+'/'+'fastqs/'+fastqR1.replace('.fastq.gz','_R1_output.fastq')
                output_reverse_fastq=result_dir+'/'+'fastqs/'+fastqR2.replace('.fastq.gz','_R2_output.fastq')
                par_filelist.append("outR1:"+output_forward_fastq)
                par_filelist.append("outR2:"+output_reverse_fastq)
                dump_forward_fastq = result_dir+'/'+'dump_fastqs/'+fastqR1.replace('.fastq.gz','_R1_dump.fastq')
                dump_reverse_fastq = result_dir+'/'+'dump_fastqs/'+fastqR2.replace('.fastq.gz','_R2_dump.fastq')
                par_filelist.append('dumpR1:'+ dump_forward_fastq)
                par_filelist.append('dumpR2:'+ dump_reverse_fastq)
                par_filelist.append('ResultTSV:'+result_dir+'/'+'intermit_csv/'+subID+'_'+sample+'_Result_outpart_'+"{:03d}".format(sp)+'.tsv')
                maininput_file = runinput_dir+'/'+'cppinput_filelist_part_'+"{:03d}".format(sp)+'.txt'
                tapripper_infile.append(maininput_file)
                cppinputfile=open(maininput_file,'w')
                cppinputfile.write('\n'.join(file_list+par_filelist))
                cppinputfile.flush() 
                cppinputfile.close()
        
        if nsplit < 2:
            if cpp_sgrna == '-g':
                cmd = script_dir+'/./singleCell_HLWFA.out -f $infile'+' -t '+ str(file_location['Numb_Threads'][0])+ ' -s no' + ' -g yes'
            else:
                cmd = script_dir+'/./tapripper -f $infile'+' -t '+ str(file_location['Numb_Threads'][0])+ ' -s no' + ' -g no'
        else:
            if cpp_sgrna == '-g':
                cmd = script_dir+'/./singleCell_HLWFA.out -f $infile'+' -t '+ str(file_location['Numb_Threads'][0])+ ' -s yes' + ' -g yes'
            else:
                cmd = script_dir+'/./tapripper -f $infile'+' -t '+ str(file_location['Numb_Threads'][0])+ ' -s yes' + ' -g no'

        #cmd = '/home/gajendra/Dropbox/Tapestri_workflow/./tapripper -f $infile'+' -t '+ str(file_location['Numb_Threads'][0])+ ' -s yes'
        cmd_list.append(cmd)
        bash_cmds= []
        bash_cmds.append('#!/bin/bash')
        bash_cmds.append('#sh bash_saminfo.sh ')
        bash_cmds.append('date')
        bash_cmds.append('for infile in '+' '.join(tapripper_infile))
        bash_cmds.append('do')
        bash_cmds.append(cmd)
        bash_cmds.append('done')
        zipfile_list=[output_forward_fastq,output_reverse_fastq,dump_forward_fastq,dump_reverse_fastq]
        zipstr = ' '.join(zipfile_list)
        
        bash_cmds.append('for tozip in '+zipstr)
        bash_cmds.append('do')
        bash_cmds.append('echo Compressing file $tozip')
        bash_cmds.append('pigz -9 -f -p12 '+'$tozip')
        bash_cmds.append('done')
        bash_cmds.append('ipython3 '+script_dir+'/count_barcode_edit_correct.py '+result_dir+'/'+'intermit_csv/ '+ subID+'_'+sample)
        fileter_countList = result_dir+'/'+'final_count/Filter_barcode_'+ subID+'_'+sample+'_list.txt'
       # bash_cmds.append('/home/gajendra/Dropbox/Tapestri_workflow/./fastqSpliter -r R1 -f '+output_forward_fastq+' -c '+fileter_countList+' -o '+result_dir+'/'+'bamFiles/'+' >'+result_dir+'/'+subID+'_'+sample+'_r1.out')
       # bash_cmds.append('/home/gajendra/Dropbox/Tapestri_workflow/./fastqSpliter -r R2 -f '+output_reverse_fastq+' -c '+fileter_countList+' -o '+result_dir+'/'+'bamFiles/'+' >'+result_dir+'/'+subID+'_'+sample+'_r2.out')
       # bash_cmds.append('pigz -9 '+result_dir+'/'+'bamFiles/A*.fastq')
       # bash_cmds.append('pigz -9 '+result_dir+'/'+'bamFiles/G*.fastq')
       # bash_cmds.append('pigz -9 '+result_dir+'/'+'bamFiles/T*.fastq')
       # bash_cmds.append('pigz -9 '+result_dir+'/'+'bamFiles/C*.fastq')
       # bash_cmds.append('bowtie2-build -f '+amp_fullseq_file)
       # bash_cmds.append('for fq in '+result_dir+'/'+'bamFiles/*R1.fastq.gz')
       # bash_cmds.append('do')
       # bash_cmds.append('bowtie2 -x '+amp_fullseq_file+' -q --end-to-end -p 4 -1 $fa -2 ${fq/R1.fastq/R2.fastq} |samtools view -bS >${fq/_R1.fastq.gz/.bam}')
        #bash_cmds.append('done')
        bash_cmds.append('rm -f '+temfastq+'*.gz')
        bash_cmds.append('date')
        bashrunfile=dirname+'/'+'bash_'+subID+'_'+sample+'.sh'
        bashwrite = open(bashrunfile,'w')
        bashwrite.write('\n'.join(bash_cmds))
        bashwrite.close() 



#bowtie2 -x /home/gajendra/AML_80analysis/Patient03/SRR14107105/runinput_dir/AMP_fullseq_CO86 -q --end-to-end -p 4 -S TTGTATCACCAATGTGCT.sam -1 TTGTATCACCAATGTGCT_R1.fastq.gz -2 TTGTATCACCAATGTGCT_R2.fastq.gz




'''
trim_seq=subprocess.Popen(cmd_list[0], stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell =True)
stdout_data, stderr_data = trim_seq.communicate()
print(stdout_data)
if trim_seq.returncode != 0:
    print("%r failed, stderr %r" % (cmd_list[0], stderr_data))
    sys.exit(0)
zipfile_list=[output_forward_fastq,output_reverse_fastq,dump_forward_fastq,dump_reverse_fastq]
for tozip in zipfile_list:
    zipcmd='gzip -f9 '+tozip
    zippro=subprocess.Popen(zipcmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell =True)
    stdout_data, stderr_data = zippro.communicate()
    if zippro.returncode != 0:
        print("%r failed, stderr %r" % (run_trim_cmd, stderr_data))
        sys.exit(0)

file_list[:]=[]
'''
