#!/usr/bin/env nextflow
params.outputDir = "working/NF_run_testing/"
params.dataDir = "/scratch/"
params.sampleFileLoc = "Tapestri_CRISPRCas9/TapCRISPRCas_sample_Edit_NF_file.xlsx"
params.workDir = "/working/"
params.genomeDir = "/home/gajendra/refgenomes/"
params.scriptDir = "/home/gajendra/Dropbox/Github_repos/scEDIT/"
params.bashFiles = "/docHome/${params.outputDir}bash*.sh"


process runDockerContainer {
    container 'interactive_cond_scedit:1.0'
    containerOptions "--rm -v ${params.workDir}:/docHome/working/ -v ${params.dataDir}:/docHome/scratch/ -v ${params.genomeDir}:/docHome/refgenomes/ -v ${params.scriptDir}:/docHome/scEDIT/"
    
    output:
    path "bashList.txt"

    script:
    """
    # activate conda env
    source activate scEDIT
    ipython3 scEDIT/src/scEDIT_TapRipper_WFAEdit_V001.py \
    /docHome/scratch/${params.sampleFileLoc} \
    /docHome/${params.outputDir}
    ls /${params.outputDir}bash*.sh > bashList.txt
    sleep 5
    """
}
process modify_bashlist{
    input:
    path infile

    script:
    """
    # Import the Docker image from the tar file
    modifed_bashlist=\$(cat ${infile}  | tr '\\n' ' ')

    #Show created list
    for bsh in \$modifed_bashlist; do echo \$bsh; done
    sleep 5
    """
}

process runbashScripts {
    container 'interactive_cond_scedit:1.0'
    containerOptions "--rm -v ${params.workDir}:/docHome/working/ -v ${params.dataDir}:/docHome/scratch/ -v ${params.genomeDir}:/docHome/refgenomes/ -v ${params.scriptDir}:/docHome/scEDIT/"
    
    input:
    path infile
    
    output:
    path "container_output.txt"

    script:
    """
    # activate conda env
    source activate scEDIT
    # get bash list
    MOD_BF=\$(cat ${infile} | tr '\\n' ' ')
    for bsh in \${MOD_BF}
    do 
    sh /docHome\${bsh}  
    done 
    > container_output.txt
    sleep 5  
    """
}
workflow {
    // Define inputs
    dockerImageTar = file("/scratch/DockerImages/interactive_cond_scedit.tar")   // Replace with your Docker tar image file
    
    // Import or load tar image
    //importDockerImage(dockerImageTar)
    
    // Run docker image and generate bash file using the  
    runDockerContainer()
    
    // Test bash file before starting the read processing  
    modLis = modify_bashlist(runDockerContainer.out)
    
    // Run bash script to process the reads and generate read counts data
    runbashScripts(runDockerContainer.out)
}
