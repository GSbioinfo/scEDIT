#!/usr/bin/env nextflow
params.outputDir = "working/NF_run_testing/"
params.dataDir = "/scratch/"
params.sampleFileLoc = "Tapestri_CRISPRCas9/TapCRISPRCas_sample_Edit_NF_file.xlsx"
params.workDir = "/working/"
params.genomeDir = "/home/gajendra/refgenomes/"
params.scriptDir = "/home/gajendra/Dropbox/Github_repos/scEDIT/"
params.bashFiles = "${params.outputDir}bash*.sh"
ch1 = Channel.of(params.bashFiles)
process importDockerImage {
    input:
    path dockerImageTar

    output:
    path "imported_image_id.txt"

    script:
    """
    # Import the Docker image from the tar file
    docker load < $dockerImageTar | awk '/Loaded image:/ {print \$3}' > imported_image_id.txt
    """
}

process runDockerContainer {
    input:
    path imageFile

    output:
    path "container_output.log"

    script:
    """
    # Get the image ID
    IMAGE_ID=\$(cat ${imageFile})
    # Run the Docker container with the provided Python script
    docker run --rm -v ${params.dataDir}:/docHome/scratch/ \
        -v ${params.workDir}:/docHome/working/ \
        -v ${params.genomeDir}:/docHome/refgenomes/ \
        -v ${params.scriptDir}:/docHome/scEDIT/ \
        \${IMAGE_ID} \
        bash -c "
        source activate scEDIT &&\
        ipython3 scEDIT/src/scEDIT_TapRipper_WFAEdit_V001.py \
        /docHome/scratch/${params.sampleFileLoc} \
        /docHome/${params.outputDir} \
        " > container_output.log
        sleep 5
    """
}

process runbashScripts {
    input:
    path flag_file
    val STR
    path imageFile
    

    output:
    path "container_output.log"

    script:
    """
    # Get the image ID
    IMAGE_ID=\$(cat ${imageFile})
    # Run the Docker container with the provided Python script
    docker run --rm -v ${params.dataDir}:/docHome/scratch/ \
        -v ${params.workDir}:/docHome/working/ \
        -v ${params.genomeDir}:/docHome/refgenomes/ \
        -v ${params.scriptDir}:/docHome/scEDIT/ \
        \${IMAGE_ID} \
        bash -c "
        source activate scEDIT &&\
        sh ${STR}  \
        " > container_output.log
        sleep 5  
    """
}
workflow {
    // Define inputs
    dockerImageTar = file("/scratch/DockerImages/interactive_cond_scedit.tar")   // Replace with your Docker tar image file
    
    //pythonScript = file("script.py")           // Replace with your Python script

    // Execute processes
    importDockerImage(dockerImageTar)
    //runDockerContainer(importDockerImage.out, pythonScript)
    step1_result = runDockerContainer(importDockerImage.out)
    ch1.view()
    runbashScripts(step1_result, ch1,importDockerImage.out)
}
