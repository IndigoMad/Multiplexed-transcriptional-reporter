# Multiplexed-transcriptional-reporter

This repository contains the quantified data and processing scripts (MATLAB) of literature "Temporal analysis of mammalian genetic circuits using multiplexed transcriptional imaging".


The **figure_scripts** folder contains the scripts for create figure from the data. Other folders contains quantified data and processing scripts of different datasets, except the data generated by simulation, which could be generated by running the simulation scripts. The datasets in those folders are listed below:

**single_barcode_transient**: the dataset of the cells with five transient transfected single barcode respectively.

**single_barcode**: the dataset of the cells with five stable transfected single barcode respectively.

**double_barcode**: the dataset of the cells with double barcode of 1:0/0:1 or 1:2/2:1.

**ABC**: the dataset of the cells with the synthetic three-gene cascade circuit

**ABCDE**: the dataset of the cells with the five-gene branched cascade circuit

**ABC_simulation**: the simulation results corresponding to the synthetic three-gene cascade circuit


The mat files in these folders are the quantified data generated by quantified software and manually classification, recording the infomations including mitosis time of cell, the fluorescent inthensity of the gene loci and a description about the detail informations of the dataset. The data are restored as four level: **cell_set**, **cell**, **spot_set** and **spot**. **Cell_sets** are sets of **cells**, which were collected under the same condition. **Cells** are the data of every single cells, including the mitosis time of the cell and several **spot_sets**. **Spot_sets** contains **spots** that classified as the same group, usually by ratio of the fluorescent intensity. **Spots** represent the fluorescent spots form by transcription reporters in the cells, namely, the gene loci that were transcribing mRNAs. A **spot** record contain the intensity of CFP and mCherry fluorscence along time.


The MATLAB scripts (filename.m) are used to process the data and to add more processed informations about bursts, regulation, etc. into the datasets. Every folders have a **total_workflow.m** as guidence of the using of these scripts. Directly running **total_workflow.m** would finish all the processing.


What's more, a **json** file is provided in every folders, these json files contains all processed data of the datasets. A user have no access to MATLAB can read this file with any other language or online/local json reader, for json is a common storage format. These files were encoded with **UTF-8**.  
