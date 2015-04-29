Installtion
=========


# Structure

## Data
In this you have your raw data, download from firehose, via firehose_get script. Here you need some storage, depending on your entities

## pre-process_data
Here are the scripts, which are required to convert raw data to TCGA-Webtools data, better known as .RData

## webapp
Here is final app, also the final data object, which requires some space



# Installation

## First you need the firehose_get script, to download data from firehose, you can with
cd /path/to/TCGA-Webtools
wget http://gdac.broadinstitute.org/runs/code/firehose_get_latest.zip
unzip firehose_get_latest.zip

## Afterwards you can download your entities of interest. This can take some time, in the meaning of hours. The lines below will download all TCGA-Webtools supported data, for Prostate adenocarcinoma (PRAD) and Breast invasive carcinoma (BRCA). For more entities please look at: http://gdac.broadinstitute.org/
cd Data
../firehose_git -b only CopyNumber_Gistic2.Level_4 analysis latest PRAD BRCA
../firehose_git -b only mRNAseq_Preprocess.Level_3 data latest PRAD BRCA
../firehose_git -b only Level_3__within_bioassay_data_set_function__data.Level_3 data latest PRAD BRCA
../firehose_git -b only Mutation_Packager data latest PRAD BRCA

## If your downloads completed, some files need to unziped
find . -type f -name '*.gz' -exec tar xvf {} \;
