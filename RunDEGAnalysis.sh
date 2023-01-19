#!/bin/bash

source /broad/software/scripts/useuse;
reuse R-4.1;

while getopts s:o:d:r:i: flag; do
    case "${flag}" in
        s) seuratObject=${OPTARG};;
        o) outputDir=${OPTARG};;
        d) dataDir=${OPTARG};;
        r) resolution=${OPTARG};;
        i) integrate=${OPTARG};;

    esac
done

Rscript $dataDir/AutomaticSeuratScripts/DEGAnalysis.R --seuratObject $seuratObject --outputDir $outputDir --integrate $integrate --resolution $resolution