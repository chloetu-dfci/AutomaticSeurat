#!/bin/bash

source /broad/software/scripts/useuse;
reuse R-4.1;

while getopts n:f:o:d:p:g:i:s: flag; do
    case "${flag}" in
        n) sampleNames=${OPTARG};;
        f) fileNames=${OPTARG};;
        o) outputDir=${OPTARG};;
        d) dataDir=${OPTARG};;
        p) pdfName=${OPTARG};;
        g) geneVec=${OPTARG};;
        i) integrate=${OPTARG};;
        s) species=${OPTARG};;
    esac
done

Rscript $dataDir/AutomaticSeuratScripts/AutomaticSeurat.R --sampleNames $sampleNames --fileNames $fileNames --outputDir $outputDir --dataDir $dataDir --pdfName $pdfName --geneVec $geneVec --integrate $integrate --species $species