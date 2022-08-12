#! /bin/sh

cellranger count --id=${sampleID}_${runID} \
                   --transcriptome=nfu_genome \
                   --fastqs=../../../../raw_data/scRNAseq_TK_headkidney_youngold/ \
                   --sample=4443_${sampleID}_${runID}_SI-GA-${sampleID}7
