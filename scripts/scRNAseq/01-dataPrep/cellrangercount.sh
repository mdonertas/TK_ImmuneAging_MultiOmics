#! /bin/sh

cellranger count --id=${sampleID} \
                   --transcriptome=nfu_genome \
                   --fastqs=../../../../raw_data/scRNAseq_TK_headkidney_youngold/ \
                   --sample=4443_${sampleID}_run636_SI-GA-${sampleID}7,4443_${sampleID}_run640_SI-GA-${sampleID}7,4443_${sampleID}_run641_SI-GA-${sampleID}7
