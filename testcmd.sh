#!/bin/bash

htseq-count -n 10 -f bam -m intersection-strict --stranded=no -a 10 -t exon -i gene_name test.bam test.gtf
