#!/bin/bash 
#BSUB -J condor_smk
#BSUB -n 1                # number of core(tasks) for parallel jobs          
#BSUB -R rusage[mem=2]    # expected resorce consumption for memory
#BSUB -W 8:00            # run time limit (hours)
#BSUB -o /data/iacobuzc/haochen/condor_pipeline_trial/condor_smk.stdout
#BSUB -eo /data/iacobuzc/haochen/condor_pipeline_trial/condor_smk.stderr

if [ -f ~/.bashrc ] ; then
    . ~/.bashrc
fi

conda activate snakemake

snakemake -s /data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/condor_pipeline.smk \
    --configfile /data/iacobuzc/haochen/condor_pipeline_trial/config_MSK_hz_lilac.yaml \
    --cores all \
    --use-conda \
    --retries 2 \
    --rerun-triggers mtime \
    --rerun-incomplete 