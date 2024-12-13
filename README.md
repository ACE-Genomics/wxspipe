# WESpipe

Parallel processing for whole exome sequencing (WES) pipeline

![Individual WES pipeline](wes_pipe.png)

## Before you go

This scripts use SLURMACE library. To install it just run the *install\_slurmace.sh* script provided here. This should install the perl module into your local PERL5 directory.

## scripts

The whole project are some scripts that use GATK and other tools to run WES pipeline. 

   * fasta2cram.pl : Transform FASTA files to aligned HG38 CRAM files
   * fasta2vcf.pl : Launch WES pipeline from FASTA to gVCF files
   * realign.pl : Realign BAM files from B37 to HG38
   * bam2vcf.pl : Launch WES pipeline from BAM (or CRAM) to gVCF files
   * vcfjoint.pl : Make a joint call from gVCF files
   * parse\_reports.pl : Parse the QC files and make a report
   * wtf\_tranches.pl : Help to inspect joint call tranches

## GET your references

First of all you will need to download your references from GATK bucket or something similar. I just downladed everything I need from  https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/. 

Now you need to know your data. FASTA or CRAM files should be together with some libraries (what GATK call _bait_ and _target_). You should bulid the _interval-list_ files (see https://gatk.broadinstitute.org/hc/en-us/articles/360036726611-BedToIntervalList-Picard).

and remember to index all the ref vcf files,

```
while read -r vcf; do gatk IndexFeatureFile -I ${vcf}; done < toindex.txt
```

So far is all ready to run the WES pipeline.

## Just go

Now you are ready to go but first you will need to edit the scripts and modify some initial variables like source and output directories, location of interval lists, etc. Then you find the script that suits your needs and try it!

### Some options

The scripts has also some basic optional input options for do some testing in your sample,

   * -c : especify a file with a subsebt of the subjects to analyze, run the script only on these subjects
   * -o : especify where to storage the output
   * -s : especify where to look for subject's fasta files
   * -g : for debugging pourposes, do not remove intermediate temporary files
   * -t : actually do not run nothing but create the full SLURM structure, usefull to inspect the slurm script that will be send into the cluster

## TO DO

   * More docs!


