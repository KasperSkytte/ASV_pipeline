# Bioinformatic pipeline for processing 16S rRNA amplicon data with exact Amplicon Sequence Variants (ASV's)
The `ASVpipeline.sh` script is a complete linux `BASH` script performing the essential steps in amplicon sequence data processing using only [usearch](http://drive5.com/usearch/).

**The script processes forward reads ONLY!**

## Outline of the steps performed

For each sample in a `samples` file in the current working directory:
 - Find samples in a folder of sequencing data (FASTQ files)
 - Filter PhiX sequences from each sample file (by [`filter_phix`](http://drive5.com/usearch/manual/cmd_filter_phix.html))
 - Perform QC filtering of reads based on Phred/Q-scores in the FastQ files (by [`fastq_filter`](http://drive5.com/usearch/manual/cmd_fastq_filter.html))
 - Concatenate all reads from all the samples into one file
 
Then:

 - Find and count unique sequences (by [`fastx_uniques`](http://drive5.com/usearch/manual/cmd_fastx_uniques.html))
 - Perform denoising and chimera filtering of each unique sequence (by [`unoise3`](http://drive5.com/usearch/manual/cmd_unoise3.html))
 - **_Search a FASTA file containing ASV's found in previous runs to reuse their ID's (closed reference) enabling cross-study comparison_**
 - Predict the taxonomy of the ASV's (by [`sintax`](http://drive5.com/usearch/manual/cmd_sintax.html))
 - Generate ASV table containing read counts of all ASV's for each sample (by [`otutab`](http://drive5.com/usearch/manual/cmd_otutab.html))
 - Cleanup temporary files

# Requirements
Only `usearch11` is needed and expected to be in `$PATH`, the rest are just general GNU Linux tools. Should work straight out of the box on all Ubuntu/Debian x86 systems and the likes. 

# Installation and usage
Just download the `ASVpipeline.sh` script, or clone this git repo, and run it with `bash ASVpipeline.sh -h`:

```
$ bash ASVpipeline.sh -h
ASV (zOTU) pipeline to process amplicon data. Forward reads only.
Version: 1.3.13
Options:
  -h    Display this help text and exit.
  -i    (Required) Path to file containing sample ID's to find and process. One sample ID per line. 
          (Default: samples)
  -d    Path to taxonomic database. This database will also be used for prefiltering ASV's (<60% id). 
          (Default: )
  -f    (Required) Path to folder containing fastq files (will be searched recursively). 
          (Default: /raw_data/sequences/Illumina/)
  -o    (Required) Output folder. 
          (Default: /user_data/ksa/software/lib/asv_pipeline/output)
  -a    Path to ASV database (a fasta file) with sequences whose names will be reused by exact matching. The ASV's (zOTU's) generated by this script will then be renamed on exact match, any unmatched ASV's will get the prefix 'newASV'. 
          (Default: )
  -k    (flag) Keep all intermediate files as well as all (decompressed) input fastq files, i.e. don't delete ANYTHING.
  -t    Max number of threads to use. 
          (Default: all available except 2)
  -v    Print version and exit.

Additional options can be set by exporting environment variables before running the script:
  - chunksize: Max number of threads to use for each job run in parallel when generating abundance table.
          (Default: 5)
  - samplesep: Separator after sample names and the rest of the filename of the fastq files.
          (Default: _)
  - newASVprefix: Prefix to use for ASVs that don't map to anything in the ASV database.
          (Default: newASV)
  - logfilename: Filename of the logfile. It will be placed in the chosen output folder.
          (Default: asvpipeline_20221014_104758.txt)
```

Initially running `chmod +x ASVpipeline.sh` can be needed if you run into a permission error.

## References

Benjamin J Callahan, Paul J McMurdie & Susan P Holmes, 2017, [Exact sequence variants should replace operational taxonomic units in marker-gene data analysis](https://www.nature.com/articles/ismej2017119/). The ISME Journal volume 11, pages 2639–2643
