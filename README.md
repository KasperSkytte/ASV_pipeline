# Bioinformatic pipeline for processing 16S rRNA amplicon data with exact Amplicon Sequence Variants (ASV's)
The `ASVpipeline.sh` script is a complete linux `BASH` script performing the essential steps in amplicon sequence data processing primarily using [usearch](http://drive5.com/usearch/). 

Outline of the steps performed

For each sample in a `samples` file in the current working directory:
 - Find samples in a folder of sequencing data
 - Filter PhiX sequences from each sample file (by [`filter_phix`](http://drive5.com/usearch/manual/cmd_filter_phix.html))
 - Perform QC filtering of reads based on Phred/Q-scores in the FastQ files (by [`fastq_filter`](http://drive5.com/usearch/manual/cmd_fastq_filter.html))
 - Concatenate all reads from all the samples into one file
 
Then:

 - Find and count unique sequences (by [`fastx_uniques`](http://drive5.com/usearch/manual/cmd_fastx_uniques.html))
 - Perform denoising and chimera filtering of each unique sequence (by [`unoise3`](http://drive5.com/usearch/manual/cmd_unoise3.html))
 - **Search a FASTA file containing ASV's found in previous runs to reuse their ID's (closed reference) enabling cross-study comparison**
 - Predict the taxonomy of the ASV's (by [`sintax`](http://drive5.com/usearch/manual/cmd_sintax.html))
 - Generate ASV table containing read counts of all ASV's for each sample (by [`otutab`](http://drive5.com/usearch/manual/cmd_otutab.html))
 - Cleanup temporary files
  
**References**

Benjamin J Callahan, Paul J McMurdie & Susan P Holmes, 2017, [Exact sequence variants should replace operational taxonomic units in marker-gene data analysis](https://www.nature.com/articles/ismej2017119/). The ISME Journal volume 11, pages 2639–2643
