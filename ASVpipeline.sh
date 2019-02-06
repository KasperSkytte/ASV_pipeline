#!/bin/bash
# Set error check
set -e
set -o pipefail

MAX_THREADS=${1:-$((`nproc`-2))}
SEQPATH=/home/kapper/dragun/space/sequences/
TAXDB=/home/kapper/Dropbox/AAU/PhD/Projects/ESV\ pipeline/runs/MiDAS3.0/midas30_20190109/output/ESVs_w_sintax.fa
ASVDB=/home/kapper/Dropbox/AAU/PhD/Projects/ASV_test/midasfull/200bp/ASVs.R1.fa
SAMPLESEP="_"

rm -rf rawdata/
rm -rf phix_filtered/
mkdir -p rawdata
mkdir -p phix_filtered
mkdir -p phix_filtered/tempdir

echoWithDate() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')]: $1"
}

echoWithDate "Running ASV pipeline (max threads: $MAX_THREADS)..."
echoWithDate "Finding samples, filtering PhiX and bad reads, truncating to 250bp..."
cat samples | tr "\r" "\n" | sed -e '$a\' | sed -e '/^$/d' -e 's/ //g' > samples_tmp.txt

NSAMPLES=$(wc -w < samples_tmp.txt)
while ((i++)); read SAMPLE
  do
  echo -ne "Processing sample: $SAMPLE ($i / $NSAMPLES)\r"
  find "$SEQPATH" -name $SAMPLE$SAMPLESEP*R1* 2>/dev/null -exec gzip -cd {} \; > rawdata/$SAMPLE.R1.fq
  usearch10 -filter_phix rawdata/$SAMPLE.R1.fq -output phix_filtered/$SAMPLE.R1.fq -threads $MAX_THREADS -quiet
  rm rawdata/$SAMPLE.R1.fq

  usearch10 -fastq_filter phix_filtered/$SAMPLE.R1.fq -fastq_maxee 1.0 -fastaout phix_filtered/tempdir/$SAMPLE.R1.QCout.fa -fastq_trunclen 250 -relabel @ -threads $MAX_THREADS -quiet
  cat phix_filtered/tempdir/$SAMPLE.R1.QCout.fa >> all.singlereads.nophix.qc.R1.fa
  rm phix_filtered/tempdir/$SAMPLE.R1.QCout.fa

  # Create concatenated fastq file of nonfiltered reads, with the sample labels
  usearch10 -fastx_relabel phix_filtered/$SAMPLE.R1.fq -prefix $SAMPLE$SAMPLESEP -fastqout phix_filtered/tempdir/$SAMPLE.R1.relabeled.fq -quiet
  cat phix_filtered/tempdir/$SAMPLE.R1.relabeled.fq >> all.singlereads.nophix.R1.fq
  rm phix_filtered/$SAMPLE.R1.fq
done < samples_tmp.txt

echoWithDate "Dereplicating reads..."
usearch10 -fastx_uniques all.singlereads.nophix.qc.R1.fa -sizeout -fastaout uniques.R1.fa -relabel Uniq -quiet

echoWithDate "Generating ASVs (zOTUs) from dereplicated reads..."
usearch10 -unoise3 uniques.R1.fa -zotus zOTUs.R1.fa -quiet

#####PUT AN IF HERE#####
#sed 's/Zotu/ASV/g' zOTUs.R1.fa > ASVs.R1.fa

echoWithDate "Finding all ASVs that match exactly with an already known ASV and renaming accordingly..."
usearch10 -search_exact zOTUs.R1.fa -db $ASVDB -maxaccepts 1 -maxrejects 0 -strand both -dbmatched ASVs.R1.fa -notmatched ASVs_nohits.R1.fa -threads $MAX_THREADS -quiet
usearch10 -fastx_relabel ASVs_nohits.R1.fa -prefix newASV -fastaout ASVs_nohits_renamed.R1.fa -quiet
#combine hits with nohits
cat ASVs_nohits_renamed.R1.fa >> ASVs.R1.fa
#####THAT SHOULD END HERE#####

echoWithDate "Predicting taxonomy of the ASVs..."
usearch10 -sintax ASVs.R1.fa -db "$TAXDB" -tabbedout ASVs.R1.sintax -strand both -sintax_cutoff 0.8 -threads $MAX_THREADS -quiet
sort -V ASVs.R1.sintax -o ASVs.R1.sintax

notused() {
echoWithDate "Adding predicted taxonomy to the ASV FASTA headers..."
R --slave << 'sintaxtofastaheader'
  suppressPackageStartupMessages({
    #Biostrings (and BiocManager which is used to install Biostrings)
    if(!require("Biostrings")) {
      if(!require("BiocManager")) {
        install.packages("BiocManager")
      }
      BiocManager::install("Biostrings", update = FALSE, ask = FALSE)
    }
    if(!require("Biostrings"))
      stop("   The Biostrings R package is not available, skipping step", call. = FALSE)
    
    #dplyr
    if(!require("dplyr")) {
      install.packages("dplyr")
      require("dplyr")
    }
    
    #plyr
    if(!require("plyr")) {
      install.packages("plyr")
      require("plyr")
    }
    
    #stringr
    if(!require("stringr")) {
      install.packages("stringr")
      require("stringr")
    }
  })
  ##### export ASVs with SINTAX taxonomy in headers #####
  sintax <- readLines("ASVs.R1.sorted.sintax")
  taxdf <- plyr::ldply(str_split(sintax, "\t"), function(x) {
    x <- x[c(1,4)]
    if(is.null(x[2]) | x[2] == "")
      x[2] <- "d:unclassified"
    return(x)
  })
  colnames(taxdf) <- c("ASV", "tax")
  
  ASVs.fa <- readBStringSet("ASVs.R1.fa")
  sintaxdf <- left_join(tibble(ASV = names(ASVs.fa)), taxdf, by = "ASV")
  sintax_header <- paste0(sintaxdf[["ASV"]], ";tax=", sintaxdf[["tax"]], ";")
  names(ASVs.fa) <- sintax_header
  writeXStringSet(ASVs.fa, "ASVs.R1.sorted_w_sintax.fa")
sintaxtofastaheader
}

echoWithDate "Generating ASV table..."
usearch10 -otutab all.singlereads.nophix.R1.fq -zotus ASVs.R1.fa -otutabout ASVtable.tsv -threads $MAX_THREADS -sample_delim $SAMPLESEP -quiet
#sort ASVtable
head -n 1 ASVtable.tsv > tmp
tail -n +2 ASVtable.tsv | sort -V >> tmp
mv tmp ASVtable.tsv

echoWithDate "Cleaning up..."
rm -rf rawdata
rm -rf phix_filtered
rm -f samples_tmp.txt
rm -f ASVs_nohits.R1.fa
rm -f ASVs_nohits_renamed.R1.fa

duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echoWithDate "Done in: $duration!"