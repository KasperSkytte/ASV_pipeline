#!/bin/bash
# Set error check
set -e
set -o pipefail

NUMTHREADS=${1:-$((`nproc`-2))}
SEQPATH="/home/kapper/dragun/space/sequences/"
TAXDB="/home/kapper/Dropbox/AAU/PhD/Projects/ASV_test/midas30_20180105.fa"
SAMPLESEP="_"

rm -rf rawdata/
rm -rf phix_filtered/
mkdir -p rawdata
mkdir -p phix_filtered

echoWithDate() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')]: $1"
}

echoWithDate "Finding samples and filtering PhiX from reads (forward only)..."
cat samples | tr "\r" "\n" | sed -e '$a\' | sed -e '/^$/d' -e 's/ //g' > samples_tmp.txt

NSAMPLES=$(wc -w < samples_tmp.txt)
while ((i++)); read SAMPLE
  do
  echo -ne "Processing sample $SAMPLE ($i / $NSAMPLES)\r"
  find $SEQPATH -name $SAMPLE$SAMPLESEP*R1* 2>/dev/null -exec gzip -cd {} \; > rawdata/$SAMPLE.R1.fq
  usearch10 -filter_phix rawdata/$SAMPLE.R1.fq -output phix_filtered/$SAMPLE.R1.fq -threads $NUMTHREADS -quiet
  rm rawdata/$SAMPLE.R1.fq
done < samples_tmp.txt

echoWithDate "Quality filtering, truncating reads to 250bp..."
mkdir phix_filtered/tempdir
# merge all sample fastq files
while ((j++)); read SAMPLE
  do
  echo -ne "Processing sample $SAMPLE ($j / $NSAMPLES)\r"
  usearch10 -fastq_filter phix_filtered/$SAMPLE.R1.fq -fastq_maxee 1.0 -fastaout phix_filtered/tempdir/$SAMPLE.R1.QCout.fa -fastq_trunclen 250 -relabel @ -threads $NUMTHREADS -quiet
  cat phix_filtered/tempdir/$SAMPLE.R1.QCout.fa >> all.singlereads.nophix.qc.R1.fa

  # Create concatenated fastq file of nonfiltered reads, with the sample labels
  usearch10 -fastx_relabel phix_filtered/$SAMPLE.R1.fq -prefix $SAMPLE$SAMPLESEP -fastqout phix_filtered/tempdir/$SAMPLE.R1.relabeled.fq -quiet
  cat phix_filtered/tempdir/$SAMPLE.R1.relabeled.fq >> all.singlereads.nophix.R1.fq
done < samples_tmp.txt

echoWithDate "Dereplicating reads..."
usearch10 -fastx_uniques all.singlereads.nophix.qc.R1.fa -sizeout -fastaout uniques.R1.fa -relabel Uniq -quiet

echoWithDate "Generating zOTUs of dereplicated reads..."
usearch10 -unoise3 uniques.R1.fa -zotus ASVs.R1.fa -quiet
sed -i 's/Zotu/ASV/g' ASVs.R1.fa

echoWithDate "Predicting taxonomy the zOTUs..."
usearch10 -sintax ASVs.R1.fa -db $TAXDB -tabbedout ASVs.R1.sintax -strand both -sintax_cutoff 0.8 -threads $NUMTHREADS -quiet
sort -V ASVs.R1.sintax > ASVs.R1.sorted.sintax

echoWithDate "Generating zOTU table..."
usearch10 -otutab all.singlereads.nophix.R1.fq -zotus ASVs.R1.fa -otutabout zotutable_notax.R1.txt -mapout ASVmapping.txt -threads $NUMTHREADS -sample_delim $SAMPLESEP
#sort table
head -n 1 zotutable_notax.R1.txt > tmp & tail -n +2 zotutable_notax.R1.txt | sort -V >> tmp & mv tmp zotutable_notax.R1.txt

echoWithDate "Cleaning up..."
rm -rf rawdata
rm -rf phix_filtered
rm -f samples_tmp.txt
rm -f ASVs.R1.sintax

duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echoWithDate "Done in: $duration!"