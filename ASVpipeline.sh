#!/bin/bash
# Set error check, so script doesn't continue if a command fails
set -e
set -o pipefail

usearch=$(which usearch10)
MAX_THREADS=${1:-$((`nproc`-2))}
SEQPATH=/space/sequences/
TAXDB=/space/databases/midas/MiDAS4.6_20191104/output/ESVs_w_sintax.fa
ASVDB=/space/databases/midas/ASVDB_250bp/ASVsV13_250bp_v2.0_20190514/ASVs.R1.fa
prefilterDB=$TAXDB
SAMPLESEP="_"

rm -rf rawdata/
rm -rf phix_filtered/
mkdir -p rawdata
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
    
    #continue only if the sample was actually found and is not empty
    if [ -s "rawdata/$SAMPLE.R1.fq" ]
      then
        #filter PhiX
        $usearch -filter_phix rawdata/$SAMPLE.R1.fq -output phix_filtered/$SAMPLE.R1.fq -threads $MAX_THREADS -quiet
        rm rawdata/$SAMPLE.R1.fq
        
        #QC
        if [ -s "phix_filtered/$SAMPLE.R1.fq" ]
          then
            $usearch -fastq_filter phix_filtered/$SAMPLE.R1.fq -fastq_maxee 1.0 -fastaout phix_filtered/tempdir/$SAMPLE.R1.QCout.fa \
              -fastq_trunclen 250 -relabel @ -threads $MAX_THREADS -quiet
            cat phix_filtered/tempdir/$SAMPLE.R1.QCout.fa >> all.singlereads.nophix.qc.R1.fa
            rm phix_filtered/tempdir/$SAMPLE.R1.QCout.fa
            
            # Create concatenated fastq file of nonfiltered reads, with the sample labels
            $usearch -fastx_relabel phix_filtered/$SAMPLE.R1.fq -prefix $SAMPLE$SAMPLESEP -fastqout phix_filtered/tempdir/$SAMPLE.R1.relabeled.fq -quiet
            cat phix_filtered/tempdir/$SAMPLE.R1.relabeled.fq >> all.singlereads.nophix.R1.fq
            rm phix_filtered/$SAMPLE.R1.fq
        fi
    fi
done < samples_tmp.txt

echoWithDate "Dereplicating reads..."
$usearch -fastx_uniques all.singlereads.nophix.qc.R1.fa -sizeout -fastaout uniques.R1.fa -relabel Uniq -quiet

echoWithDate "Generating ASVs (zOTUs) from dereplicated reads..."
$usearch -unoise3 uniques.R1.fa -zotus zOTUs.R1.fa

echoWithDate "Filtering ASVs that are <60% similar to reference reads..."
if [ -s "$prefilterDB" ]
  then
    $usearch -usearch_global zOTUs.R1.fa -db $prefilterDB \
      -strand both -id 0.6 -maxaccepts 1 -maxrejects 8 -matched prefilt_out.fa -threads $MAX_THREADS -quiet
    mv prefilt_out.fa zOTUs.R1.fa
  else
  	echo "Could not find prefilter reference database, continuing without prefiltering..."
fi

echoWithDate "Searching ASVs against already known ASVs (exact match) and renaming accordingly..."
if [ -s "$ASVDB" ]
  then
    $usearch -search_exact zOTUs.R1.fa -db $ASVDB -maxaccepts 0 -maxrejects 0 -top_hit_only -strand both \
      -dbmatched ASVs.R1.fa -notmatched ASVs_nohits.R1.fa -threads $MAX_THREADS -quiet
    $usearch -fastx_relabel ASVs_nohits.R1.fa -prefix newASV -fastaout ASVs_nohits_renamed.R1.fa -quiet
    #combine hits with nohits
    cat ASVs_nohits_renamed.R1.fa >> ASVs.R1.fa
  else
  	echo "Could not find ASV database, continuing without renaming ASVs..."
    sed 's/Zotu/ASV/g' zOTUs.R1.fa > ASVs.R1.fa
fi

echoWithDate "Predicting taxonomy of the ASVs..."
if [ -s "$TAXDB" ]
  then
    $usearch -sintax ASVs.R1.fa -db "$TAXDB" -tabbedout ASVs.R1.sintax -strand both -sintax_cutoff 0.8 -threads $MAX_THREADS -quiet
    sort -V ASVs.R1.sintax -o ASVs.R1.sintax
  else
    echo "Could not find taxonomy database, continuing without assigning taxonomy..."    
fi

echoWithDate "Generating ASV table..."
$usearch -otutab all.singlereads.nophix.R1.fq -zotus ASVs.R1.fa -otutabout ASVtable.tsv -threads $MAX_THREADS -sample_delim $SAMPLESEP
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
rm -f all.singlereads.nophix.R1.fq
rm -f all.singlereads.nophix.qc.R1.fa

duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echoWithDate "Done in: $duration"
