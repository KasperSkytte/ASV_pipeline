#!/bin/bash

VERSIONNUMBER=5.0.3.beta1.0

###################################################################################################
#
#  Amplicon DNA workflow
#
#  Version 5.0.3.beta1.0
#
#  This workflow script generates OTU tables from raw bacterial V13 and V4
#  16S rRNA and fungal ITS 1 amplicon data.
#
#  It is currently only supported for internal use at Aalborg University.
#
#  Author: Erika Yashiro, Ph.D.
#
#  Last modified: 2 November, 2018
#
###################################################################################################


# Check if user is running in Bash or Shell
if [ ! -n "$BASH" ]
    then
    echo "Please rerun this script with bash (bash), not shell (sh)."
    echo ""
    echo "Rerun the script with the -h or -help for more information."
    exit 1
fi

# Set error check
set -e
set -o pipefail


#########################################################
# HELP
#########################################################

#if [[ $1 =~ ^(-h|-help)$ ]]  => only works with bash, not shell.
Help_Function () {
# Can move the "set" parameters to after the help function so that initial command can be done with either shell or bash. the set parameters are only usable in bash, not shell.

    
    echo ""
    echo "############################################################################"
    echo "#"
    echo "#  amplicon.workflow.v$VERSIONNUMBER.sh"
    echo "#"
    echo "#  This workflow script generates OTU tables from raw bacterial V13 and V4 "
    echo "#  16S rRNA and fungal ITS 1 amplicon data."
    echo "#"
    echo "#  It is currently only supported for internal use at Aalborg University."
    echo "#"
    echo "#  Erika Yashiro, Ph.D."
    echo "#"
    echo "#  Last modified: 2 November, 2018"
    echo "#"
    echo "############################################################################"
    echo ""
    echo ""
    echo "To run the script's full pipeline: "
    echo "   1. Make sure that you create an empty directory, where you have just the samples file."
    echo "   2. Type in the terminal:      bash amplicon.workflow.v5.0.sh"
    echo "                             or  AmpProc5"   
    echo "   3. Be prepared to answer the questions asked by the script." 
    echo "      Answers are case-sensitive!"
    echo "       - Whether you want OTU and ZOTU tables"
    echo "       - Whether you want single-end and/or paired-end read processing"
    echo "       - Which ribosomal region you have"
    echo "       - Which reference database to use for taxonomy prediction"
    echo ""
    echo "To obtain the README file, Run Step 2. and then type quit. The README file contains the description of all the output files from the workflow, citation of external tools used, and version history."
    echo ""
    echo "To rerun the taxonomy prediction a postiori using a different reference database, run the script with the following arguments."
    echo "  -i    Input file. Must be otus.fa or zotus.fa (or single read variants)."
    echo "  -t    Reference database number for taxonomy prediction."
    echo "              1  - MiDAS v2.1.3"
    echo "              2  - SILVA LTP v128"
    echo "              3  - RDP training set v16"
    echo "              4  - UNITE v7.2 (2017-12-01)"
    echo ""
    echo "To incorporate the new taxonomy output into an OTU table, run the script: otutab_sintax_to_ampvis.v1.1.sh"
    echo ""
    echo "To change the number of CPUs used by the USEARCH stages of the pipeline, adjust the number of threads when running the script."
    echo "The fasttree typically uses about 20 cores at its maximum run and this cannot be adjusted."
    echo ""
    echo "NOTE: The reference databases have been slightly reformatted in order to be compatible with the usearch10 algorithm."
    echo ""
    echo "Have a nice day!"
    echo ""
    echo ""
    exit 0
}


#########################################################
# THREADS
#########################################################

# Define the maximum number of threads (cpus) to use.
NUMTHREADS=5

#########################################################
# OTHER PARAMS
#########################################################

# Define the location of the sequences folders
SEQPATH="/space/sequences/"

# Make /tmp/$USER directory as needed
mkdir -p /tmp/$USER

#########################################################
# FUNCTIONS
#########################################################

Find_reads_phix_Function () {

# Check that samples file exists, if yes, make sure carriage return is not used.
echo ""
echo "Checking the presence of a \"samples\" file"

if [ -f "samples" ]
    then
        # remove carriage returns, add newline at the end of file if it's not there, remove empty lines and white spaces.
        cat samples | tr "\r" "\n" | sed -e '$a\' | sed -e '/^$/d' -e 's/ //g' > samples_tmp.txt
        echo ""
        echo "    Done"
        date
    
    else
        echo ""
        echo "Please make sure that you have a samples list file named \"samples\" in your current directory, and then run this script again."
        echo "You can also run the help function with the option -h or -help."
        echo "    Exiting script."
        date
        echo ""
        exit 1
fi    


# Check that the working directories, or files with the same names, don't exist.
echo ""
echo "Checking the presence of previous samples directories in the current directory"

if [ -f "rawdata" ]
    then
        echo " A file called rawdata already exists. Do you want to remove it and continue with the script? (yes/no)"
        read ANSWER
        if [ $ANSWER = "yes" ]
            then 
            echo "    Removing rawdata file."
            rm -f rawdata
            else 
            if [ $ANSWER = "no" ]
                then
                echo "    Exiting script."
                date
                echo ""
                exit 0
                else
                echo "Sorry I didn't understand you."
                echo "    Exiting script."
                date
                echo ""
                exit 1
            fi
        fi
fi

if [ -d "rawdata" ]
    then
        echo ""
        echo "The rawdata/ directory already exists. Do you want to replace it? (yes/no)"
        read ANSWER
        if [ $ANSWER = "yes" ]
            then 
            echo "    Removing rawdata/ directory."
            rm -rf rawdata/
            else 
            if [ $ANSWER = "no" ]
                then
                echo "    Exiting script."
                date
                echo ""
                exit 0
                else
                echo "Sorry I didn't understand you."
                echo "    Exiting script."
                date
                echo ""
                exit 1
            fi
        fi
fi


if [ -f "phix_filtered" ]
    then
        echo " A file called phix_filtered already exists. Do you want to remove it and continue with the script? (yes/no)"
        read ANSWER
        if [ $ANSWER = "yes" ]
            then 
            echo "    Removing phix_filtered file."
            rm -f phix_filtered
            else 
            if [ $ANSWER = "no" ]
                then
                echo "    Exiting script."
                date
                echo ""
                exit 0
                else
                echo "Sorry I didn't understand you."
                echo "    Exiting script."
                date
                echo ""
                exit 1
            fi
        fi
fi

if [ -d "phix_filtered" ]
    then 
        echo ""
        echo "The phix_filtered/ directory already exists. Do you want to replace it? (yes/no)"
        read ANSWER
        if [ $ANSWER = "yes" ]
            then 
            echo "    Removing phix_filtered/ directory."
            rm -rf phix_filtered/
            else 
            if [ $ANSWER = "no" ]
                then
                echo "    Exiting script."
                date
                echo ""
                exit 0
                else
                echo "Sorry I didn't understand you."
                echo "    Exiting script."
                date
                echo ""
                exit 1
            fi
        fi
fi

echo ""
echo "    Done"
date

echo ""
echo "Retrieving sequenced files and removing PhiX contamination."

# Make new working directories
  mkdir -p rawdata
  mkdir -p phix_filtered

# Find the samples from the samples file
# copy sample sequence files to current directory,
# Filter PhiX
# Path to sequences folders: $SEQPATH = /space/sequences/
  NSAMPLES=$(wc -w < samples_tmp.txt)
  while ((i++)); read SAMPLES
  do
      echo -ne "Processing sample: $SAMPLES ($i / $NSAMPLES)\r"
      # Retrieve sequenced reads
      a="_";
      NAME=$SAMPLES;
      #find /space/sequences/ -name $NAME*R1* 2>/dev/null -exec gzip -cd {} \;
      #find /space/sequences/ -name $NAME*R1* 2>/dev/null -exec cp {} samplegz/ \;
      #find /space/sequences/ -name $NAME* 2>/dev/null -exec cp {} samplegz/ \;
      find $SEQPATH -name $NAME$a*R1* 2>/dev/null -exec gzip -cd {} \; > rawdata/$NAME.R1.fq
      find $SEQPATH -name $NAME$a*R2* 2>/dev/null -exec gzip -cd {} \; > rawdata/$NAME.R2.fq
      # Filter phix
      usearch10 -filter_phix rawdata/$NAME.R1.fq -reverse rawdata/$NAME.R2.fq -output phix_filtered/$NAME.R1.fq -output2 phix_filtered/$NAME.R2.fq -threads $NUMTHREADS -quiet
  done < samples_tmp.txt

#rm -rf rawdata/


    
}


Merge_Function () {

# Merge paired end reads
# Add sample name to read label (-relabel option)
# Pool samples together
# $usearch -fastq_mergepairs ../data/${Sample}*_R1.fq -fastqout $Sample.merged.fq -relabel $Sample.

#while read SAMPLES
#    do
#    NAME=$SAMPLES;
#    #usearch10 -fastq_mergepairs phix_filtered/$NAME.R1.fq -reverse phix_filtered/$NAME.R2.fq -fastqout phix_filtered/$NAME.merged.fq -relabel $NAME -quiet
#    cat phix_filtered/$NAME.merged.fq >> mergeout.fq
#    done < samples_tmp.txt

echo ""
echo "Merging paired end reads"

usearch10 -fastq_mergepairs phix_filtered/*.R1.fq -reverse phix_filtered/*.R2.fq -fastqout mergeout.fq -relabel @ -fastq_maxdiffs 15 -threads $NUMTHREADS -quiet

}

Fastqc_Function () {

# Quality filter
# Note: there is no maxlength in usearch
# V13=425 minlength
# V4=200 minlength
# ITS=200 minlength
# Single reads=250

# INFILE=fastq file of all reads after phix removal
# SEQLEN=minimum length cutoff for all reads.
# output file: QCout.fa
INFILE=$1
SEQLEN=$2

echo ""
echo "Quality filtering and removing consensus reads less than $SEQLEN bp"

usearch10 -fastq_filter $INFILE -fastq_maxee 1.0 -fastaout QCout.fa -fastq_minlen $SEQLEN -quiet -threads $NUMTHREADS

}

Fastqc_singlereads_Function () {

# Quality filter
# Label reads to samples
# Truncate reads to 250bp
# Remove reads less than 250bp

echo ""
echo "Quality filtering, truncating reads to 200bp, and removing reads less than 200bp."

# make temporary directory
mkdir phix_filtered/tempdir

# QC
# merge all sample fastq files
NSAMPLES=$(wc -w < samples_tmp.txt)
while ((j++)); read SAMPLES
    do
    echo -ne "Processing sample: $SAMPLES ($i / $NSAMPLES)\r"
    NAME=$SAMPLES
    usearch10 -fastq_filter phix_filtered/$NAME.R1.fq -fastq_maxee 1.0 -fastaout phix_filtered/tempdir/$NAME.R1.QCout.fa -fastq_trunclen 200 -relabel @ -threads $NUMTHREADS -quiet
    #(-fastq_minlen 250 )
    usearch10 -fastq_filter phix_filtered/$NAME.R2.fq -fastq_maxee 1.0 -fastaout phix_filtered/tempdir/$NAME.R2.QCout.fa -fastq_trunclen 200 -relabel @ -threads $NUMTHREADS -quiet
    cat phix_filtered/tempdir/$NAME.R1.QCout.fa >> all.singlereads.nophix.qc.R1.fa
    cat phix_filtered/tempdir/$NAME.R2.QCout.fa >> all.singlereads.nophix.qc.R2.fa

    # Create concatenated fastq file of nonfiltered reads, with the sample labels
    usearch10 -fastx_relabel phix_filtered/$NAME.R1.fq -prefix $NAME. -fastqout phix_filtered/tempdir/$NAME.R1.relabeled.fq -quiet
    usearch10 -fastx_relabel phix_filtered/$NAME.R2.fq -prefix $NAME. -fastqout phix_filtered/tempdir/$NAME.R2.relabeled.fq -quiet
    cat phix_filtered/tempdir/$NAME.R1.relabeled.fq >> all.singlereads.nophix.R1.fq
    cat phix_filtered/tempdir/$NAME.R2.relabeled.fq >> all.singlereads.nophix.R2.fq
    done < samples_tmp.txt
    
#rm -r phix_filtered/tempdir
    
}



Dereplicate_Function () {

# Find unique read sequences and abundances => Dereplicating
echo ""
echo "Dereplicating reads"

# INFILE=all.merged.nophix.qc.fa and file variants.
# output: DEREPout.fa, which is typically renamed as uniques.fa afterwards.
INFILE=$1
usearch10 -fastx_uniques $INFILE -sizeout -fastaout DEREPout.fa -relabel Uniq -quiet

}



#(EY) Remove primer sequences, use cutadapt? leave blank for now.
    # cutadapt: http://cutadapt.readthedocs.io/en/stable/guide.html
    # other things out there: Bamclipper https://www.nature.com/articles/s41598-017-01703-6
    # edge effects: https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-1073
    # http://www.usadellab.org/cms/?page=trimmomatic
    # https://github.com/ezorita/seeq

# (EY) trim reads to minlen350? => minlen200 and stay generic for v13 and v4


Prefilter_60pc_Function () {

# Prefilter reads <60% ID to reference dataset from full Silva
echo ""
echo "Prefiltering reads that are <60% similar to reference reads"
echo ""

# INFILE=uniques.fa or otus.fa
# output: prefilt_out.fa
INFILE=$1
REF_DATABASE="/space/users/ey/Documents/Amplicon_databases/gg_13_8_otus97/97_otus.fasta"

#usearch10 -closed_ref $INFILE -db /space/users/ey/Documents/gg_13_8_otus97/97_otus.fasta -strand both -id 0.6 -mapout closed_mapped.txt
#usearch10 -usearch_global $INFILE -db $REF_DATABASE -strand both -id 0.6 -maxaccepts 1 -maxrejects 256 -matched prefilt_out.tmp -threads $NUMTHREADS -quiet
usearch10 -usearch_global $INFILE -db $REF_DATABASE -strand both -id 0.6 -maxaccepts 1 -maxrejects 8 -matched prefilt_out.fa -threads $NUMTHREADS

# relabel the unique, prefiltered reads so that the reads are in numerical order.
#usearch10 -fastx_relabel prefilt_out.tmp -prefix Prefilt -fastaout prefilt_out.fa -keep_annots
#usearch10 -sortbysize prefilt_out.tmp -fastaout prefilt_out.fa -quiet

#rm prefilt_out.tmp

echo ""

}


Cluster_otus_Function () {

# Make 97% OTUs and filter chimeras, use de novo
# cluster_otus: OTU clustering with chimera filtering (UPARSE-OTU algorithm)
echo ""
echo "Making 97% OTUs and filter chimeras"

# INFILE=prefilt_out.fa or uniques.fa
# output: otus.fa
INFILE=$1
usearch10 -cluster_otus $INFILE -otus otus.fa -relabel OTU -minsize 2 -quiet

}


Unoise3_Function () {

# Denoise: predict biological sequences and filter chimeras
echo ""
echo "Creating ZOTUs of dereplicated reads file using UNOISE3"

# INFILE=prefilt_out.fa
# output: zotus.fa
INFILE=$1
usearch10 -unoise3 $INFILE -zotus zotus.fa -quiet

}


Make_otutable_Function () {

# Make OTU table
# Ensure that taxonomy file is available
# OTU table without taxonomy information

echo ""
echo "Making an OTU table"

# FASTAFILE=all.merged.nophix.qc.fa (output from Fastqc_singlereads_Function)
# OTUSFILE=otus.fa (output from Cluster_otus_function)
# SINTAX=sintax_out.otus.txt  (output from Predict_taxonomy_Function)
FASTAFILE=$1
OTUSFILE=$2
SINTAX=$3
usearch10 -otutab $FASTAFILE -otus $OTUSFILE -otutabout otutable_notax.txt -id 0.97 -threads $NUMTHREADS -quiet -sample_delim .

bash /space/users/ey/Documents/Scripts/otutab_sintax_to_ampvis.v1.1.sh -i otutable_notax.txt -t $SINTAX -r $REFDATABASE

rm otutable_notax.txt
mv otutable_notax_sorted.txt otutable_notax.txt
mv otutable_$REFDATABASE.txt otutable.txt

}


Make_zotutable_Function () {

# Make ZOTU table
# Ensure that taxonomy file is available

echo ""
echo "Making a zOTU table"
 
# FASTAFILE=all.merged.nophix.qc.fa (output from Fastqc_Function)
# OTUSFILE=zotus.fa (output from Unoise3_Function)
# SINTAX=sintax_out.otus.txt  (output from Predict_taxonomy_Function)
FASTAFILE=$1
ZOTUSFILE=$2
SINTAX=$3
sed 's/Zotu/Otu/g' $ZOTUSFILE > zotus.tmp
usearch10 -otutab $FASTAFILE -zotus zotus.tmp -otutabout zotutable_notax.txt -id 0.97 -threads $NUMTHREADS -quiet -sample_delim . #this needs to be same as with filter_phix+fastq_filter
sed -i 's/Otu/Zotu/g' zotutable_notax.txt
rm zotus.tmp

bash /space/users/ey/Documents/Scripts/otutab_sintax_to_ampvis.v1.1.sh -i zotutable_notax.txt -t $SINTAX -r $REFDATABASE

rm zotutable_notax.txt
mv otutable_notax_sorted.txt zotutable_notax.txt
mv otutable_$REFDATABASE.txt zotutable.txt
#sed -i 's/Otu/Zotu/g' zotutable.txt
#sed -i 's/Otu/Zotu/g' zotutable_notax.txt

}

Predict_taxonomy_Function () {

# Predict taxonomy, set to multithreads on 8 cores
#usearch10 -sintax all.merged.qc.nophix.uniques.otus.fa -db /space/databases/midas/MiDAS_S123_2.1.3.sintax.fasta -strand both -tabbedout all.merged.qc.nophix.uniques.otus.sintax.txt -sintax_cutoff 0.8 -threads 8 ‑notrunclabels 
INFILE=$1
ELEMENT=$2
# INFILE is the otus.fa file
# $ELEMENT is either otus or zotus

echo ""
echo "Predicting taxonomy (Classifying the $ELEMENT) using SINTAX"

# Determine reference database for taxonomy predictions
if [ $REFDATABASE == 1 ]
    then
    REFDATAPATH="/space/users/ey/Documents/Amplicon_databases/midas_database/MiDAS_S123_2.1.3.sintax.cleaned.20180103.udb"
    echo "using MiDAS reference database."
fi

if [ $REFDATABASE == 2 ]
    then
    REFDATAPATH="/space/users/ey/Documents/Amplicon_databases/SILVA_LTP/LTPs128_SSU_unaligned.sintax.udb"
    echo "using SILVA LTP reference database."
fi

if [ $REFDATABASE == 3 ]
    then
    REFDATAPATH="/space/users/ey/Documents/Amplicon_databases/RDP_training_set/rdp_16s_v16s_sp_sintax.cleaned.20180103.udb"
    echo "using RDP reference database."
fi

if [ $REFDATABASE == 4 ]
    then
    REFDATAPATH="/space/users/ey/Documents/Amplicon_databases/UNITE_vers_7.2/utax_reference_dataset_01.12.2017.cleaned.20180103.udb"
    echo "using UNITE reference database."
fi

# Run usearch
    usearch10 -sintax $INFILE -db $REFDATAPATH -strand both -tabbedout sintax_out.txt -sintax_cutoff 0.8 -threads $NUMTHREADS -quiet


}


Taxonomy_reports_Function () {
# Taxonomy summary reports
# INFILE1=sintax_out.txt
# INFILE2=otutable_notax.txt
# OUTFILE=output file prefix

INFILE1=$1
INFILE2=$2
OUTFILE=$3

echo ""
echo "Generating $OUTFILE taxonomy summary"

# if no taxonomy consensus was found for an otu and the 4th column is blank, then sintax_summary causes error, so add something.
sed -i 's/+\t$/+\td:__unknown__/g' $INFILE1
sed -i 's/-\t$/-\td:__unknown__/g' $INFILE1

usearch10 -sintax_summary $INFILE1 -otutabin $INFILE2 -rank p -output $OUTFILE.phylum_summary.txt -quiet
echo ""
echo "    Output phlyum summary: $OUTFILE.phylum_summary.txt"

usearch10 -sintax_summary $INFILE1 -otutabin $INFILE2 -rank c -output $OUTFILE.sintax.class_summary.txt -quiet
echo "    Output class summary: $OUTFILE.sintax.class_summary.txt"

usearch10 -sintax_summary $INFILE1 -otutabin $INFILE2 -rank o -output $OUTFILE.sintax.order_summary.txt -quiet
echo "    Output order summary: $OUTFILE.sintax.order_summary.txt"

usearch10 -sintax_summary $INFILE1 -otutabin $INFILE2 -rank f -output $OUTFILE.sintax.family_summary.txt -quiet
echo "    Output family summary: $OUTFILE.sintax.family_summary.txt"

usearch10 -sintax_summary $INFILE1 -otutabin $INFILE2 -rank g -output $OUTFILE.sintax.genus_summary.txt -quiet
echo "    Output genus summary: $OUTFILE.sintax.genus_summary.txt"


}

Betadiv_Function () {
# Build phylogenetic tree, probably for use in calculating unifrac. Therefore, using fasttree.
# As of December 2017, Current version installed on Dragon is FastTree v2.1.7.
# Need the otus.fa as INFILE,and and otutable_notax.txt as input files.
# export OMP_NUM_THREADS=16
#   This will make the fasttreeMP run only 16 threads instead of all CPUs on the server.

# Note that this kind of error message could occur if there is only 1 otu that passed the zotu filter: /usr/local/lib/python2.7/dist-packages/cogent/maths/unifrac/fast_tree.py:369: RuntimeWarning: invalid value encountered in double_scalars
#  (branch_lengths*logical_or(i,j)).sum())


# INFILE=otus.fa
# OTUTABLE=otutable_notax.txt
# ELEMENT=suffix for type of otu/zotu
INFILE=$1
OTUTABLE=$2
ELEMENT=$3
REP_ALIGNED_PATH="/space/users/ey/Documents/Amplicon_databases/core_set_aligned.fasta.imputed"
USER_PATH=`echo $PWD`

echo ""
echo "Building beta diversity matrices."

# Calculate number of reads per sample
echo ""
echo "   Normalizing the OTU table to 1000"
OTUTABLE2=`echo $OTUTABLE | sed 's/.txt$//g'`
usearch10 -alpha_div $OTUTABLE -output $OTUTABLE2.number_reads_per_sample.txt -metrics reads -quiet

  # Check that at least one sample has at least 1000 reads total.
SAMPLESIZE=`awk -F "\t" 'NR>1{ if ($2 > 1000) {print "OVER1000";exit} }' $OTUTABLE2.number_reads_per_sample.txt`
  # Check how many samples have at least 1000 reads total.
SAMPLENUM=`awk -F "\t" 'NR>1{ if ($2 > 1000) {print "OVER1000"} }' $OTUTABLE2.number_reads_per_sample.txt | wc -l`


if [ "$SAMPLESIZE" = "OVER1000" ]
  then
  # Normalize OTU table to 1000 reads per sample
  usearch10 -otutab_trim $OTUTABLE -min_sample_size 1000 -output $OTUTABLE2.tmp -quiet
  usearch10 -otutab_norm $OTUTABLE2.tmp -sample_size 1000 -output $OTUTABLE2.norm1000.txt -quiet
  rm $OTUTABLE2.tmp
  echo ""
  echo "    Output of normalized OTU table: $OTUTABLE2.norm1000.txt"
  else
  echo ""
  echo "   Cannot normalize OTU table to 1000 reads per sample because none of the samples have >1000 reads." 
  echo "   Using only non-normalized OTU table."
fi

if [[ $AMPREGION =~ ^(V4|V13)$ ]]
    then
    echo ""
    echo "   Aligning the bacterial sequenced reads using PyNAST with QIIME native parameters."
    # Using PyNAST in Unifrac 1.8.0
    align_seqs.py -i $USER_PATH/$INFILE -m pynast -t $REP_ALIGNED_PATH -o $USER_PATH/aligned_seqs_$ELEMENT/

    echo ""
    echo "   Generating FastTree maximum likelihood tree of the bacterial sequenced reads with QIIME native parameters"
    echo ""
    
    # Qiime 1.8 / 1.9 default params is fasttree default params.
    # Set the number of threads for fasttreeMP to 16 cores
    #export OMP_NUM_THREADS=16
    
    #fasttree
    INFILE2=`echo $INFILE | sed 's/.fa$//g'`
    fasttreeMP -nt aligned_seqs_$ELEMENT/${INFILE2}_aligned.fasta > aligned_seqs_$ELEMENT/$INFILE.$ELEMENT.tre
    
    # Reset the OMP threads
    #export OMP_NUM_THREADS=""
    
    # or: make_phylogeny.py -i $USER_PATH/aligned_seqs/${INFILE2}_aligned.fasta -o $USER_PATH/$INFILE.tre

    echo ""
    echo "   Warning: The R package Ampvis uses the Generalized UniFrac instead of the original weighted and unweighted UniFrac equations implemented in QIIME version 1.x.x."
    echo ""
    echo "   Generating beta diversity matrices: Bray Curtis, original version of weighted & unweighted UniFrac from Fasttree tree"

    # Create betadiv folder
    mkdir beta_div_$ELEMENT

    # Convert classic otu table to biom format
    biom convert -i $OTUTABLE -o $OTUTABLE.biom --table-type="OTU table" --to-hdf5

    # Run Qiime 1.8.0 beta_diversity script for UniFrac
    beta_diversity.py -i $OTUTABLE.biom -m weighted_unifrac,unweighted_unifrac -o beta_div_$ELEMENT/ -t aligned_seqs_$ELEMENT/$INFILE.$ELEMENT.tre
    # Change file names of output matrices
    mv beta_div_$ELEMENT/weighted_unifrac_$OTUTABLE.txt beta_div_$ELEMENT/$ELEMENT.weighted_unifrac.txt
    mv beta_div_$ELEMENT/unweighted_unifrac_$OTUTABLE.txt beta_div_$ELEMENT/$ELEMENT.unweighted_unifrac.txt

    # Run Usearch for Bray Curtis
    usearch10 -beta_div $OTUTABLE -metrics bray_curtis -filename_prefix beta_div_$ELEMENT/$ELEMENT. -quiet
    
   if [ "$SAMPLESIZE" = "OVER1000" ] && [ "$SAMPLENUM" -gt 1 ]
      then
      # Convert normalized otu table to biom format
      biom convert -i $OTUTABLE2.norm1000.txt -o $OTUTABLE2.norm1000.biom --table-type="OTU table" --to-hdf5

      # Run Qiime script for UniFrac matrices
      beta_diversity.py -i $OTUTABLE2.norm1000.biom -m weighted_unifrac,unweighted_unifrac -o beta_div_norm1000_$ELEMENT/ -t aligned_seqs_$ELEMENT/$INFILE.$ELEMENT.tre
      # Change file name of output matrices
      mv beta_div_norm1000_$ELEMENT/weighted_unifrac_$OTUTABLE2.norm1000.txt beta_div_norm1000_$ELEMENT/$ELEMENT.weighted_unifrac.txt
      mv beta_div_norm1000_$ELEMENT/unweighted_unifrac_$OTUTABLE2.norm1000.txt beta_div_norm1000_$ELEMENT/$ELEMENT.unweighted_unifrac.txt

      # Run Usearch for Bray Curtis matrix
      usearch10 -beta_div $OTUTABLE2.norm1000.txt -metrics bray_curtis -filename_prefix beta_div_norm1000_$ELEMENT/$ELEMENT. -quiet
      else
      echo ""
      echo "   Note: Beta diversity matrices from normalized OTU table could not be generated."
   fi
    
    
    echo ""
    echo "    Output files of alignment in aligned_seqs_$ELEMENT/"
    echo "    Output files of beta diversity in beta_div_$ELEMENT/"

   if [ "$SAMPLESIZE" = "OVER1000" ] && [ "$SAMPLENUM" -gt 1 ]
      then
    echo "    Output files of beta diversity of normalized data in beta_div_norm1000_$ELEMENT/"
   fi
fi

if [ $AMPREGION = "ITS" ]
    then

    # Check that there is more than 2 OTUs in the OTU table
    NUMOTUS=`grep -c ">" $INFILE`

    if [ $NUMOTUS -le 1 ]
      then
      echo ""
      echo "   There is not enough $ELEMENT to generate beta diversity output. No matrices generated"
      else
       # Build clustering tree for fungi
       echo ""
       echo "  Using maximum linkage clustering to build OTU tree"
       mkdir beta_div_$ELEMENT

       usearch10 -cluster_agg $INFILE -treeout beta_div_$ELEMENT/$ELEMENT.cluster.tre -id 0.80 -linkage max -quiet
    
       echo ""
       echo "   Warning: Fungal ITS regions are too variable for proper phylogenetic tree. Therefore the maximum linkage tree will be used for generating phlyogeny-based beta diversity matrices."
       echo ""
       echo "   Generating beta diversity matrices: Bray Curtis, weighted & unweighted UniFrac from clustering tree"

       
       # Calculate beta diversity matrices
       usearch10 -beta_div $OTUTABLE -metrics bray_curtis,unifrac,unifrac_binary -tree beta_div_$ELEMENT/$ELEMENT.cluster.tre -filename_prefix beta_div_$ELEMENT/$ELEMENT. -quiet
    
      if [ $SAMPLESIZE = "OVER1000" ] && [ "$SAMPLENUM" -gt 1 ]
         then
         # Calculate beta diveristy matrices for normalized otu table of normalized otu table is large enough.
         usearch10 -beta_div $OTUTABLE2.norm1000.txt -metrics bray_curtis,unifrac,unifrac_binary -tree beta_div_$ELEMENT/$ELEMENT.cluster.tre -filename_prefix beta_div_$ELEMENT/$ELEMENT.norm1000_ -quiet
         else
         echo ""
         echo "   Note: Beta diversity matrices from normalized OTU table could not be generated."
      fi
        
       echo ""
       echo "    Output files of clustering: $ELEMENT.cluster.tre"
       echo "    Output files of beta diversity in beta_div_$ELEMENT/"
   fi
fi

}


Cleanup_Function () {

rm -r rawdata
rm -r phix_filtered
rm samples_tmp.txt
rm *.nophix.*
rm -f prefilt_out.*
rm sintax_out.*
rm -f *.biom
rm -f uniques.*
rm -f otus.R1.tmp
rm -f otus.R2.tmp
rm -f otus.tmp

mkdir taxonomy_summary
mv *_summary.txt taxonomy_summary/.

}

#########################################################
# ARGUMENTS
#########################################################


# Arguments: help/h, presence of samples file, or go with a postiori taxonomy options.
if [[ "$1" =~ ^(-help|-h)$ ]]
    then
    Help_Function
    else
    # if there are arguments present
    if [ $1 ]
        then
        echo ""
        echo "Running: 16S workflow version $VERSIONNUMBER"
        echo "A postiori taxonomy assignment"
        echo "To incorporate new taxonomy into OTU table, run otutab_sintax_to_ampvis.v1.1.sh"
        date
        while getopts :i:t: option
            do
            case "${option}"
            in
            i) OTUINFILE=${OPTARG} ;;
            t) TAX=${OPTARG} ;;
            #o) OTUOUTFILE=${OPTARG} ;;
            #h) Help_Function ;;
            \?) echo ""
                echo "Invalid option: -$OPTARG" >&2 
                echo "Check the help manual using -help or -h"
                echo "Exiting script."
                echo ""
                exit 1 ;;
            :) echo ""
               echo "Option -$OPTARG requires an argument"
               echo "Check the help manual using -help or -h"
               echo "Exiting script."
               echo ""
               exit 1 ;;
            esac
        done
        
        if [ ! -f $OTUINFILE ]
            then 
            echo ""
            echo "Input file $OTUINFILE does not exist. Exiting script."
            echo ""
            exit 1
        fi
        
       # if [[ $TAX !=~ ^(1|2|3|4)$ ]]
	if [ $TAX -lt 1 ] || [ $TAX -gt 4 ]
            then
            echo ""
            echo "Taxonomy needs to be selected (1,2,3, or 4). Check -help or -h for more information. Exiting script."
            echo ""
            exit 1
        fi
        
        REFDATABASE=$TAX
        if [ $TAX == 1 ]
         then
         TAXFILE="midas"
        fi
        
        if [ $TAX == 2 ]
         then
         TAXFILE="silva"
        fi
        
        if [ $TAX == 3 ]
         then
         TAXFILE="rdp"
        fi
        
        if [ $TAX == 4 ]
         then
         TAXFILE="unite"
        fi
        # Run taxonomy prediction with specified reference database
        Predict_taxonomy_Function $OTUINFILE OTUs
        # input file radical, remove file extension
        FILERAD=`echo $OTUINFILE | sed -e 's/\.fa$//g' -e 's/\.fas$//g' -e 's/\.fasta$//g'`
        mv sintax_out.txt $FILERAD.sintax.$TAXFILE.txt
        echo "Output files of the OTU/ZOTU taxonomy assignment: $FILERAD.sintax.$TAXFILE.txt"
        date
        echo ""
        exit 0
        
        else
        if [ ! -f "samples" ]
        then 
        echo "samples file does not exist. Check -help or -h for more information. Exiting script."  
        echo ""
        exit 1
	fi
    fi
fi
    

#########################################################
# MAIN WORKFLOW
#########################################################

# Define ZOTUS = yes/no
# Define single read = yes/no/both (SINGLEREADS)
# Define amplicon region = V13/V4/ITS (AMPREGION)
# Define reference database = 1/2/3/4 (REFDATABASE)


#########################################################
# Questions
#########################################################

clear
echo ""
echo "Running: 16S workflow version $VERSIONNUMBER"
echo ""
echo "WARNING: Please note that this version uses a different workflow than v.4.3, so you are advised to rerun this script on your older datasets if you wish to proceed with comparative analyses with the older data."
echo ""
echo "Use the -h or -help options to run the Help function"
echo ""
echo "Read the README file for more information on how to cite this workflow."
echo "Description of all the output files are also found in the README file."
date
echo ""

# Copy the README file
cp /space/users/ey/Documents/Scripts/amplicon_workflow/README_amplicon.workflow_v$VERSIONNUMBER.txt .

# Define ZOTUs = yes/no
echo ""
echo "In addition to an OTU table, do you want to generate a ZOTU table using UNOISE3? (yes/no/quit)"
echo "        yes  - Generate ZOTU table and OTU table"
echo "        no   - Generate only OTU table"
echo "        quit - Quit the script so I can go and read up on UNOISE3"
read ZOTUS

# Check that the question answers are script readable.
# note: arg between quotes means it can be nothing without producing error.
if [[ ! "$ZOTUS" =~ ^(yes|no|quit)$ ]]
    then
    echo ""
    echo "OTU / ZOTU table: $ZOTUS is an invalid argument."
    echo "Sorry I didn't understand what you wrote. Please also make sure that you have the correct upper/lower case letters."
    echo "You can also run the help function with the option -help or -h."
    echo "    Exiting script"
    date
    echo ""
    exit 1
    else
    if [ $ZOTUS = "quit" ]
        then
        echo ""
        echo "    Exiting script"
        date
        echo ""
        exit 0
    fi
fi


# Define single read process = SR/PE/both (SINGLEREADS)
echo ""
echo "Do you want process single-end reads or paired-end reads? (yes/no/both)"
echo "Select to run single end reads if you want to control for species that have very long variable regions between the forward and reverse primers. These species can be rare in the community but will get filtered out if the paired end reads do not overlap."
echo "        SR  - Process only single reads"
echo "        PE   - Process only paired-end reads"
echo "        both - process both single reads and paired-end reads"
read SINGLEREADS

# Define amplicon region = V13/V4/ITS (AMPREGION)
echo ""
echo "What genomic region does your PCR amplicons amplify? (V13/V4/ITS)"
echo "        V13 - Bacterial 16S rRNA hypervariable regions 1 & 3"
echo "        V4  - Bacterial 16S rRNA hypervariable region 4"
echo "        ITS - Fungal ribosomal ITS 1 region"
read AMPREGION

# Define reference database for taxonomy prediction of OTUs (REFDATABASE)
echo ""
echo "Which reference database do you want to use for taxonomy prediction? ()"
echo "        1  - MiDAS v2.1.3"
echo "        2  - SILVA LTP v128"
echo "        3  - RDP training set v16"
echo "        4  - UNITE v7.2 (2017-12-01) [Select for Fungal analysis!]"
read REFDATABASE

# Define number of threads to use (NUMTHREADS)
echo ""
echo "How many CPUs do you want to use at maximum run?"
echo "Dragon (up to 80), Coco and Peanut (up to 12)"
echo "If unsure, then start the run with 5 threads"
read NUMTHREADS


# Check that the question answers are script readable.
# note: arg between quotes means it can be nothing without producing error.

if [[ ! "$SINGLEREADS" =~ ^(SR|both|PE)$ ]]
    then
    echo ""
    echo "Single reads / SR+PE: $SINGLEREADS invalid argument."
    echo "Sorry I didn't understand what you wrote. Please also make sure that you have the correct upper/lower case letters."
    echo "You can also run the help function with the option -help or -h."
    echo "    Exiting script"
    date
    echo ""
    exit 1
fi

if [[ ! "$AMPREGION" =~ ^(V13|V4|ITS)$ ]]
    then
    echo ""
    echo "Amplicon region: $AMPREGION invalid argument."
    echo "Sorry I didn't understand what you wrote. Please also make sure that you have the correct upper/lower case letters."
    echo "You can also run the help function with the option -help or -h."
    echo "    Exiting script"
    date
    echo ""
    exit 1
fi

#if [[ ! "$REFDATABASE" =~ ^(1|2|3|4)$ ]]
if [[ "$REFDATABASE" -lt 1 ]] || [[ "$REFDATABASE" -gt 4 ]]
    then
    echo ""
    echo "Reference database: $REFDATABASE invalid argument."
    echo "Sorry I didn't understand what you wrote. Please make sure that you select the correct database."
    echo "You can also run the help function with the option -help or -h."
    echo "    Exiting script"
    date
    echo ""
    exit 1
fi



#########################################################
# DATA PROCESSING
#########################################################


# Finding your samples and copying them to the current directory
# Filter potential phiX contamination
echo ""
echo "Finding your samples and copying them to the current directory, and filtering potential phiX contamination"
Find_reads_phix_Function
echo "     Output raw sequences in phix_filtered/"
date


################
# SINGLE READS
################
  
if [[ $SINGLEREADS =~ ^(SR|both)$ ]]
    then
    # Quality filtering
    Fastqc_singlereads_Function
    echo "    Output files of QC:" 
    echo "      all.singlereads.nophix.qc.R1.fa"
    echo "      all.singlereads.nophix.qc.R2.fa"
    date

    # Dereplicate uniques
    Dereplicate_Function all.singlereads.nophix.qc.R1.fa
    mv DEREPout.fa uniques.R1.fa
    Dereplicate_Function all.singlereads.nophix.qc.R2.fa
    mv DEREPout.fa uniques.R2.fa
    echo "    Output files of derep: uniques.R1.fa and uniques.R2.fa"
    date
  
#############          
 #   # Prefilter to remove anomalous reads
 #   # Prefiltering at 60% ID not implemented for ITS.
 #   if [ $AMPREGION = "ITS" ]
 #       then
 #           mv uniques.R1.fa prefilt_out.R1.fa
 #           mv uniques.R2.fa prefilt_out.R2.fa
 #       else
 #           # Prefilter bacteria to remove anomalous reads
 #           Prefilter_60pc_Function uniques.R1.fa
 #           mv prefilt_out.fa prefilt_out.R1.fa
 #           Prefilter_60pc_Function uniques.R2.fa
 #           mv prefilt_out.fa prefilt_out.R2.fa
 #           echo "    Output files of prefiltering: prefilt_out.R1.fa and prefilt_out.R2.fa"
 #   fi
 #############
           
    # Cluster OTUs + taxonomy assignment
 #   Cluster_otus_Function prefilt_out.R1.fa
    Cluster_otus_Function uniques.R1.fa
    mv otus.fa otus.R1.tmp
#    Cluster_otus_Function prefilt_out.R2.fa
    Cluster_otus_Function uniques.R2.fa
    mv otus.fa otus.R2.tmp

   # Prefilter to remove anomalous reads
   # Prefiltering at 60% ID not implemented for ITS.
   if [ $AMPREGION = "ITS" ]
    then
      mv otus.R1.tmp otus.R1.fa
      mv otus.R2.tmp otus.R2.fa
    else
      Prefilter_60pc_Function otus.R1.tmp
      mv prefilt_out.fa otus.R1.fa
      Prefilter_60pc_Function otus.R2.tmp
      mv prefilt_out.fa otus.R2.fa
   fi
    echo "    Output files of OTU clustering: otus.R1.fa and otus.R2.fa"
           
    Predict_taxonomy_Function otus.R1.fa OTUsR1
    mv sintax_out.txt sintax_out.otus.R1.txt
    Predict_taxonomy_Function otus.R2.fa OTUsR2
    mv sintax_out.txt sintax_out.otus.R2.txt
    echo "    Output files of OTU taxonomy assignment: sintax_out.otus.R1.txt and sintax_out.otus.R2.txt"
    
    if [ $ZOTUS = "yes" ]
        then
        # Cluster zOTUs + taxonomy assignment
#        Unoise3_Function prefilt_out.R1.fa
        Unoise3_Function uniques.R1.fa
        mv zotus.fa zotus.R1.tmp
#        Unoise3_Function prefilt_out.R2.fa
        Unoise3_Function uniques.R2.fa
        mv zotus.fa zotus.R2.tmp

     # Prefilter to remove anomalous reads
     # Prefiltering at 60% ID not implemented for ITS.
     if [ $AMPREGION = "ITS" ]
      then
        mv zotus.R1.tmp zotus.R1.fa
        mv zotus.R2.tmp zotus.R2.fa
      else
        Prefilter_60pc_Function zotus.R1.tmp
        mv prefilt_out.fa zotus.R1.fa
        Prefilter_60pc_Function zotus.R2.tmp
        mv prefilt_out.fa zotus.R2.fa
        rm zotus.R1.tmp zotus.R2.tmp
     fi
        echo "    Output ZOTU files: zotus.R1.fa and zotus.R2.fa"
        date
                
        Predict_taxonomy_Function zotus.R1.fa ZOTUsR1
        mv sintax_out.txt sintax_out.zotus.R1.txt
        Predict_taxonomy_Function zotus.R2.fa ZOTUsR2
        mv sintax_out.txt sintax_out.zotus.R2.txt
        echo "    Output ZOTU taxonomy files: sintax_out.zotus.R1.txt and sintax_out.zotus.R2.txt"
        date
    fi
            
    # Build OTU table
    Make_otutable_Function all.singlereads.nophix.R1.fq otus.R1.fa sintax_out.otus.R1.txt
    mv otutable.txt otutable.R1.txt
    mv otutable_notax.txt otutable_notax.R1.txt
    Make_otutable_Function all.singlereads.nophix.R2.fq otus.R2.fa sintax_out.otus.R2.txt
    mv otutable.txt otutable.R2.txt
    mv otutable_notax.txt otutable_notax.R2.txt
    echo "    Output file of OTU table: otutable_notax.R1.txt and otutable_notax.R1.txt"
    echo "    Output file for final otu table: otutable.R1.txt and otutable.R2.txt"
    date

    # Build zOTU table
    if [ $ZOTUS = "yes" ]
        then
        Make_zotutable_Function all.singlereads.nophix.R1.fq zotus.R1.fa sintax_out.zotus.R1.txt
        mv zotutable.txt zotutable.R1.txt
        mv zotutable_notax.txt zotutable_notax.R1.txt
        Make_zotutable_Function all.singlereads.nophix.R2.fq zotus.R2.fa sintax_out.zotus.R2.txt
        mv zotutable.txt zotutable.R2.txt
        mv zotutable_notax.txt zotutable_notax.R2.txt
    
        echo "    Output file of zOTU table: zotutable_notax.txt"
        echo "    Output file for final zOTU table: zotutable.txt"
        date
    fi

    # Taxonomy reports from SINTAX output
    Taxonomy_reports_Function sintax_out.otus.R1.txt otutable_notax.R1.txt otusR1
    Taxonomy_reports_Function sintax_out.otus.R2.txt otutable_notax.R2.txt otusR2

    if [ $ZOTUS = "yes" ]
        then
        Taxonomy_reports_Function sintax_out.zotus.R1.txt zotutable_notax.R1.txt zotusR1
        Taxonomy_reports_Function sintax_out.zotus.R2.txt zotutable_notax.R2.txt zotusR2
    fi
fi

if [ $SINGLEREADS = "SR" ]
    then
    echo ""
    echo "Removing temporary files and directories"
    Cleanup_Function
    echo ""
    echo "Single read data processing is done. Enjoy."
    date
    echo ""
    exit 0
    
    else
    if [ $SINGLEREADS = "both" ]
      then
      echo ""
      echo "Single read data processing is done."
      date
    fi
fi


################
# PE READS
################

echo ""
echo "Starting workflow for paired-end read data."

# Merge paired ends
Merge_Function
mv mergeout.fq all.merged.nophix.fq
echo "    Output file of merging: all.merged.nophix.fq"
date

# Quality filtering
if [ $AMPREGION = "V13" ]
    then
    Fastqc_Function all.merged.nophix.fq 425
fi

if [ $AMPREGION = "V4" ]
    then
    Fastqc_Function all.merged.nophix.fq 200
fi

if [ $AMPREGION = "ITS" ]
    then
    Fastqc_Function all.merged.nophix.fq 200
fi

mv QCout.fa all.merged.nophix.qc.fa
echo "    Output file of QC: all.merged.nophix.qc.fa"
date

# Dereplicate to uniques
Dereplicate_Function all.merged.nophix.qc.fa
mv DEREPout.fa uniques.fa
echo "    Output file of derep: uniques.fa"
date

########################
## Prefiltering at 60% ID not implemented for ITS.
## Prefiltering at 60% ID for V4 and V13.
#if [ $AMPREGION = "ITS" ]
#    then
#        mv uniques.fa prefilt_out.fa
#    else
#        # Prefilter to remove anomalous reads
#        Prefilter_60pc_Function uniques.fa
#        #mv prefilt_out.fa all.merged.nophix.qc.uniques.prefilt.fa
#        echo "    Output file of prefiltering: prefilt_out.fa"
#fi
#########################

# Cluster OTUs + taxonomy assignment
#Cluster_otus_Function prefilt_out.fa
Cluster_otus_Function uniques.fa
mv otus.fa otus.tmp

if [ $AMPREGION = "ITS" ]
  then
    mv otus.tmp otus.fa
  else
    # Prefilter to remove any other anomalous reads
    Prefilter_60pc_Function otus.tmp
    mv prefilt_out.fa otus.fa
fi
echo "    Output file of OTU clustering: otus.fa"
date

Predict_taxonomy_Function otus.fa OTUs
mv sintax_out.txt sintax_out.otus.txt
echo "    Output file of taxonomy: sintax_out.otus.txt"
date

# Build OTU table
Make_otutable_Function all.merged.nophix.fq otus.fa sintax_out.otus.txt
echo "    Output file of OTU table: otutable_notax.txt"
echo "    Output file for final otu table: otutable.txt"
date

# Taxonomy reports from SINTAX output
Taxonomy_reports_Function sintax_out.otus.txt otutable_notax.txt otus
date

# Align sequences, build tree, and generate beta diversity matrices that can be fed into Ampvis or R base.
Betadiv_Function otus.fa otutable_notax.txt OTUS

#echo ""
#echo "    Output files of alignment in aligned_seqs_OTUS/"
#echo "    Output files of beta diversity in beta_div_OTUS/"
date

# Cluster zOTUs + taxonomy assignment
if [ $ZOTUS = "yes" ]
    then
#    Unoise3_Function prefilt_out.fa 
#    #mv zotus.fa zotus.fa
     Unoise3_Function uniques.fa
     mv zotus.fa zotus.tmp

     if [ $AMPREGION = "ITS" ]
       then
         mv zotus.tmp zotus.fa
       else
         # Prefilter to remove any other anomalous reads
         Prefilter_60pc_Function zotus.tmp
         mv prefilt_out.fa zotus.fa
     fi

    echo "    Output zotu file: zotus.fa"
    date

    Predict_taxonomy_Function zotus.fa ZOTUs
    mv sintax_out.txt sintax_out.zotus.txt
    echo "    Output file of ZOTU taxonomy files: sintax_out.zotus.txt"
    date

    # Build zOTU table
    Make_zotutable_Function all.merged.nophix.fq zotus.fa sintax_out.zotus.txt
    echo "    Output file of zOTU table: zotutable_notax.txt"
    echo "    Output file for final zOTU table: zotutable.txt"
    date

    # Taxonomy reports from SINTAX output
    Taxonomy_reports_Function sintax_out.zotus.txt zotutable_notax.txt zotus
    date

    # Align sequences, build tree, and generate beta diversity matrices that can be fed into Ampvis or R base.
    Betadiv_Function zotus.fa zotutable_notax.txt ZOTUS
    date
fi




#########################################################
# TEMPORARY FILE REMOVAL
#########################################################

# Remove temporary files and directories
echo ""
echo "Removing temporary files and directories"

Cleanup_Function

echo ""
echo "Paired-end read processing is done. Enjoy."
date
echo ""


