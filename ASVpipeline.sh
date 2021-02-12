#!/usr/bin/env bash
# This BASH script is based on the template from https://github.com/kasperskytte/bash_template
# License is MIT, which means the author takes no responsibilities, but you can use it for anything you want

#exit when a command fails (use "|| true" to allow a command to fail)
set -o errexit
# exit when a pipe fails
set -o pipefail
#disallow undeclared variables
set -o nounset
#disallow clobbering (overwriting) of files
#set -o noclobber
#print exactly what gets executed (useful for debugging)
#set -o xtrace

#variables
VERSION="1.2.1"
max_threads=$(($(nproc)-2))
fastq="/space/sequences/Illumina/"
taxdb=""
asvdb=""
prefilterdb="$taxdb"
input="samples"
output=$(pwd)

#default error message if bad usage
usageError() {
  echo "Error: $1" 1>&2
  echo ""
  eval "$0 -h"
}

#function to add timestamps to progress messages
scriptMessage() {
  #check user arguments
  if [ ! $# -eq 1 ]
  then
    echo "Error: function must be passed exactly 1 argument" >&2
    exit 1
  fi
  echo " *** [$(date '+%Y-%m-%d %H:%M:%S')] $(basename "$0"): $1"
}

while getopts ":hi:d:f:o:a:t:v" opt
do
case ${opt} in
  h )
    echo "ASV (zOTU) pipeline to process amplicon data. Forward reads only."
    echo "Version: $VERSION"
    echo "Options:"
    echo "  -h    Display this help text and exit."
    echo -e "  -i    (Required) Path to file containing sample ID's to find and process. One sample ID per line. \n          (Default: ${input})"
    echo -e "  -d    Path to taxonomic database. This database will also be used for prefiltering ASV's (<60% id). \n          (Default: ${taxdb})"
    echo -e "  -f    (Required) Path to folder containing fastq files (will be searched recursively). \n          (Default: ${fastq})"
    echo -e "  -o    (Required) Output folder. \n          (Default: ${output})"
    echo -e "  -a    Path to ASV database (a fasta file) with sequences whose names will be reused by exact matching. The ASV's (zOTU's) generated by this script will then be renamed on exact match, any unmatched ASV's will get the prefix 'newASV'. \n          (Default: ${asvdb})"
    echo -e "  -t    Max number of threads to use. \n          (Default: all available except 2)"
    echo "  -v    Print version and exit."
    exit 1
    ;;
  i )
    input="$OPTARG"
    ;;
  d )
    taxdb="$OPTARG"
    if [ ! -s "$taxdb" ]
    then
      usageError "File '${taxdb}' does not exist"
      exit 1
    fi
    ;;
  f )
    fastq="$OPTARG"
    ;;
  o )
    output="$OPTARG"
    ;;
  a )
    asvdb="$OPTARG"
    if [ ! -s "$asvdb" ]
    then
      usageError "File '${asvdb}' does not exist"
      exit 1
    fi
    ;;
  t )
    max_threads="$OPTARG"
    ;;
  v )
    echo "Version: $VERSION"
    exit 0
    ;;
  \? )
    usageError "Invalid Option: -$OPTARG"
    exit 1
    ;;
  : )
    usageError "Option -$OPTARG requires an argument"
    exit 1
    ;;
esac
done
shift $((OPTIND -1)) #reset option pointer

#check for usearch11
if [ -z "$(which usearch11)" ]
then
  usageError "usearch11 was not found in \$PATH. Please make sure it's installed and findable in \$PATH as exactly: 'usearch11'."
  exit 1
fi
# check options
if [ ! -s "$input" ]
then
  usageError "File '${input}' does not exist or is empty"
  exit 1
fi
if [ ! -d "$fastq" ]
then
  usageError "Directory '${fastq}' does not exist"
  exit 1
fi

mkdir -p "$output"

#wrap everything in a function to allow writing stderr+stdout to log file
main() {
  echo "#################################################"
  echo "Script: $(realpath "$0")"
  echo "System time: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "Script version: ${VERSION}"
  echo "Current user name: $(whoami)"
  echo "Current working directory: $(pwd)"
  echo "Input file with sample ID's: $(realpath "$input")"
  echo "Output folder: $(realpath -m "$output")"
  echo "fastq folder with Illumina reads: $(realpath "$fastq")"
  if [ -s "$taxdb" ]
  then
    echo "Taxonomic database: $(realpath "$taxdb")"
  fi
  if [ -s "$asvdb" ]
  then
    echo "ASV database: $(realpath "$asvdb")"
  fi
  if [ -s "$prefilterdb" ]
  then
    echo "Prefiltering database: $(realpath "$prefilterdb")"
  fi
  echo "Max. number of threads: ${max_threads}"
  echo "Log file: $(realpath -m "$logFile")"
  echo "#################################################"
  echo

  tempdir="${output}/temp"
  if [ -d "$tempdir" ]
  then
    scriptMessage "Folder '${tempdir}' already exists, removing it before continuing..."
    rm -rf "$tempdir"
  fi
  mkdir -p "$tempdir"
  
  rawdata="${tempdir}/rawdata"
  mkdir -p "$rawdata"

  #clean samples file
    tr "\r" "\n" < "$input" |\
    sed -e '$a\' |\
    sed -e '/^$/d' -e 's/ //g' > "${tempdir}/samples.txt"

  nsamples=$(wc -w < "${tempdir}/samples.txt")
  scriptMessage "Finding and unpacking ${nsamples} sample(s)..."
  local i=0
  while ((i++)); read -r sample
  do
    echo "($i/$nsamples) $sample"
    sample="${sample}_"
    find "$fastq" \
      -type f \
      -name "*${sample}*R1*.f*q.gz" \
      -exec gzip -cd {} \; > "${rawdata}/${sample}R1.fq" 2> /dev/null
    if [ ! -s "${rawdata}/${sample}R1.fq" ]
    then
      echo " - R1 file not found or empty"
    fi
    find "$fastq" \
      -type f \
      -name "*${sample}*R2*.f*q.gz" \
      -exec gzip -cd {} \; > "${rawdata}/${sample}R2.fq" 2> /dev/null
    if [ ! -s "${rawdata}/${sample}R2.fq" ]
    then
      echo " - R2 file not found or empty"
    fi
  done < "${tempdir}/samples.txt"

  scriptMessage "Merging forward and reverse reads..."
  usearch11 -fastq_mergepairs "${rawdata}"/*_R1*.fq \
    -reverse "${rawdata}"/*_R2*.fq \
    -fastq_maxdiffs 10 \
    -relabel @ \
    -fastqout "${tempdir}"/all_merged.fq \
    -threads "$max_threads" \
    -quiet

  scriptMessage "Filtering reads containing PhiX spike-in..."
  all_merged_nophix="${tempdir}/all_merged_nophix.fq"
  usearch11 -filter_phix "${tempdir}/all_merged.fq" \
    -output "$all_merged_nophix" \
    -threads "$max_threads" \
    -quiet

  scriptMessage "Orienting sequences..."
  if [ -s "$prefilterdb" ]
  then
    usearch11 -orient "$all_merged_nophix" \
      -db "$prefilterdb" \
      -fastqout "${tempdir}/all_merged_nophix_oriented.fq" \
      -threads "$max_threads" \
      -quiet
    mv "${tempdir}/all_merged_nophix_oriented.fq" "$all_merged_nophix"
  else
    scriptMessage " - Could not find prefilter reference database, continuing without orienting"
  fi

  scriptMessage "Filtering reads with expected error > 1.0..."
  usearch11 -fastq_filter \
    "$all_merged_nophix" \
    -fastq_maxee 1.0 \
    -fastaout "${tempdir}/all_filtered.fa" \
    -threads "$max_threads" \
    -quiet

  scriptMessage "Dereplicating reads..."
  usearch11 -fastx_uniques \
    "${tempdir}/all_filtered.fa" \
    -sizeout \
    -fastaout "${tempdir}/all_uniques.fa" \
    -relabel Uniq \
    -quiet

  scriptMessage "Generating ASVs (zOTUs) from dereplicated reads..."
  usearch11 -unoise3 "${tempdir}/all_uniques.fa" \
    -zotus "${tempdir}/zOTUs.fa"

  scriptMessage "Filtering ASVs that are <60% similar to reference reads..."
  if [ -s "$prefilterdb" ]
  then
    usearch11 -usearch_global "${tempdir}/zOTUs.fa" \
      -db "$prefilterdb" \
      -strand both \
      -id 0.6 \
      -maxaccepts 1 \
      -maxrejects 8 \
      -matched "${tempdir}/prefilt_out.fa" \
      -threads "$max_threads" \
      -quiet
    mv "${tempdir}/prefilt_out.fa" "${tempdir}/zOTUs.fa"
  else
    scriptMessage " - Could not find prefilter reference database, continuing without prefiltering..."
  fi

  scriptMessage "Searching ASVs against already known ASVs (exact match) and renaming accordingly..."
  if [ -s "$asvdb" ]
  then
    usearch11 -search_exact "${tempdir}/zOTUs.fa" \
      -db "$asvdb" \
      -maxaccepts 0 \
      -maxrejects 0 \
      -strand both \
      -dbmatched "${output}/ASVs.fa" \
      -notmatched "${tempdir}/ASVs_nohits.fa" \
      -threads "$max_threads" \
      -quiet
    usearch11 -fastx_relabel "${tempdir}/ASVs_nohits.fa" \
      -prefix newASV \
      -fastaout "${tempdir}/ASVs_nohits_renamed.fa" \
      -quiet
    #combine hits with nohits
    cat "${tempdir}/ASVs_nohits_renamed.fa" >> "${output}/ASVs.fa"
  else
    scriptMessage " - Could not find ASV database, continuing without renaming ASVs..."
    sed 's/Zotu/ASV/g' "${tempdir}/zOTUs.fa" > "${output}/ASVs.fa"
  fi

  scriptMessage "Predicting taxonomy of the ASVs..."
  if [ -s "$taxdb" ]
  then
    usearch11 -sintax "${output}/ASVs.fa" \
      -db "$taxdb" \
      -tabbedout "${tempdir}/ASVs.sintax" \
      -strand both \
      -sintax_cutoff 0.8 \
      -threads "$max_threads" \
      -quiet
    sort -V "${tempdir}/ASVs.sintax" -o "${output}/ASVs.sintax"
  else
    scriptMessage " - Could not find taxonomy database, continuing without assigning taxonomy..."    
  fi

  scriptMessage "Generating ASV table..."
  samplesep="$(head -n 1 "$all_merged_nophix" | grep -o '[^0-9][0-9]*$' | grep -o '[^0-9]')"
  echo "Guessed sample separator based on the first read to be '${samplesep}'"
  usearch11 -otutab "$all_merged_nophix" \
    -zotus "${output}/ASVs.fa" \
    -otutabout "${tempdir}/ASVtable.tsv" \
    -threads "$max_threads" \
    -sample_delim "$samplesep"
  #sort ASVtable
  head -n 1 "${tempdir}/ASVtable.tsv" > "${output}/ASVtable.tsv"
  tail -n +2 "${tempdir}/ASVtable.tsv" | sort -V >> "${output}/ASVtable.tsv"

  #print elapsed time since script was invoked
  duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
  scriptMessage "Done in: $duration!"
}

#clear log file first if it happens to already exist
logFile="${output}/log.txt"
true > "$logFile"
main |& tee -a "$logFile"
