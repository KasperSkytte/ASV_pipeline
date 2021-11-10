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

#set DK timezone if not set already
if [ -z "$(env | grep '^TZ=')" ]
then
  TZ="Europe/Copenhagen"
fi

#variables
VERSION="1.3.1"
maxthreads=$(($(nproc)-2))
fastq="/space/sequences/Illumina/"
taxdb=""
asvdb=""
prefilterdb="$taxdb"
samplesep="_"
input="samples"
output="$(pwd)/output"
logfilename="asvpipeline_$(date '+%Y%m%d_%H%M%S').txt"
keepfiles="no"
chunksize=5

#default error message if bad usage
usageError() {
  echo "Error: $1" 1>&2
  echo ""
  eval "bash $0 -h"
}

#function to add timestamps to progress messages
scriptMessage() {
  #check user arguments
  if [ ! $# -eq 1 ]
  then
    echo "Error: function must be passed exactly 1 argument" >&2
    exit 1
  fi
  echo " *** [$(date '+%Y-%m-%d %H:%M:%S')] script message: $1"
}

while getopts ":hi:d:f:o:a:kt:v" opt
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
    echo -e "  -k    (flag) Keep all intermediate files as well as all (decompressed) input fastq files, i.e. don't delete ANYTHING."
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
  k )
    keepfiles="yes"
    ;;
  t )
    maxthreads="$OPTARG"
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

#check for GNU parallel
if [ -z "$(which parallel)" ]
then
  usageError "parallel was not found in \$PATH. Please make sure it's installed and findable in \$PATH as exactly: 'parallel'."
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
if [ -n "$(ls -A ${output})" ]
then
  echo "The directory ${output} is not empty, please clear. Exiting..."
  exit 1
fi

logfilepath="${output}/${logfilename}"

#wrap everything in a function to allow writing stderr+stdout to log file
main() {
  echo "#################################################"
  echo "Script: $(realpath "$0")"
  echo "System time: $(date '+%Y-%m-%d %H:%M:%S') (${TZ})"
  echo "Script version: ${VERSION}"
  echo "Current user name: $(whoami)"
  echo "Current working directory: $(pwd)"
  echo "Input file with sample ID's: $(realpath "$input")"
  echo "Output folder: $(realpath -m "$output")"
  echo "fastq folder with Illumina reads: $(realpath "$fastq")"
  if [ -s "$taxdb" ]
  then
    taxdb=$(realpath "${taxdb}")
  else
    taxdb="<not provided>"
  fi
  echo "Taxonomic database: ${taxdb}"
  if [ -s "$asvdb" ]
  then
    asvdb=$(realpath "${asvdb}")
  else
    asvdb="<not provided>"
  fi
  echo "ASV database: ${asvdb}"
  if [ -s "$prefilterdb" ]
  then
    prefilterdb=$(realpath "${prefilterdb}")
  else
    prefilterdb="<not provided>"
  fi
  echo "Prefiltering database: ${prefilterdb}"
  echo "Keep all intermediate/temporary files: ${keepfiles}"
  echo "Max. number of threads: ${maxthreads}"
  echo "Log file: $(realpath -m "$logfilepath")"
  echo "#################################################"
  echo

  tempdir="${output}/temp"
  rawdata="${tempdir}/rawdata"
  phix_filtered="${tempdir}/phix_filtered"
  phix_filtered_temp="${phix_filtered}/tempdir"
  mkdir -p "$output" "$rawdata" "$phix_filtered_temp"

  scriptMessage "Finding samples, filtering PhiX and bad reads, truncating to 250bp..."
  #clean samples file
  # shellcheck disable=SC1003
  cat < "$input" |\
    tr "\r" "\n" |\
    sed -e '$a\' |\
    sed -e '/^$/d' -e 's/ //g' >\
    "${tempdir}/samples.txt"

  nsamples=$(wc -w < "${tempdir}/samples.txt")
  local i=0
  while ((i++)); read -r sample
  do
    echo -ne "Processing sample ($i/$nsamples): $sample\r"
    find "$fastq" \
      -name "${sample}${samplesep}*R1*" 2>/dev/null \
      -exec gzip -cd {} \; >\
       "${rawdata}/${sample}.R1.fq"
    
    #continue only if the sample was actually found and is not empty
    if [ -s "${rawdata}/${sample}.R1.fq" ]
    then
      #filter PhiX
      usearch11 -filter_phix \
        "${rawdata}/${sample}.R1.fq" \
        -output "${phix_filtered}/${sample}.R1.fq" \
        -threads "$maxthreads" \
        -quiet
      if [ "$keepfiles" == "no" ]
      then
        rm "${rawdata}/${sample}.R1.fq"
      fi
      
      #QC
      if [ -s "${phix_filtered}/${sample}.R1.fq" ]
      then
        usearch11 -fastq_filter \
          "${phix_filtered}/${sample}.R1.fq" \
          -fastq_maxee 1.0 \
          -fastaout "${phix_filtered_temp}/${sample}.R1.QCout.fa" \
          -fastq_trunclen 250 \
          -relabel @ \
          -threads "$maxthreads" \
          -quiet
        cat "${phix_filtered_temp}/${sample}.R1.QCout.fa" >>\
          "${tempdir}/all.singlereads.nophix.qc.R1.fa"
        if [ "$keepfiles" == "no" ]
        then
          rm "${phix_filtered_temp}/${sample}.R1.QCout.fa"
        fi
        
        # Create concatenated fastq file of nonfiltered reads, with the sample labels
        usearch11 -fastx_relabel \
          "${phix_filtered}/${sample}.R1.fq" \
          -prefix ${sample}${samplesep} \
          -fastqout "${phix_filtered_temp}/${sample}.R1.relabeled.fq" \
          -quiet
        cat "${phix_filtered_temp}/${sample}.R1.relabeled.fq" >>\
          "${tempdir}/all.singlereads.nophix.R1.fq"
        if [ "$keepfiles" == "no" ]
        then
          rm "${phix_filtered}/${sample}.R1.fq"
        fi
      fi
    else
      echo -e "\n  sample not found or empty file"
    fi
    echo -ne "\e[K"
  done < "${tempdir}/samples.txt"

  scriptMessage "Dereplicating reads..."
  usearch11 -fastx_uniques \
    "${tempdir}/all.singlereads.nophix.qc.R1.fa" \
    -sizeout \
    -fastaout "${tempdir}/uniques.R1.fa" \
    -relabel Uniq \
    -quiet

  scriptMessage "Generating ASVs (zOTUs) from dereplicated reads..."
  usearch11 -unoise3 \
    "${tempdir}/uniques.R1.fa" \
    -zotus "${tempdir}/zOTUs.R1.fa"

  scriptMessage "Filtering ASVs that are <60% similar to reference reads..."
  if [ -s "$prefilterdb" ]
  then
    usearch11 -usearch_global \
      "${tempdir}/zOTUs.R1.fa" \
      -db "$prefilterdb" \
      -strand both \
      -id 0.6 \
      -maxaccepts 1 \
      -maxrejects 8 \
      -matched "${tempdir}/prefilt_out.fa" \
      -threads "$maxthreads" \
      -quiet
    mv "${tempdir}/prefilt_out.fa" "${tempdir}/zOTUs.R1.fa"
  else
    echo "Could not find prefilter reference database, continuing without prefiltering..."
  fi

  scriptMessage "Searching ASVs against already known ASVs (exact match) and renaming accordingly..."
  if [ -s "$asvdb" ]
  then
    usearch11 -search_exact \
      "${tempdir}/zOTUs.R1.fa" \
      -db "$asvdb" \
      -maxaccepts 0 \
      -maxrejects 0 \
      -strand both \
      -dbmatched "${output}/ASVs.R1.fa" \
      -notmatched "${tempdir}/ASVs_nohits.R1.fa" \
      -threads "$maxthreads" \
      -quiet
    usearch11 -fastx_relabel \
      "${tempdir}/ASVs_nohits.R1.fa" \
      -prefix newASV \
      -fastaout "${tempdir}/ASVs_nohits_renamed.R1.fa" \
      -quiet
    #combine hits with nohits
    cat "${tempdir}/ASVs_nohits_renamed.R1.fa" >> "${output}/ASVs.R1.fa"
  else
    echo "Could not find ASV database, continuing without renaming ASVs..."
    sed 's/Zotu/ASV/g' "${tempdir}/zOTUs.R1.fa" > "${output}/ASVs.R1.fa"
  fi

  scriptMessage "Predicting taxonomy of the ASVs..."
  if [ -s "$taxdb" ]
  then
    usearch11 -sintax \
      "${output}/ASVs.R1.fa" \
      -db "$taxdb" \
      -tabbedout "${output}/ASVs.R1.sintax" \
      -strand both \
      -sintax_cutoff 0.8 \
      -threads "$maxthreads" \
      -quiet
    sort -V "${output}/ASVs.R1.sintax" -o "${output}/ASVs.R1.sintax"
  else
    echo "Could not find taxonomy database, continuing without assigning taxonomy..."    
  fi

  scriptMessage "Generating ASV table..."
  #usearch11 -otutab does not scale linearly with the number of threads
  #much faster to split into smaller chunks and run in parallel using
  # GNU parallel and then merge tables afterwards
  jobs=$((( "${maxthreads}" / "${chunksize}" - 1)))
  if [ $jobs -gt 1 ]
  then
    echo "Splitting into $jobs jobs using max $chunksize threads each..."
    splitfolder="${tempdir}/split_asvtable"
    mkdir -p "$splitfolder"

    #split all unfiltered reads
    usearch11 -fastx_split "${tempdir}/all.singlereads.nophix.R1.fq" \
      -splits $jobs \
      -outname "${splitfolder}/all.singlereads.nophix.R1_@" \
      -quiet

    #run a usearch11 -otutab command for each file
    find "$splitfolder" -type f -name '*all.singlereads.nophix.R1_*' |\
      parallel usearch11 -otutab {} \
        -zotus "${output}/ASVs.R1.fa" \
        -otutabout {}_asvtab.tsv \
        -threads $chunksize \
        -sample_delim "_" \
        -quiet

    #generate a comma-separated list of filenames to merge
    asvtabslist=""
    while IFS= read -r -d '' asvtab
    do
      #exclude table if empty, ie only contains one line with "#OTU ID"
      if [ "$(head -n 2 "$asvtab" | wc -l)" -lt 2 ]
      then
        continue
      fi
      if [ -z "$asvtabslist" ]
      then
        asvtabslist="$asvtab"
      else
        asvtabslist="$asvtabslist,$asvtab"
      fi
    done < <(find "$splitfolder" -type f -iname '*_asvtab.tsv' -print0)

    #merge asvtables
    usearch11 -otutab_merge "$asvtabslist" -output "${output}/ASVtable.tsv" -quiet
  else
    #dont run in parallel if maxthreads <= 2*chunksize
    usearch11 -otutab \
      "${tempdir}/all.singlereads.nophix.R1.fq" \
      -zotus "${output}/ASVs.R1.fa" \
      -otutabout "${output}/ASVtable.tsv" \
      -threads "$maxthreads" \
      -sample_delim "$samplesep"
  fi

  #sort ASVtable
  head -n 1 "${output}/ASVtable.tsv" > "${tempdir}/tmp"
  tail -n +2 "${output}/ASVtable.tsv" | sort -V >> "${tempdir}/tmp"
  mv "${tempdir}/tmp" "${output}/ASVtable.tsv"

  if [ "$keepfiles" == "no" ]
  then
    scriptMessage "Cleaning up intermediate files..."
    rm -rf "$tempdir"
  fi

  #print elapsed time since script was invoked
  duration=$(printf '%02dh:%02dm:%02ds\n' $((SECONDS/3600)) $((SECONDS%3600/60)) $((SECONDS%60)))
  scriptMessage "Done. Time elapsed: $duration!"
}

#clear log file first if it happens to already exist
true > "$logfilepath"
main |& tee -a "$logfilepath"
