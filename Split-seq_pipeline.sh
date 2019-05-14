#!/usr/bin/env bash
#
# This is a script to process Split-seq reads. It is a modified version of \
# Drop-seq_alignment.sh provided alongside of  Dropseqtools from Steve McCarroll's lab.
#
# Author: Rebekka Wegmann, Snijder lab, ETH Zurich
# Original version Copyright 2017 Broad Institute

# MIT License
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#


tmpdir=`pwd`
outdir=`pwd`
genomedir=
reference=
pipeline=0
echo_prefix=
dropseq_root=$(dirname $0)/external_tools/Drop-seq_tools-2.1.0
star_executable=STAR
estimated_num_cells=500
progname=`basename $0`
splitseq_root=$(dirname $0)
barcode_dir=$splitseq_root/data/barcode_lists

function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <unmapped-queryname-sorted.bam>
Perform Split-seq tagging, barcode filtering, alignment and digital expression matrix calculation

-g <genomedir>      : Directory of STAR genome directory.  Required.
-r <referencefasta> : Reference fasta of the Drop-seq reference metadata bundle.  Required.
-d <dropseq_root>   : Directory containing Drop-seq executables.  Default: Subdirectory of the splitseq toolbox.
-o <outputdir>      : Where to write output bam.  Default: current directory.
-t <tmpdir>         : Where to write temporary files.  Default: current directory.
-s <STAR_path>      : Full path of STAR.  Default: STAR is found via PATH environment variable.
-b <barcode_dir>    : Full path to directory where the list of expected barcodes is stored. Default: subdirectory of the splitseq toolbox. 
-n <num_cells>      : Estimated number of cells in the library. Only affects visualization of barcode filtering results. Default: 500.
-p                  : Reduce file I/O by pipeline commands together.  Requires more memory and processing power.
-e                  : Echo commands instead of executing them.  Cannot use with -p.
EOF
}

function error_exit() {
    echo "ERROR: $1
    " >&2
    usage
    exit 1
}

function check_set() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "$value" ]]
    then error_exit "$name has not been specified.  $flag flag is required"
    fi
}

set -e # exit immediately if any of the commands in the script fail
set -o pipefail # Fail if any of the commands in a pipeline fails

#getopts parses input options. Options followed by a : expect an input argument. The : at the very beginning prevents standard error messages.
while getopts ":d:t:o:pg:r:es:n:b:" options; do
  case $options in
    d ) dropseq_root=$OPTARG;;
    t ) tmpdir=$OPTARG;;
    o ) outdir=$OPTARG;;
    n ) estimated_num_cells=$OPTARG;;
    b ) barcode_dir=$OPTARG;;
    p ) pipeline=1;;
    g ) genomedir=$OPTARG;;
    r ) reference=$OPTARG;;
    s ) star_executable=$OPTARG;;
    e ) echo_prefix="echo";;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $(($OPTIND - 1))

#Checking inputs
if [[ "$pipeline" == 1 && -n "$echo_prefix" ]]
then error_exit "-p and -e cannot be used together"
fi

check_set "$genomedir" "Genome directory" "-g"
check_set "$reference" "Reference fasta"  "-r"

if (( $# != 1 ))
then error_exit "Incorrect number of arguments"
fi

if [[ "$star_executable" != "STAR" ]]
then if [[ ! ( -x $star_executable && -f $star_executable ) ]]
     then error_exit "STAR executable $star_executable passed via -s does not exist or is not executable"
     fi
elif which STAR > /dev/null
then echo > /dev/null
else error_exit "STAR executable must be on the path"
fi

#Create output directories if they do not exist
if [[ ! -d $outdir ]]
then mkdir $outdir
fi

if [[ ! -d $tmpdir ]]
then mkdir $tmpdir
fi

reference_suffix=$(echo $reference | sed s/.*\\./\\./) #reference can be .fa or .fasta
reference_basename=$(basename $reference $reference_suffix)
refflat=$(dirname $reference)/$reference_basename.refFlat
gene_intervals=$(dirname $reference)/$reference_basename.genes.intervals
exon_intervals=$(dirname $reference)/$reference_basename.exon.intervals
rRNA_intervals=$(dirname $reference)/$reference_basename.rRNA.intervals
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar

unmapped_bam=$1
tagged_unmapped_bam=${tmpdir}/unaligned_tagged_BC_filtered.bam
aligned_sam=${tmpdir}/star.Aligned.out.sam
aligned_sorted_bam=${tmpdir}/aligned.sorted.bam
files_to_delete="${aligned_sorted_bam} ${aligned_sam} ${tagged_unmapped_bam}"

## Stage 1: pre-alignment tag
# Extract UMI (bases 1-10 of read2)
tag_molecules="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Molecular.bam_summary.txt \
    BASE_RANGE=1-10 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 INPUT=${unmapped_bam}"

# Extract the 3 cellular barcodes
tag_cells_1="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular_1.bam_summary.txt \
BASE_RANGE=87-94 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XD NUM_BASES_BELOW_QUALITY=1"

tag_cells_2="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular_2.bam_summary.txt \
BASE_RANGE=49-56 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XE NUM_BASES_BELOW_QUALITY=1"

tag_cells_3="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular_3.bam_summary.txt \
BASE_RANGE=11-18 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=true TAG_NAME=XF NUM_BASES_BELOW_QUALITY=1" #setting discard_read=true will make sure read 2 is discarded after the last tagging step, resulting in a tagged, single read bam file

get_primer_type="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_with_primer_type.bam_summary.txt \
BASE_RANGE=95-100 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=true TAG_NAME=XF NUM_BASES_BELOW_QUALITY=3"

#discard all reads where any one of the barcode regions has at least 1 base with quality < 10
filter_bam="${dropseq_root}/FilterBam TAG_REJECT=XQ"

## Stage 2: Trim reads
#trim away adapter (template switching oligo)
trim_starting_sequence="${dropseq_root}/TrimStartingSequence OUTPUT_SUMMARY=${outdir}/adapter_trimming_report.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5"

#trim anything that follows >= 6 contiguous As (assuming this is the polyA tail)
trim_poly_a="${dropseq_root}/PolyATrimmer OUTPUT_SUMMARY=${outdir}/polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6"

## Stage 3: Filter barcodes
# The lower option collapses barcodes that are in the same well, teh upper one trats each barcode as its own. Uncomment the one you want.
# filter_barcodes="python ${splitseq_root}/src/Splitseq_barcode_filtering.py -d ${outdir} -o ${tagged_unmapped_bam} -n ${estimated_num_cells} -b ${barcode_dir}"
filter_barcodes="python ${splitseq_root}/src/Splitseq_barcode_filtering.py -d ${outdir} -o ${tagged_unmapped_bam} -n ${estimated_num_cells} -b ${barcode_dir} --collapse_wells"

# Stage 4: alignment
sam_to_fastq="java -Xmx500m -jar ${picard_jar} SamToFastq INPUT=${tagged_unmapped_bam}"
star_align="$star_executable --genomeDir ${genomedir} --runThreadN 5 --quantMode GeneCounts --outFileNamePrefix ${tmpdir}/star."

# Stage 5: Merge and tag BAM
# sort aligned reads in queryname order (STAR does not necessarily emit reads in the same order as the input)
sort_aligned="java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar ${picard_jar} \
SortSam INPUT=${aligned_sam} OUTPUT=${aligned_sorted_bam} SORT_ORDER=queryname TMP_DIR=${tmpdir}"

# merge and tag aligned reads
merge_bam="java -Xmx4000m -jar ${picard_jar} MergeBamAlignment REFERENCE_SEQUENCE=${reference} UNMAPPED_BAM=${tagged_unmapped_bam} \
ALIGNED_BAM=${aligned_sorted_bam} INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false"

# This is obsolete and was replaced with the two functions below in v.2.0.0
#tag_with_gene_exon="${dropseq_root}/TagReadWithGeneExon O=${outdir}/star_gene_exon_tagged.bam ANNOTATIONS_FILE=${refflat} TAG=GE"

# This one is a more flexible version of ta with gene exon, introduced in version 2.0.0 of Drop-seq tools
tag_with_gene_interval="${dropseq_root}/TagReadWithInterval INTERVALS=${gene_intervals} TAG=XG"
tag_with_gene_function="${dropseq_root}/TagReadWithGeneFunction O=${outdir}/gene_function_tagged.bam ANNOTATIONS_FILE=${refflat}"

if (( $pipeline == 1 ))
then
     # Stage 1 and 2
     $tag_molecules OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \
       $tag_cells_1 INPUT=/dev/stdin OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \
       $tag_cells_2 INPUT=/dev/stdin OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \
       $tag_cells_3 INPUT=/dev/stdin OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \
       # the below is only valid if sequenced with 150 PE, should be commented by default
       $get_primer_type INPUT=/dev/stdin OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \
       $filter_bam INPUT=/dev/stdin OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \
       $trim_starting_sequence INPUT=/dev/stdin OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \
       $trim_poly_a INPUT=/dev/stdin OUTPUT=/dev/stdout | \
       
       
     # Stage 3
     $filter_barcodes

     # Stage 4
     $sam_to_fastq FASTQ=/dev/stdout | \
       $star_align --readFilesIn /dev/stdin

     # Stage 5
     $sort_aligned
     $merge_bam OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \

     # $tag_with_gene_exon I=/dev/stdin
     $tag_with_gene_interval I=/dev/stdin O=/dev/stdout | \
     $tag_with_gene_function I=/dev/stdin
else
     # Stage 1
     $echo_prefix $tag_molecules OUTPUT=$tmpdir/unaligned_tagged_Molecular.bam
     $echo_prefix $tag_cells_1 INPUT=$tmpdir/unaligned_tagged_Molecular.bam OUTPUT=$tmpdir/unaligned_tagged_MC1.bam
     $echo_prefix $tag_cells_2 INPUT=$tmpdir/unaligned_tagged_MC1.bam OUTPUT=$tmpdir/unaligned_tagged_MC1C2.bam
     $echo_prefix $tag_cells_3 INPUT=$tmpdir/unaligned_tagged_MC1C2.bam OUTPUT=$tmpdir/unaligned_tagged_MC1C2C3.bam
     $echo_prefix $filter_bam INPUT=$tmpdir/unaligned_tagged_MC1C2C3.bam OUTPUT=$tmpdir/unaligned_tagged_filtered.bam

     # Stage 2
     $echo_prefix $trim_starting_sequence INPUT=$tmpdir/unaligned_tagged_filtered.bam OUTPUT=$tmpdir/unaligned_tagged_trimmed_smart.bam
     $echo_prefix $trim_poly_a INPUT=$tmpdir/unaligned_tagged_trimmed_smart.bam OUTPUT=$tmpdir/unaligned_mc_tagged_polyA_filtered.bam
     
     files_to_delete="$files_to_delete $tmpdir/unaligned_tagged_Molecular.bam $tmpdir/unaligned_tagged_MC1.bam $tmpdir/unaligned_tagged_MC1C2.bam $tmpdir/unaligned_tagged_MC1C2C3.bam \
                        $tmpdir/unaligned_tagged_filtered.bam $tmpdir/unaligned_tagged_trimmed_smart.bam"
    
     # Stage 3
     $echo_prefix $filter_barcodes -i $tmpdir/unaligned_mc_tagged_polyA_filtered.bam

     # exit 0 #break here, for testing only


     # Stage 4
     $echo_prefix $sam_to_fastq FASTQ=$tmpdir/unaligned_tagged_BC_filtered.fastq
     $echo_prefix $star_align --readFilesIn $tmpdir/unaligned_tagged_BC_filtered.fastq
     files_to_delete="$files_to_delete $tmpdir/unaligned_tagged_BC_filtered.fastq"

     # # Stage 5
     $echo_prefix $sort_aligned
     $echo_prefix $merge_bam OUTPUT=$tmpdir/merged.bam
     # $echo_prefix $tag_with_gene_exon INPUT=$tmpdir/merged.bam
     $echo_prefix $tag_with_gene_interval I=$tmpdir/merged.bam O=$tmpdir/gene_tagged.bam TMP_DIR=${tmpdir}
     $echo_prefix $tag_with_gene_function INPUT=$tmpdir/merged.bam

     files_to_delete="$files_to_delete $tmpdir/merged.bam $tmpdir/gene_tagged.bam"

     $echo_prefix rm $files_to_delete
fi



## Stage 6: create DGE matrix

# counting exonic reads only
dge="${dropseq_root}/DigitalExpression I=${outdir}/gene_function_tagged.bam O=${outdir}/DGE_matrix.txt.gz READ_MQ=10 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100"
$echo_prefix $dge

# counting both intronic and exonic reads
dge_with_introns="${dropseq_root}/DigitalExpression I=${outdir}/gene_function_tagged.bam O=${outdir}/DGE_matrix_with_introns.txt.gz READ_MQ=10 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 LOCUS_FUNCTION_LIST=INTRONIC"
$echo_prefix $dge_with_introns

# collect RNAseq metrics with PICARD
rnaseq_metrics="java -jar ${picard_jar} CollectRnaSeqMetrics I=${outdir}/gene_function_tagged.bam O=${outdir}/rnaseq_metrics.RNA_Metrics REF_FLAT=${refflat} STRAND=FIRST_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=${rRNA_intervals}"
$echo_prefix $rnaseq_metrics

