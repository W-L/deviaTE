#! /bin/bash

deviaTE=$(dirname $0)


# parse arguments with getopts --------
# defaults
LIB="$deviaTE/library/te_library"
Q_TR=15
MIN_RL=1
MIN_AL=1
Q_ENC='sanger'
THREADS=1

usage(){
    printf $(basename $0)" usage:\n"
    grep "^    .)  #" $0 | tr -d ')#' | sed 's/^..../&-/'
    exit 1
}

# check arguments starting with -
while getopts ":i:l:q:r:a:y:t:h" opt; do
    case $opt in
    i)  # input fastq-file; required
        INP=$OPTARG;;
    l)  # reference library; default provided
        LIB=$OPTARG;;
    q)  # quality threshold for trimming; 15
        Q_TR=$OPTARG;;
    r)  # min read length for trimming; 40
        MIN_RL=$OPTARG;;
    a)  # min length for high-scoring pairs; 1
        MIN_AL=$OPTARG;;
    y)  # quality encoding; either 'sanger' or 'illumina'; 'sanger'
        Q_ENC=$OPTARG;;
    t)  # num. of threads; 1
        THREADS=$OPTARG;;
    \? | h) 
        usage;;
    :)  printf "Option -$OPTARG requires an argument." >&2
        exit 1;;
    esac
done

# if no input is given
if [ -z "$INP" ]
then
    usage
fi

# set location of scripts/progs
TRIMMING="$deviaTE/prog/popoolation/basic-pipeline/trim-fastq.pl"
FILTER_LENGTH="$deviaTE/bin/deviaTE_filter_al.py"
FUSING="$deviaTE/bin/deviaTE_fuse.py"
MAPPING=$(source $deviaTE/CONFIG | tail -n1)
SAMTOOLS=$(source $deviaTE/CONFIG | head -n1)

# set tmp filenames
LOG=$INP.log
TRIM=$INP.trim.tmp
MAP=$INP.sam.tmp
FILT=$INP.filt.tmp
FUSE=$INP.fused
#SORT=$INP.sort.bam
#INDEX=$INP.sort.bam.bai


:>$LOG
printf "starting to prepare $INP \n\n" >>$LOG

# trimming
trim="perl $TRIMMING --input1 $INP --output $TRIM --quality-threshold $Q_TR --min-length $MIN_RL --fastq-type $Q_ENC 2>&1 | tee -a $LOG"
printf "$trim\n\n" >>$LOG
eval $trim

# count total read length for normalization
trl=$(grep -A1 ^@ $TRIM | egrep -v '^@|-' | tr -d '\n\r' | wc -m)
printf "#total_read_length: $trl\n\n" >>$LOG

# index reference
index="$MAPPING index $LIB 2>&1 | tee -a $LOG"
printf "$index\n\n" >>$LOG
eval $index

# mapping
mapping="$MAPPING bwasw -t $THREADS $LIB $INP >$MAP 2>>$LOG"
printf "\n$mapping\n\n" >>$LOG
eval $mapping

# filtering unmapped etc.
filtering="$SAMTOOLS view -h $MAP -o tmp -@ $THREADS -F 0x4 -F 0x400 -F 0x200 2>&1 | tee -a $LOG"
printf "\n$filtering\n\n" >>$LOG
eval $filtering

# filtering length
if [ $MIN_AL != 1 ]; then
    len_filt="python3 $FILTER_LENGTH $MIN_AL $FILT 2>&1 | tee -a $LOG"
    printf "$len_filt\n\n" >>$LOG
    eval $len_filt
else
    mv="mv tmp $FILT 2>&1 | tee -a $LOG"
    printf "$mv\n\n" >>$LOG
    eval $mv
fi

# fuse split reads
fuse="python3 $FUSING --input $FILT --output $FUSE 2>&1 | tee -a $LOG"
printf "$fuse\n\n" >>$LOG
eval $fuse

## sort and index bam
#sort="$SAMTOOLS view -b $FUSE -@ $THREADS | $SAMTOOLS sort -o $SORT -@ $THREADS 2>&1 | tee -a $LOG"
#printf "$sort\n\n" >>$LOG
#eval $sort
#
#indx="$SAMTOOLS index $SORT $INDEX -@ $THREADS 2>&1 | tee -a $LOG"
#printf "$indx\n\n" >>$LOG
#eval $indx

# remove tmp files
rm $INP.*.tmp*
printf "prep done\n\n" >>$LOG
