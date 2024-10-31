#######
#USAGE:
INPUT=$1
OUTPUT=$2
THREADS=$3
NREADS=$4

mkdir -p $OUTPUT

nreads_per_thread="$((NREADS / THREADS * 2))"

mkdir -p $OUTPUT/split_fa
cat $INPUT | seqtk seq -A | awk -v x="$nreads_per_thread" -v y="$OUTPUT/split_fa" 'BEGIN {n_seq=0;} {if(n_seq%x==0){file=sprintf("%s/split_fasta_%d.fasta",y,n_seq);} print >> file; n_seq++; next;} { print >> file; }'
