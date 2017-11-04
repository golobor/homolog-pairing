ADDCHR=false
ISMALE=false
while getopts "cm" flag; do
case "$flag" in
    c) ADDCHR=true;;
    m) ISMALE=true;;
esac
done

shift $(expr $OPTIND - 1 )

if [ "$#" -ne 1 ]; then
    echo "USAGE: cat INPUT | bed_dm3_to_bigwig_dm3-057-439-jj-snps.sh OUTFILE"
    exit 0
fi

OUTFILE=$1

TMPDIR=`mktemp -d`
trap "rm -rf $TMPDIR" EXIT

TMPCHROMSIZESFILE=${TMPDIR}/dm3-057-439-jj-snps-male.reduced.chromsizes
echo "chr2L_057	23011544" >> ${TMPCHROMSIZESFILE}
echo "chr2R_057	21146708" >> ${TMPCHROMSIZESFILE}
echo "chr2L_439	23011544" >> ${TMPCHROMSIZESFILE}
echo "chr2R_439	21146708" >> ${TMPCHROMSIZESFILE}
echo "chr3L_057	24543557" >> ${TMPCHROMSIZESFILE}
echo "chr3R_057	27905053" >> ${TMPCHROMSIZESFILE}
echo "chr3L_439	24543557" >> ${TMPCHROMSIZESFILE}
echo "chr3R_439	27905053" >> ${TMPCHROMSIZESFILE}
echo "chrX_057	22422827" >> ${TMPCHROMSIZESFILE}
if [ "$ISMALE" = false ]; then
    echo "chrX_439	22422827" >> ${TMPCHROMSIZESFILE}
fi

TMP057BEDFILE=${TMPDIR}/057_sorted_filtered_addedchr.bed.lz4
TMP439BEDFILE=${TMPDIR}/439_sorted_filtered_addedchr.bed.lz4
TMPSORTEDBEDFILE=${TMPDIR}/sorted_filtered_addedchr.bed
TMPPIPE=${TMPDIR}/057pipe
mkfifo ${TMPPIPE}
#TMP057BEDFILE=/tmp/057_filtered_addedchr.bed.lz4
#TMP439BEDFILE=/tmp/439_filtered_addedchr.bed.lz4
#TMPSORTEDBEDFILE=/tmp/sorted_filtered_addedchr.bed

if [ "$ADDCHR" = true ]; then
    { sed 's/^\([^\t]*\)\t/chr\1_057\t/' | lz4c -cz ; } >${TMP057BEDFILE} <${TMPPIPE} &
    PID=$!
    cat | grep -P '^2L\t|^2R\t|^3L\t|^3R\t|^X\t' | tee ${TMPPIPE} \
        | { if [ "$ISMALE" = true ]; then grep -v '^X'; else cat; fi } \
        | sed 's/^\([^\t]*\)\t/chr\1_439\t/' | lz4c -cz >${TMP439BEDFILE} 
    wait ${PID}
else
    { sed 's/^\([^\t]*\)\t/\1_057\t/' | lz4c -cz ; } >${TMP057BEDFILE} <${TMPPIPE} &
    PID=$!
    cat | grep -P '^chr2L\t|^chr2R\t|^chr3L\t|^chr3R\t|^chrX\t' | tee ${TMPPIPE} \
        | { if [ "$ISMALE" = true ]; then grep -v '^chrX'; else cat; fi } \
        | sed 's/^\([^\t]*\)\t/\1_439\t/' | lz4c -cz >${TMP439BEDFILE} 
    wait ${PID}
fi

{ lz4c -dc ${TMP057BEDFILE} ; lz4c -dc ${TMP439BEDFILE} ; } | sort -k1,1 -k2,2n > ${TMPSORTEDBEDFILE} 

#| cut -f1,2,3,5 \

bedGraphToBigWig ${TMPSORTEDBEDFILE} ${TMPCHROMSIZESFILE} ${OUTFILE} 

