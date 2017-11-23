INFILE=$1
OUTFILE=$2

converter_path=./bed_dm3_to_bigwig_dm3-057-439-jj-snps.sh
tmppipe=$(mktemp -u)
mkfifo -m 600 "$tmppipe"
bash ${converter_path} ${OUTFILE} <${tmppipe} &
PID=$!
bigWigToBedGraph ${INFILE} ${tmppipe}
wait ${PID}
rm ${tmppipe}

