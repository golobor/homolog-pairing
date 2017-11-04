converter_path=./bed_dm3_to_bigwig_dm3-057-439-jj-snps.sh
for path in /net/levsha/share/lab/DrosophilaWulab/tracks/chip/dm3/*.fc.signal.bw
do  
    filename=${path##*/}
    newfilename=${filename/fc.signal.bw/dm3-057-439.fc.signal.bw}
    newfilename=${newfilename/treat_pileup.bw/dm3-057-439.treat_pileup.bw}
    newfilename=${newfilename/control_lambda.bw/dm3-057-439.control_lambda.bw}
    newpath=/net/levsha/share/lab/DrosophilaWulab/tracks/chip/dm3-057-439/${newfilename}
    if [ ! -f ${newpath} ]; then
        tmppipe=$(mktemp -u)
        mkfifo -m 600 "$tmppipe"
        echo "convert ${path} to ${newpath}"
        bash ${converter_path} ${newpath} <${tmppipe} &
        PID=$!
        bigWigToBedGraph ${path} ${tmppipe}
        wait ${PID}
        rm ${tmppipe}
    fi
done
