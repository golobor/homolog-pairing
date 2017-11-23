converter_path=./bed_dm3_to_bigwig_dm3-057-439-jj-snps.sh
for path in /net/levsha/share/lab/DrosophilaWulab/tracks/rnaseq/dm3/*.bw
do  
    filename=${path##*/}
    newfilename=${filename/.unique.bw/.dm3-057-439.unique.bw}
    newfilename=${filename/.all.bw/.dm3-057-439.all.bw}
    newpath=/net/levsha/share/lab/DrosophilaWulab/tracks/rnaseq/dm3-057-439/${newfilename}
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
