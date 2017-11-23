for path in /net/levsha/share/lab/DrosophilaWulab/tracks/chip/dm3/*.fc.signal.bw
do  
    filename=${path##*/}
    newfilename=${filename/fc.signal.bw/dm3-057-439.fc.signal.bw}
    newfilename=${newfilename/treat_pileup.bw/dm3-057-439.treat_pileup.bw}
    newfilename=${newfilename/control_lambda.bw/dm3-057-439.control_lambda.bw}
    newpath=/net/levsha/share/lab/DrosophilaWulab/tracks/chip/dm3-057-439/${newfilename}
done

