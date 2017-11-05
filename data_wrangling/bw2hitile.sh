parallel  --ungroup -P3 'clodius aggregate bigwig -a dm3-057-439 {}' ::: *.bw
parallel  -P1  'hgserver ingest -n higlass-jjwu2 --type hitile --display-name "ChIP__{= s/.*\///; s/\..*// =}, fc, dm3-057-439" {}' ::: tracks/chip/dm3-057-439/*fc.signal.hitile
