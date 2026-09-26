# Preprocess every MMC folder recon_006_MR_* and submit one PBS job per level
# (cray_0/1/2.qsub).

for f in recon_006_MR_*; do
    band=$(echo "$f")
    cd "$band"
    $wxpp -i ssrecon_wv_0.pin
    qsub cray_0.qsub

    if [ -a "ssrecon_wv_1.pin" ]
	then
	$wxpp -i ssrecon_wv_1.pin
	qsub cray_1.qsub
    fi

    if [ -a "ssrecon_wv_2.pin" ]
	then
	$wxpp -i ssrecon_wv_2.pin
	qsub cray_2.qsub
    fi

    cd ../
done
