for f in recon_008_MR_*; do
    band=$(echo "$f")
    cd "$band"
    $wxpp -i ssrecon_wv.pin
    qsub cray.qsub
    cd ../
done
