# Preprocess every MC folder recon_008_MR_* and submit it with cray.qsub.

for f in recon_008_MR_*; do
    band=$(echo "$f")
    cd "$band"
    $wxpp -i ssrecon_wv.pin
    qsub cray.qsub
    cd ../
done
