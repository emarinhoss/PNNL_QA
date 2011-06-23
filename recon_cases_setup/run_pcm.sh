for f in recon_001*; do
    band=$(echo "$f")

    cd "$band"
    $wxpp -i recon_pcm.pin
    msub batch_pcm.msub
    cd ../

done
