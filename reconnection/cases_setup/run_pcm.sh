# Preprocess every PCM folder recon_001* and submit it with batch_pcm.msub (Moab).

for f in recon_001*; do
    band=$(echo "$f")

    cd "$band"
    $wxpp -i recon_pcm.pin
    msub batch_pcm.msub
    cd ../

done
