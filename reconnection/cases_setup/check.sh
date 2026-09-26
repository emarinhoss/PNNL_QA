# Report which PCM folders recon_001* finished (last output recon_pcm_40.h5 exists).

for f in recon_001*; do
    band=$(echo "$f")
    cd "$band"
    if [ -a "recon_pcm_40.h5" ]; then
        pwd
        echo "40 exists"
    else
        pwd
        echo "40 DOES NOT EXITS"
    fi
    cd ../
done
