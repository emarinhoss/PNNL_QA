for f in advect_003*; do
    band=$(echo "$f")

    cd "$band"
    $wxpp -i advect_pcm.pin
    $warpxs -i advect_pcm.inp
    cd ../

done