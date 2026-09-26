# Preprocess ($wxpp) and run ($warpxs) every MC folder advect_002_*.

for f in advect_002*; do
    band=$(echo "$f")

    cd "$band"
    $wxpp -i advect_mc.pin
    $warpxs -i advect_mc.inp
    cd ../

done