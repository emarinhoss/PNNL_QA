# Preprocess ($wxpp) and run ($warpxs) every MMC folder advect_001_*,
# including the level-1 and level-2 inputs where they exist.

for f in advect_001*; do
    band=$(echo "$f")
    cd "$band"
    $wxpp -i advect_0.pin
    $warpxs -i advect_0.inp

    if [ -a "advect_1.pin" ]
	then
	$wxpp -i advect_1.pin
	$warpxs -i advect_1.inp
    fi

    if [ -a "advect_2.pin" ]
	then
	$wxpp -i advect_2.pin
	$warpxs -i advect_2.inp
    fi

    cd ../
done
