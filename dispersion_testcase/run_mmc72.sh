for f in advect_mmc_72_*; do
    band=$(echo "$f")
    cd "$band"
    $wxpp -i advect_0.pin
    $warpxser -i advect_0.inp

    if [ -e "advect_1.pin" ]
	then
	$wxpp -i advect_1.pin
	$warpxser -i advect_1.inp
    fi

    if [ -e "advect_2.pin" ]
	then
	$wxpp -i advect_2.pin
	$warpxser -i advect_2.inp
    fi

    cd ../
done
