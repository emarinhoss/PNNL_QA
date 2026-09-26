for f in recon_001_MR_*; do
    band=$(echo "$f")
    cd "$band"
    $warpxs -i ssrecon_wv.inp
    cd ../
done
