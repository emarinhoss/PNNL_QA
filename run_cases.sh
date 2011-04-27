for f in recon_004_MR_*; do
    band=$(echo "$f")
    cd "$band"
    $warpxs -i ssrecon_wv.inp
#    gedit ssrecon_wv.inp &
    cd ../
done
