# Early study: run the solver ($warpxs) in every folder recon_001_MR_*.

for f in recon_001_MR_*; do
    band=$(echo "$f")
    cd "$band"
    $warpxs -i ssrecon_wv.inp
    cd ../
done
