# Report which MC folders recon_008_MR_* finished (last output ssrecon_wv_40.h5 exists).

for f in recon_008_MR_*; do
    band=$(echo "$f")
    cd "$band"
    if [ -a "ssrecon_wv_40.h5" ]; then
	pwd
	echo "ended"
	cd ../
    else
	pwd
	echo "didn't end"
	cd ../
	#rm -rf "$band"
    fi
done
