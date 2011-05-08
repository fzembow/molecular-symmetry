
# get the compounds files
echo "getting files"
wget -O data.sdf.gz ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_000000001_000025000.sdf.gz
wget -O data-3d.sdf.gz ftp://ftp.ncbi.nih.gov/pubchem/Compound_3D/01_conf_per_cmpd/SDF/00000001_00025000.sdf.gz

# uncompress the files
echo "uncompressing files"
gunzip data.sdf.gz
gunzip data-3d.sdf.gz

# run the classifier
echo "classifying files"
python2.6 extract.py data.sdf data-3d.sdf

# clean up
echo "done, cleaning up"
rm data.sdf
rm data-3d.sdf