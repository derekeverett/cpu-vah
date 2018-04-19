rm -R output
mkdir output
rm cpu-vah
make clean
make
./cpu-vah --config rhic-conf/ -o output -h
