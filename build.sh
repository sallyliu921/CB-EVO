cd lib/abc
make -j128 
cd ../..

mkdir build
cd build
cmake ..
make -j128
ln -s ../lib/abc/abc ./abc 
cp ../lib/Cgen/ccirc/ccirc ./