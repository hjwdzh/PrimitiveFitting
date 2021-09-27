mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
cd ..
./bin/Release/region_growing_3d
./bin/Release/region_growing_2d -c data/color.jpg -d data/depth.pgm -o mask.png

