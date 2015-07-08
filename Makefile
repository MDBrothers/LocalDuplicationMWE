CXX_FLAGS = -std=c++11 -g
THIRD_PARTY = -I/usr/local/hdf5/include -I/usr/local/trilinos/include -I/usr/local/Trilinos/include -L/usr/local/hdf5/lib -L/usr/local/trilinos/lib -L/usr/local/Trilinos/lib -lflann -lmpi -lhdf5 -lepetra -lepetraext -lteuchoscore -lteuchoscomm -lteuchosparameterlist -lnox -lnoxepetra -lloca -llocaepetra -I/usr/local/boost/include

compile_tests: 
	data_test

map_test: just_the_map.cpp
	mpicxx $(CXX_FLAGS) -o map_test just_the_map.cpp $(THIRD_PARTY)

hash_test: hashes_for_enums.cpp
	mpicxx $(CXX_FLAGS) -o hash_test hashes_for_enums.cpp $(THIRD_PARTY)

data_test: libdata.so data_test.cpp data.cpp libplot3d.so
	mpicxx $(CXX_FLAGS) -o data_test data_test.cpp -L./ -ldata -lplot3d $(THIRD_PARTY)

libneighborhoodkernels.so: kernels.hpp neighborhood_kernels.cpp 
	mpicxx $(CXX_FLAGS) -shared -fPIC neighborhood_kernels.cpp -o libneighborhoodkernels.so

libplot3d.so:
	mpicxx $(CXX_FLAGS) -shared -fPIC plot3d.cpp -o libplot3d.so -lGLU -lGL

libdata.so: data.cpp data.hpp
	mpicxx $(CXX_FLAGS) -shared -fPIC data.cpp -o libdata.so $(THIRD_PARTY)

clean:
	rm *.so main *_test
