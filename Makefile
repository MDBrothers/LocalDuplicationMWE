CXX_FLAGS = -std=c++11 -Og
THIRD_PARTY = -I/usr/local/trilinos/include -I/usr/local/Trilinos/include -L/usr/local/trilinos/lib -L/usr/local/Trilinos/lib -lflann -lmpi -lhdf5 -lepetra -lepetraext -lteuchoscore -lteuchoscomm -lteuchosparameterlist -lnox -lnoxepetra -lloca -llocaepetra -I/usr/local/boost/include


all: libplot3d.so libdata.so libperidynamicsProblem.so libproblemInterface.so libLOCAImplicit.so
	mpicxx $(CXX_FLAGS) -o main main.cpp -I/usr/local/trilinos/include -I/usr/local/Trilinos/include -L./ -L/usr/local/trilinos/lib -L/usr/local/Trilinos/lib -lplot3d -ldata -lperidynamicsProblem -lproblemInterface -lLOCAImplicit -lGLU -lGL -lflann -lmpi -lhdf5 -lepetra -lepetraext -lteuchoscore -lteuchoscomm -lteuchosparameterlist -lnox -lnoxepetra -lloca -llocaepetra -I/usr/local/boost/include

compile_tests: 
	node_kernels_test bond_kernels_test neighborhood_kernels_test data_test

node_kernels_test: libnodekernels.so node_kernels.cpp node_kernels_test.cpp kernels.hpp
	mpicxx $(CXX_FLAGS) -o node_kernels_test node_kernels_test.cpp -L./ -lnodekernels

bond_kernels_test: libbondkernels.so bond_kernels.cpp bond_kernels_test.cpp kernels.hpp
	mpicxx $(CXX_FLAGS) -o bond_kernels_test bond_kernels_test.cpp -L./ -lbondkernels

neighborhood_kernels_test: bond_kernels_test libneighborhoodkernels.so neighborhood_kernels.cpp neighborhood_kernels_test.cpp
	mpicxx $(CXX_FLAGS) -o neighborhood_kernels_test neighborhood_kernels_test.cpp -L./ -lbondkernels -lneighborhoodkernels

data_test: libdata.so data_test.cpp data.cpp libplot3d.so
	mpicxx $(CXX_FLAGS) -o data_test data_test.cpp -L./ -ldata -lplot3d $(THIRD_PARTY)

peridynamics_problem_test: libperidynamicsProblem.so libdata.so
	mpicxx $(CXX_FLAGS) -o 

problem_interface_test: libperidynamicsProblem.sp libdata.so libproblemInterface.so
	mpicxx $(CXX_FLAGS) -o problem_interface_test problem_interface_test.cpp -L./ -lperidynamicsProblem -ldata -lproblemInterface

libnodekernels.so: kernels.hpp node_kernels.cpp
	mpicxx $(CXX_FLAGS) -shared -fPIC node_kernels.cpp -o libnodekernels.so

libbondkernels.so: kernels.hpp bond_kernels.cpp
	mpicxx $(CXX_FLAGS) -shared -fPIC bond_kernels.cpp -o libbondkernels.so

libneighborhoodkernels.so: kernels.hpp neighborhood_kernels.cpp 
	mpicxx $(CXX_FLAGS) -shared -fPIC neighborhood_kernels.cpp -o libneighborhoodkernels.so

libplot3d.so:
	mpicxx $(CXX_FLAGS) -shared -fPIC plot3d.cpp -o libplot3d.so -lGLU -lGL

libdata.so: data.cpp data.hpp
	mpicxx $(CXX_FLAGS) -shared -fPIC data.cpp -I/usr/local/Trilinos/include -I/usr/local/trilinos/include -I/usr/local/boost/include -L/usr/local/trilinos/lib -L/usr/local/Trilinos/lib -o libdata.so -lepetra -lepetraext -lteuchoscore -lteuchoscomm -lteuchosparameterlist -lnox

libperidynamicsProblem.so: libdata.so
	mpicxx $(CXX_FLAGS) -shared -fPIC peridynamicsProblem.cpp -I/usr/local/Trilinos/include -I/usr/local/trilinos/include -I/usr/local/boost/include -L/usr/local/trilinos/lib -L/usr/local/Trilinos/lib -L./ -o libperidynamicsProblem.so -ldata -lepetra -lepetraext -lteuchoscore -lteuchoscomm -lteuchosparameterlist -lnox

libproblemInterface.so: libdata.so libperidynamicsProblem.so
	mpicxx $(CXX_FLAGS) -shared -fPIC problemInterface.cpp -I/usr/local/Trilinos/include -I/usr/local/trilinos/include -I/usr/local/boost/include -L/usr/local/trilinos/lib -L/usr/local/Trilinos/lib -L./ -o libproblemInterface.so -ldata -lperidynamicsProblem -lepetra -lepetraext -lteuchoscore -lteuchoscomm -lteuchosparameterlist -lnox -lnoxepetra -lloca

libLOCAImplicit.so: libproblemInterface.so libperidynamicsProblem.so libdata.so
	mpicxx $(CXX_FLAGS) -shared -fPIC LOCAImplicit.cpp -I/usr/local/Trilinos/include -I/usr/local/trilinos/include -I/usr/local/boost/include -L/usr/local/trilinos/lib -L/usr/local/Trilinos/lib -L./ -o libLOCAImplicit.so -lperidynamicsProblem -lproblemInterface -lepetra -lepetraext -lteuchoscore -lteuchoscomm -lteuchosparameterlist -lnox -lnoxepetra -lloca -llocaepetra

clean:
	rm *.so main *_test
