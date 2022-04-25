# The Dll

A shared-link library can also be built as a lower-level interface to core BMDS functionality.

## Mac (homebrew)

Use gcc instead of xcode:

```bash
brew install gcc automake cmake make libtool
ln -s /usr/local/bin/glibtoolize /usr/local/bin/libtoolize
export "CC=/usr/local/bin/gcc-11"
export "CXX=/usr/local/bin/g++-11"
export "LDFLAGS=-L/usr/local/lib/"
export "CPPFLAGS=-I/usr/local/include"
```

Next, install dependencies:

```bash
export "CC=/usr/local/bin/gcc-11"
export "CXX=/usr/local/bin/g++-11"
export "LDFLAGS=-L/usr/local/lib/"
export "CPPFLAGS=-I/usr/local/include"
wget https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
wget https://github.com/stevengj/nlopt/archive/v2.6.2.tar.gz
# build eigen
tar -xvf ./eigen-3.3.7.tar.gz
cd eigen-3.3.7 && mkdir build && cd build && cmake .. && make install
cd ../.. && rm -rf eigen-3.3.7/ && rm eigen-3.3.7.tar.gz
# build nlopt
tar -xvf ./v2.6.2.tar.gz
cd ./nlopt-2.6.2/ && mkdir build && cd build && cmake .. && make install
cd ../..  && rm -rf nlopt-2.6.2/ && rm v2.6.2.tar.gz
```

Finally, build rbmds:

```bash
export "CC=/usr/local/bin/gcc-11"
export "CXX=/usr/local/bin/g++-11"
export "LDFLAGS=-L/usr/local/lib/"
export "CPPFLAGS=-I/usr/local/include"
# build bmds
rm -rf ../BMD_DLL_COMPILE
mkdir ../BMD_DLL_COMPILE

# disable openmp support on mac by default
sed -i 's/.\/configure EIGEN_INCLUDE/.\/configure \-\-disable\-openmp EIGEN_INCLUDE/g' ./create_dll_compile.bash
./create_dll_compile.bash
cd ../BMD_DLL_COMPILE
make -j6  # -j flag is num processors; hardware-specific
mkdir -p ../RBMDS/build/mac
cp .libs/* ../RBMDS/build/mac
```

### Linker troubleshooting

This doesn't seem required in all installs, but just in case:

```bash
# see if anything is missing here....
otool -L bmds/bin/BMDS330/libDRBMD.dylib
bmds/bin/BMDS330/libDRBMD.dylib:
    /usr/local/lib/libDRBMD.0.dylib (compatibility version 1.0.0, current version 1.0.0)
    /usr/local/lib/qt/libnlopt.0.dylib (compatibility version 0.0.0, current version 0.10.0)
    /usr/local/opt/gsl/lib/libgsl.25.dylib (compatibility version 26.0.0, current version 26.0.0)
    /usr/local/opt/gsl/lib/libgslcblas.0.dylib (compatibility version 1.0.0, current version 1.0.0)
    /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1281.100.1)
    /usr/local/lib/gcc/10/libgcc_s.1.dylib (compatibility version 1.0.0, current version 1.0.0)

# sometimes this works, but you can hard-code the path
install_name_tool -change \
    @rpath/libnlopt.0.dylib \
    /usr/local/lib/libnlopt.0.dylib \
    bmds/bin/BMDS330/libDRBMD.dylib
```

## Python debian docker container

To build in a standard debian python container:

```bash
# disable openmp support on mac by default
sed -i 's/.\/configure EIGEN_INCLUDE/.\/configure \-\-disable\-openmp EIGEN_INCLUDE/g' ./create_dll_compile.bash

# build in a fresh container
docker build -t rbmds-python:latest -f dll/python/Dockerfile --platform="linux/x86_64" .

# demo a successful execution inside container
docker run --rm -it rbmds-python:latest python /app/dll/python/dll_test.py

# copy shared libraries for reuse
rm -rf build/linux/python
mkdir -p build/linux/python
chmod 777 build/linux/python
docker run --rm -v $(pwd)/build/linux/python:/tmp rbmds-python:latest sh -c "cp /usr/local/lib/*.so /tmp"

# delete container now that we've captured the libraries
docker rmi rbmds-python:latest
```

The final step copies the required shared libraries into a fresh python container to demo running
as a non-root user w/o the build system. There are a few system-level packages that were required
(see docker file for more details).
