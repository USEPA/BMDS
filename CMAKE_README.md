


># CMAKE build instructions


## Developer setup
### Requirements
Python 3.10 or higher
GSL 2.7 or greater
EIGEN 3.4.0 or greater
NLOPT 2.6.2 or greater
CMAKE 3.10 or greater

### Linux compiler
gcc 9 or greater

### Windows compiler
Visual Studio 19 or greater

### Environmental Variables needed
*Note: some may not be needed if installed in standard path*

| Variable Name | Desc | Ex |
|:----------|:-------|:----|
|EIGEN_DIR|path to eigen directory|/home/username/libs/eigen-3.4.0|
|GSL_DIR|path to GSL directory|/usr/local/apps/gsl-2.7|
|NLOPT_DIR|path to nlopt lib directory|/home/username/libs/nlopt-2.6.2/lib64|
|CMAKE_C_COMPILER|path to c compiler|/usr/local/apps/gcc/9.1.0/bin/gcc|
|CMAKE_CXX_COMPILER|path to c++ compiler|/usr/local/apps/gcc/9.1.0/bin/g++|
|PYTHON_EXECUTABLE|path to python executable |/home/username/mambaforge/bin/python|
|PYTHON_LIBRARY_DIR|path to python site-packages directory| /home/username/venv/Lib/site-packages|

## Build process


### clone project
```bash
git clone git@github.com:USEPA/bmds.git
cd bmds
```

### build c++ src
*from base bmds directory*
```bash
cd src
mkdir build
cd build
cmake ..
make   (for Linux/Mac)
cmake --build . --config Release  (for Windows)
```

### build c++ test
*from base bmds directory*
```bash
cd src/tests
mkdir build
cd build
cmake ..
make   (for Linux/Mac)
cmake --build . --config Release  (for Windows)
```

### build python shared object
*from base bmds directory*
```bash
mkdir build
cd build
cmake ..
make   (for Linux/Mac)
cmake --build . --config Release  (for Windows)
```

### build pybmds
follow instructions in [README.md](README.md)


