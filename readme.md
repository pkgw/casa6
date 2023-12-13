# CASA 

This repository contains the sources of the [CASA](https://casa.nrao.edu/) package.

CASA, the Common Astronomy Software Applications package, is the primary data processing software for the Atacama Large Millimeter/submillimeter Array [ALMA](https://www.almaobservatory.org/en/home/) and NSF's Karl G. Jansky Very Large Array [VLA](https://public.nrao.edu/visit/very-large-array/), and is frequently used also for other radio telescopes.

The CASA software can process data from both single-dish and aperture-synthesis telescopes, and one of its core functionalities is to support the data reduction and imaging pipelines for ALMA, VLA and the VLA Sky Survey [VLASS](https://public.nrao.edu/vlass/).

A full description of all the functionality is descried in the [User documentation](https://casadocs.readthedocs.io/en/stable/).

## User installation

CASA is available as Python wheels from a [PyPI repository](https://casa-pip.nrao.edu/).
These wheels are compatible with Python 3.8 and Python 3.10 in Linux and macOS platforms.

The recommended way to install CASA is using a virtual environment:
```
$: python3 -m venv myvenv
$: source myvenv/bin/activate
(myvenv) $: pip install --upgrade pip wheel
(myvenv) $: pip install casatasks
```

This will install the following Python packages:

  * [casatools](casatools/readme.md) -- the C++ tools with a minimal Python layer
  * [casatasks](casatasks/readme.md) -- a pure Python layer which provides a higher level of abstraction

However, users are free to install CASA's wheels as they like.

## Source code organization

[casatasks](casatasks/readme.md) contains the high level components that the users use for their data processing.
[casatools](casatools/readme.md) offers a lower level Python interface that is used mainly by the casatasks but can also be used by experienced users.
[casacpp](casatools/src/code) implements the low level C++ algorithms that are used by casatools.

## Developer builds

Here is a summary of the process:

1. Install the software prerequisites. Depending on your platform this might be as easy as install system packages or can be more elaborated.
1. Install casacore and libsakura
1. Install CASA C++ (aka casacpp).
1. Install casatools
1. Install casatasks

There is a Makefile in the casa distribution that simplifies all but the first step. See section "Installation using Makefile"

In the instructions below the following variables are being used as placeholders for different directories:

- $CASAINSTALL has been chosen as a placeholder for the directory where the compiled software will be installed. Define the variable or replace it accordingly.
- $CASASRC has been chosen as a placeholder for the directory which contains a git clone of the CASA6 repository.
- $CASABUILD has been chosen as a placeholder for the directory where the build is performed. This can be deleted after the whole installation has finished.
- $CASATESTDIR has been chosen as a placeholder for a directory that will be used to run the casa tests.

What follows is a detailed explanation of these steps:

### Install the software  prerequisites

The following software is needed to install casa, casacore and libsakura:

- cmake
- C++ compiler (in the case of gcc at least version 4.9 and optionally support for OpenMP)
- Fortran compiler
- OpenMPI (optional)
- readline library
- ncurses library
- blas library
- lapack library
- fftw3 library
- wcslib library
- libxml library
- libxslt library
- gsl library
- sqlite library
- protobuf library + protobuf compiler
- gRPC + gRPC protobuf plugin
- python3. Supported versions: 3.8 and 3.10.
- numpy
- flex command
- bison command
- pkg-config command
- curl command
- tar command

Requirements for GPU build (development still in progress):

- C++ compiler (gcc 9+)
- CUDA
- cuFFT
- ATLAS
- Kokkos (libkokkosscore)

To install these dependencies one can either install them manually using any native method (which might not be trivial) or use package managers. The later is way easier and we have collected the instructions for the following operating systems:

- RHEL equivalent (tested on RockyLinux 8\)
- Fedora (tested on 36)
- Ubuntu (tested on 22.04)
- Debian (tested on 11)
- Manylinux2014/RHEL 7 pkg deps
- macOS (tested on 11 and 12)

What follows is the detailed instructions to install the prerequisites for all these platforms.

#### Installing prerequisites in RHEL equivalent (tested on RockyLinux 8 only)

Note: on RHEL8.5 powertools is provided by: 

    codeready-builder-for-rhel-8-x86_64-rpms                                 Red Hat CodeReady Linux Builder for RHEL 8 x86_64 (RPMs)

Note: on NRAO systems CodeReadyBuilder is packaged as part of nrao-rhel-8.repo and EPEL is packaged as nrao-epel-8.repo and both should be installed on all RHEL8 NRAO system.

Run as root or as a used with sudo rights the following commands:

```
    # Make sure that the needed repos are there
    $ dnf -y install epel-release
    $ dnf install -y dnf-plugins-core
    $ dnf config-manager --set-enabled powertools

    # Packages needed for casacore development 
    $ sudo dnf -y install git cmake gcc-c++ gcc-gfortran gtest-devel ccache readline-devel ncurses-devel blas-devel lapack-devel cfitsio-devel fftw-devel wcslib-devel python38 python38-devel python38-numpy flex bison tar curl

    # Packages needed for libsakura development
    $ sudo dnf -y install eigen3-devel fftw-devel

    # Additional packages needed for casa development 
    $ dnf -y install java-1.?.0-openjdk-devel chrpath python38-wheel python38-numpy swig pkgconf-pkg-config xerces-c-devel libxml2-devel libxslt-devel gsl-devel sqlite-devel wcslib-devel openmpi-devel xorg-x11-server-Xvfb
 
    # Set python 3.8 as the default python
    $ sudo alternatives --set python /usr/bin/python3.8
    $ sudo alternatives --set python3 /usr/bin/python3.8

    $ dnf -y install https://casa.nrao.edu/download/devel/grpc/el8-1.18/protobuf-3.6.1-3.el8.x86_64.rpm
    $ dnf -y install https://casa.nrao.edu/download/devel/grpc/el8-1.18/protobuf-compiler-3.6.1-3.el8.x86_64.rpm
    $ dnf -y install https://casa.nrao.edu/download/devel/grpc/el8-1.18/protobuf-devel-3.6.1-3.el8.x86_64.rpm
    $ dnf -y install https://casa.nrao.edu/download/devel/grpc/el8-1.18/grpc-1.18.0-2.el8.x86_64.rpm
    $ dnf -y install https://casa.nrao.edu/download/devel/grpc/el8-1.18/grpc-devel-1.18.0-2.el8.x86_64.rpm
    $ dnf -y install https://casa.nrao.edu/download/devel/grpc/el8-1.18/grpc-plugins-1.18.0-2.el8.x86_64.rpm
```

#### Installing prerequisites in Fedora (tested only in Fedora 36)

Run as root or as a used with sudo rights the following commands:
```
    # Packages needed for casacore development 
    $ sudo dnf -y install git cmake gcc-c++ gcc-gfortran gtest-devel ccache readline-devel ncurses-devel blas-devel lapack-devel cfitsio-devel fftw-devel wcslib-devel python3 python3-devel python3-numpy python-unversioned-command flex bison tar curl

    # Packages needed for libsakura development
    $ sudo dnf -y install eigen3-devel fftw-devel

    # Additional packages needed for casa development 
    $ sudo dnf -y install java-*-openjdk-devel chrpath perl-File-Fetch python3-scipy python3-matplotlib python3-certifi python3-pytest-xvfb python3-pytest python3-wheel python3-build python3-numpy swig pkgconf-pkg-config xerces-c-devel libxml2-devel libxslt-devel gsl-devel sqlite-devel protobuf-compiler grpc grpc-devel grpc-plugins wcslib-devel openmpi-devel xorg-x11-server-Xvfb redhat-lsb-core patchelf ImageMagick
```
#### Installing prerequisites in Ubuntu  (tested only in Ubuntu 22.04)

Run as root or as a used with sudo rights the following commands:
```
    # Packages needed for casacore development 
    $ sudo apt -y install git g++ cmake gfortran libreadline-dev libncurses5-dev libblas-dev liblapack-dev libfftw3-dev libcfitsio-dev wcslib-dev python3-numpy python-is-python3 flex bison tar curl

    # Packages needed for libsakura development
    $ sudo apt -y install libfftw3-dev libeigen3-dev

    # Additional packages needed for CASA development 
    $ sudo apt -y install default-jdk python3-scipy python3-matplotlib python3-certifi python3-pytest-xvfb python3-pytest python3-wheel python3-venv python3-build python3-numpy swig ccache pkg-config libxerces-c-dev libxml2-dev libxslt1-dev libgsl-dev libsqlite3-dev protobuf-compiler-grpc libgrpc-dev grpc-proto libgrpc++-dev wcslib-dev libopenmpi-dev openmpi-bin xvfb patchelf imagemagick libxml2-utils
```
#### Installing prerequisites in Debian  (tested only in Debian 11)

Run as root or as a used with sudo rights the following commands:
```
    # Packages needed for casacore development 
    $ sudo apt -y install git g++ cmake gfortran libreadline-dev libncurses5-dev libblas-dev liblapack-dev libfftw3-dev libcfitsio-dev wcslib-dev python3-numpy python-is-python3 flex bison tar curl

    # Packages needed for libsakura development
    $ sudo apt -y install libfftw3-dev libeigen3-dev

    # Additional packages needed for CASA development 
    $ sudo apt -y install default-jdk python3-scipy python3-matplotlib python3-certifi python3-pytest-xvfb python3-pytest python3-wheel python3-venv python3-build python3-numpy swig ccache pkg-config libxerces-c-dev libxml2-dev libxslt1-dev libgsl-dev libsqlite3-dev protobuf-compiler-grpc libgrpc-dev grpc-proto libgrpc++-dev wcslib-dev libopenmpi-dev openmpi-bin xvfb patchelf imagemagick libxml2-utils
```
#### Installing prerequisites in manylinux2014/RHEL7

Run as root or as a used with sudo rights the following commands:
```
    # Manylinux2014 Docker image
    # quay.io/pypa/manylinux2014_x86_64:2022-04-24-d28e73e

    # Update yum and add epel
    $ yum update -y
    $ yum install -y epel-release git cmake cmake3 gcc-c++ gcc-gfortran gtest-devel ccache readline-devel ncurses-devel blas-devel lapack-devel cfitsio-devel fftw-devel wcslib-devel flex bison tar curl libxml2-devel libxslt-devel java-1.8.0-openjdk-devel xorg-x11-server-Xvfb
    $ yum install -y xerces-c-devel

    # Get the newer gsl version from the casa03 repo
    $ yum install -y https://casa.nrao.edu/download/repo/el7/x86_64/casa03-gsl-2.5-1.el7.x86_64.rpm https://casa.nrao.edu/download/repo/el7/x86_64/casa03-gsl-devel-2.5-1.el7.x86_64.rpm

    # Install multithreaded MPI from the casa03 repo
    $ yum install -y https://casa.nrao.edu/download/repo/el7/x86_64/casa03-openmpi-1.10.4-1.el7.x86_64.rpm

    # Install Utilities and packages needed for libsakura installation
    $ yum install -y eigen3-devel fftw-devel man which time file tree rsync bc vim wget perf cppcheck valgrind kcachegrind patch mlocate git git-lfs gdb xterm screen tmux sysstat zsh ncdu numactl ncompress fuse fuse-devel libffi-devel
yum groupinstall -y "Development Tools"

    # Create local repository so that yum can solve the  grpc/protobuf depencies.
    $ mkdir /grpc
    $ yum --assumeyes install createrepo
    $ createrepo /grpc
    $ chmod -R o-w+r /grpc
    $ echo -en '[local]\nname=localgrpc\nbaseurl=file:///grpc\nenabled=1\ngpgcheck=0' > /etc/yum.repos.d/localgrpc.repo
    $ dnf -y install https://cbs.centos.org/kojifiles/packages/grpc/1.20.1/1.el7/x86_64/grpc-1.20.1-1.el7.x86_64.rpm 
    $ dnf -y install https://cbs.centos.org/kojifiles/packages/grpc/1.20.1/1.el7/x86_64/grpc-devel-1.20.1-1.el7.x86_64.rpm 
    $ dnf -y install https://cbs.centos.org/kojifiles/packages/grpc/1.20.1/1.el7/x86_64/grpc-plugins-1.20.1-1.el7.x86_64.rpm
    $ dnf -y install https://cbs.centos.org/kojifiles/packages/protobuf/3.6.1/4.el7/x86_64/protobuf-3.6.1-4.el7.x86_64.rpm
    $ dnf -y install https://cbs.centos.org/kojifiles/packages/protobuf/3.6.1/4.el7/x86_64/protobuf-compiler-3.6.1-4.el7.x86_64.rpm
    $ dnf -y install https://cbs.centos.org/kojifiles/packages/protobuf/3.6.1/4.el7/x86_64/protobuf-devel-3.6.1-4.el7.x86_64.rpm

    # Python 3.8
    $ wget https://www.python.org/ftp/python/3.8.7/Python-3.8.7.tar.xz
    $ tar -xJf Python-3.8.7.tar.xz
    $ cd Python-3.8.7 && ./configure && make && make install
    $ rm -rf  Python-3.8.7
    $ rm  Python-3.8.7.tar.xz

    $ scl enable devtoolset-10 bash
```
Ensure that the PATH points to the relevant software:
```
    $ export PATH=/opt/casa/03/bin:/usr/lib64/ccache:/opt/rh/devtoolset-10/root/usr/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
```
#### Installing prerequisites in macOS using macports (tested only on 11 and 12 Intel)

As a prerequisite [XCode](https://developer.apple.com/xcode/) must be installed, as well as [macports](https://macports.org). It could be that the prerequisites could be installed from Homebrew, but no attempt has been made to try it out.

The XCode installation already fulfills the requirement of the compiler (point 2. in the list of requirements).

When using any of the below commands, make sure that macports is in the PATH. That's part of the macports installation instructions so it could be that the PATH is already correct:
```
    $ export PATH=/opt/local/bin:$PATH
```
Run as root or as a used with sudo rights the following commands:
```
    # Packages needed for casacore development 
    $ sudo port install git cmake gcc12 +gfortran cfitsio wcslib 

    # Packages needed for libsakura development
    $ sudo port install fftw-3 fftw-3-single eigen3

    # Additional packages needed for CASA development 
    $ sudo port install ccache  py38-build py38-pip py38-numpy swig-python xercesc3 pkgconfig protobuf3-cpp grpc gsl libxslt openmpi-clang libxml2 fftw-3 fftw-3-single

    # Select python default version
    $ sudo port select --set python python38
    $ sudo port select --set python3 python38
    $ sudo port select --set pip pip38
    $ sudo port select --set pip3 pip38
```

(note that we exclude boost related packages that are optional in the casacore build instructions: libboost-dev, libboost-python-dev. We do not need those as CASA does not use the python bindings from casacore. In addition, if one wanted to compile all the boost related functionality of casacore (there is additional code that uses boost for Arrays and Dysco tests), the following packages would be needed: libboost-filesystem-dev libboost-test-dev libboost-system-dev).


### Setting up ccache

ccache significantly speeds up the compilation time for successive builds after the first build.

Depending on your OS the installation might be different. The recommended way is to have a link named after the compiler (like c++) that points to ccache. That link *must* be in the PATH before the proper compiler. For most Linux distributions one of the following commands should be enough:
```
    $ export PATH=/usr/lib64/ccache:$PATH # Redhat-like
    $ export PATH=/usr/lib/ccache:$PATH   # Debian-like
```
On NRAO workstations: Make sure that your cache directory is not in your home directory. The home directory is using the filer and tends to be slow for this purpose.

    $ export CCACHE_DIR=<not_your_home_directory>

Also, to avoid issues with some ccache versions, it is recommended to set the maximum size of ccache to some large value (please ensure that the ccache directory has enough space):

    $ export CCACHE_MAXSIZE=50G

### Install casacore and libsakura

To ease the instructions, define the variables CASAINSTALL, CASASRC and CASATESTDIR described above. The use of these variables is optional, but they help to follow the installation procedure as described here. The instructions to set the variables depend on the shell being used:
```
    $ export CASAINSTALL=/installation/path/  # (bash, zsh, POSIX shell)
    $ export CASASRC=/source/to/casa6/repo    # (bash, zsh, POSIX shell)
    $ export CASATESTDIR=/testdir/path/       # (bash, zsh, POSIX shell)
    $ export CASABUILD=/temporary/build/path/ # (bash, zsh, POSIX shell)

    $ setenv CASAINSTALL /installation/path/  # (csh, tcsh)
    $ setenv CASASRC /source/to/casa6/repo    # (csh, tcsh)
    $ setenv CASATESTDIR /testdir/path/       # (csh, tcsh)
    $ setenv CASABUILD /temporary/build/path/ # (csh, tcsh)
```

> macOS 13 /Ventura users: Add this to all of the CMake commands:
```
    -DCMAKE_CXX_FLAGS="-Qunused-arguments -flat_namespace" \
```

To install libsakura:

1. Get the sources
```
    $ cd $CASABUILD
    $ curl -L https://github.com/tnakazato/sakura/archive/refs/tags/libsakura-5.1.3.tar.gz | gunzip | tar -xvf -
```
1. Compile and install with cmake (you might change the build directory or the make options).
```
    $ cd sakura-libsakura*/libsakura
    $ mkdir build
    $ cd build
    $ cmake  \
        -DCMAKE_INSTALL_PREFIX=$CASAINSTALL \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DBUILD_DOC:BOOL=OFF \
        -DPYTHON_BINDING:BOOL=OFF \
        -DSIMD_ARCH=GENERIC \
        -DENABLE_TEST:BOOL=OFF \
         ..
    $ make install -j `getconf _NPROCESSORS_ONLN`
```

1. This will install libraries under ` $CASAINSTALL/lib ` and header files under ` $CASAINSTALL/include `

To install casacore:

Install measures data:
```
    $ mkdir -p $CASAINSTALL/data
    $ curl ftp://ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar | tar -C $CASAINSTALL/data -xzf -
```
Compile and install with cmake (you might change the build directory or the make options):
```
    $ mkdir $CASABUILD/casacore
    $ cd $CASABUILD/casacore
```
> NOTE! Macs have an extra cmake directive compared to linux, so be certain to use the cmake command below that is valid for your OS
```
    # LINUX ONLY! use the following cmake command 
    $ cmake \
        -DCMAKE_INSTALL_PREFIX=$CASAINSTALL \
        -DDATA_DIR=$CASAINSTALL/data \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DUSE_OPENMP=ON \
        -DUSE_THREADS=ON \
        -DBUILD_FFTPACK_DEPRECATED=ON \
        -DBUILD_TESTING=ON \
        -DBUILD_PYTHON3=OFF \
        -DBUILD_DYSCO=ON \
        -DPORTABLE=ON \
        -DUSE_PCH=OFF \
        -DUseCcache=1 \
        $CASASRC/casatools/casacore

    # manylinux2014 Note the extra gsl flags 
    cmake \
        ${extra_flags} \
        -DCMAKE_INSTALL_PREFIX=$INSTALLPREFIX \
        -DDATA_DIR=$INSTALLPREFIX/data \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DUSE_OPENMP=ON \
        -DUSE_THREADS=ON \
        -DBUILD_FFTPACK_DEPRECATED=ON \
        -DBUILD_TESTING=ON \
        -DBUILD_PYTHON3=OFF \
        -DBUILD_DYSCO=ON \
        -DPORTABLE=ON \
        -DUSE_PCH=OFF \
        -DUseCcache=1 \
        -DCMAKE_CXX_FLAGS="-I /opt/casa/03/include -L /opt/casa/03/lib/ -lgsl -lgslcblas" \
        -DCMAKE_EXE_LINKER_FLAGS="-L /opt/casa/03/lib/" \
        -DCMAKE_MODULE_LINKER_FLAGS="-Wl,-rpath,-L /opt/casa/03/lib/" \
        -DCMAKE_SHARED_LINKER_FLAGS="-L /opt/casa/03/lib/" \
        $CASASRC/casatools/casacore

    # MAC ONLY! use the following cmake command (see additional note below)
    # Note: we use '/opt/local/' as location where gcc is installed (MacPorts install prefix)
    $ export FC=/opt/local/bin/gfortran-mp-11  # (gcc version can be 11, 12, etc.)
    $ cmake \
        -DCMAKE_INSTALL_PREFIX=$CASAINSTALL \
        -DDATA_DIR=$CASAINSTALL/data \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DUSE_OPENMP=ON \
        -DUSE_THREADS=ON \
        -DBUILD_FFTPACK_DEPRECATED=ON \ 
        -DBUILD_TESTING=ON \
        -DBUILD_PYTHON3=OFF \
        -DBUILD_DYSCO=ON \
        -DPORTABLE=ON \
        -DUSE_PCH=OFF \
        -DUseCcache=1 \
        -DPRIVATE_LIBS="-framework Accelerate -lm -ldl -Wl,-rpath,/opt/local/lib/libgcc/" \
        $CASASRC/casatools/casacore

    # All OS
    $ make install -j `getconf _NPROCESSORS_ONLN`
```

1. This will install libraries under ` $CASAINSTALL/lib ` which are named like `libcasa_*` and header files under ` $CASAINSTALL/include `

>    NOTE: here we use a build of type "RelWithDebInfo" (with code optimization and debug line info. Other alternatives include Debug, Release, etc.

>    NOTE: The modular system allows to install casacore from other repository which is not the git submodule configured in casatools. This adds flexibility to get a customized casacore if needed. However, it is up to the developer to ensure that the version installed is the one that needs to be used for the subsequent casa C++ build. For instance, one could use a single casacore installation for several casa branches if that's practical or desired.

>    NOTE (macOS): In the macOS cmake command line we have to add '-DPRIVATE_LIBS="-framework Accelerate -lm -ldl"' as a temporary workaround until a fix can be added in casacore and used from casa. The problem is explained in (first comment, temporary workaround for an issue with LAPACK).
>    The additional linker rpath to lib/libgcc is needed for gcc 11 and 12, but not for older gcc versions (9). For gcc 9, The "-DPRIVATE_LIBS" can be set to simply -DPRIVATE_LIBS="-framework Accelerate -lm -ldl". The -rpath flag is needed in more modern versions of gcc to load libgfortran, or otherwise one would have to add that path to libgfortran in the DYLD_FALLBACK_LIBRARY_PATH environment variable.


### Install CASA C++ (aka casacpp)

1. Compile and install with cmake
```
    $ mkdir $CASABUILD/casacpp
    $ cd $CASABUILD/casacpp
    $ export PATH=/usr/lib64/openmpi/bin/:$PATH   # This is not needed in Debian or Ubuntu and is optional in the other platforms if no MPI support at the C++ level is needed

    # LINUX ONLY! (not manylinux2014, not macOS) use this cmake command
    $ PKG_CONFIG_PATH=$CASAINSTALL/lib/pkgconfig cmake \
         -DCMAKE_INSTALL_PREFIX=$CASAINSTALL \
         -DPKG_CONFIG_USE_CMAKE_PREFIX_PATH=$CASAINSTALL \
          $CASASRC/casatools/src/code

    $ manylinux2014 only!
    $ export PKG_CONFIG_PATH=$CASAINSTALL/lib/pkgconfig:/opt/casa/03/lib/pkgconfig
    $ cmake -DCMAKE_INSTALL_PREFIX=$CASAINSTALL \
          -DPKG_CONFIG_USE_CMAKE_PREFIX_PATH=$CASAINSTALL:/opt/casa/03/lib/pkgconfig \
          -DCMAKE_CXX_FLAGS="-I /opt/casa/03/include -L /opt/casa/03/lib/ -lgsl -lgslcblas -Wl,-rpath,/opt/casa/03/lib" \
          $CASASRC/casatools/src/code

    # MAC ONLY! use this cmake command
    # Note: we use '/opt/local/include' as location where to find the WCSLIB includes.
    $ PKG_CONFIG_PATH=$CASAINSTALL/lib/pkgconfig cmake \
         -DCMAKE_INSTALL_PREFIX=$CASAINSTALL \
         -DPKG_CONFIG_USE_CMAKE_PREFIX_PATH=$CASAINSTALL \
         -DCMAKE_CXX_FLAGS="-isystem /opt/local/include" \
          $CASASRC/casatools/src/code

    $ make install -j `getconf _NPROCESSORS_ONLN`
```

1. This will install libraries under ` $CASAINSTALL/lib ` which are named like `libcasacpp_*` and header files under ` $CASAINSTALL/include `

>    NOTE: In the macOS cmake command line we have to add '-DCMAKE_CXX_FLAGS="-isystem /opt/local/include"' as a temporary workaround to give the include path of wcslib, until a fix can be used from casacore. This will no longer be needed once the pointer to casacore is updated to include the fix in casa (). Alternatively, before running the casacore cmake one would need to  modify the following line inside $CASASRC/casatools/casacore/casacore.pc.in.

#### Run casacpp unit tests (optional)

1.    Go to the build directory.
```
    $ cd $CASABUILD/casacpp
```
1.    Export the variable CASADATA to point to the contents of the casatestdata repository
```
    $ export CASADATA=/path/to/casatestdata
```
1.    Run the tests
```
    $ ctest -T test --output-on-failure
```
1.    It is also possible to list all available tests. This command can be run from $CASABUILD/casacpp and will list all available tests. If it is run from within a subdirectory of it (for instance $CASABUILD/casacpp/singledish), then it will list the unit tests associated to that module:
```
    ctest -N
```
1.    To run a individual test:
```
    ctest -T test -R nameOfTest
```
1.    To run valgrind on a given test:
```
    ctest -T memcheck -R nameOfTest
```
1.    To run the debugger on a given test you can use the following one-liner :
```
    gdb `ctest -V -N -R nameOfTest | grep "Test command" | sed -r 's/.*Test command:(.*)/\1/g' `
```
### Create CASA casatools wheel

Please note that this procedure might be affected by the PYTHONPATH variable. Consider unsetting it.
1.    Create build directory
```
    $ mkdir $CASABUILD/casatools
    $ cd $CASABUILD/casatools 
```
1.    REQUIRED in Rocky Linux 8, Ubuntu 22.04 and MacOS ! OPTIONAL for other platforms. Create a virtual environment that has the needed packages:
```
    $ mkdir build_env
    $ python3 -m venv build_env
    $ . ./build_env/bin/activate
    $ pip install build setuptools wheel
```
1.    Remove the output directory to avoid confusion with old created wheels:
```
    $ rm -rf $CASAINSTALL/dist
```
1.    Create casatools wheel. Note that this assumes that the virtual environment in previous step (if required) is still active.
```
    $ PKG_CONFIG_PATH=$CASAINSTALL/lib/pkgconfig python3 -m build -o $CASAINSTALL/dist $CASASRC/casatools # Rocky Linux 8, Ubuntu 22.04 and MacOS with a venv

    $ PKG_CONFIG_PATH=$CASAINSTALL/lib/pkgconfig python3 -m build -o $CASAINSTALL/dist $CASASRC/casatools # Ubuntu and Fedora without a venv

    # WARNING: Add -C="--build-option=--mod-closure" to make wheels portable. This is required to build ManyLinux compatible wheels. However it will make the build slower and is not neccessary if you don't plan to distribute your wheel to someone else or install the wheel in a different computer than yours. The command line would then be:
    # PKG_CONFIG_PATH=$CASAINSTALL/lib/pkgconfig python3 -m build -n -o $CASAINSTALL/dist -C="--build-option=--mod-closure" $CASASRC/casatools


    # manylinux2014 only!
    PKG_CONFIG_PATH=/data/install/lib/pkgconfig:/opt/casa/03/lib/pkgconfig python3.8 -m build -n -o /data/install/dist ..

    # macOS
    PKG_CONFIG_PATH=$CASAINSTALL/lib/pkgconfig python3 -m build -n -o $CASAINSTALL/dist -C="--build-option=--mod-closure" $CASASRC/casatools
```
1. This will create an output wheel inside the `$CASAINSTALL/dist` directory. That wheel depends on the libraries installed under `$CASAINSTALL/lib`.

1. Optional: Convert Casatools wheel to ManyLinux compatible format
```
    # Add to your build_env
    pip install --upgrade pip auditwheel patchelf

    # ManyLinux 2.28
    auditwheel repair --exclude libglibmm-2.4.so.1 --exclude libglibmm-2.4.so.1.3.0 --exclude libblas.so.3 --exclude libffi.so.6 --exclude libfftw3.so.3 --exclude libfftw3f.so.3 --exclude libgmp.so.10 --exclude libgnutls.so.30 --exclude libgpr.so.7 --exclude libgrpc.so.7 --exclude libhogweed.so.4 --exclude libhwloc.so.15 --exclude libidn2.so.0 --exclude liblzma.so.5 --exclude libnettle.so.6 --exclude libpcre.so.1 --exclude libprofiler.so.0 --exclude libprotobuf.so.17 --exclude libssl.so.1.1 --exclude libtasn1.so.6 --exclude libtinfo.so.6 --exclude libunistring.so.2 --exclude libunwind.so.8 --exclude libxml2.so.2 --exclude libquadmath.so.0 --plat=manylinux_2_28_x86_64 $CASAINSTALL/dist/casatools*.whl

    # ManyLinux 2014
    auditwheel repair --exclude libglibmm-2.4.so.1 --exclude libglibmm-2.4.so.1.3.0 --exclude libblas.so.3 --exclude libffi.so.6 --exclude libfftw3.so.3 --exclude libfftw3f.so.3 --exclude libfftw3_threads.so.3 --exclude libfftw3f_threads.so.3 --exclude libfftw3f.so.3  --exclude libgmodule-2.0.so.0 --exclude libgpr.so.7 --exclude libgrpc.so.7 --exclude libgrpc++.so.1 --exclude libgssapi_krb5.so.2 --exclude libhwloc.so.5 --exclude libkeyutils.so.1 --exclude libltdl.so.7 --exclude liblzma.so.5 --exclude libnuma.so.1 --exclude libpcre.so.1 --exclude libprofiler.so.0 --exclude libprotobuf.so.17 --exclude libquadmath.so.0 --exclude libsigc-2.0.so.0 --exclude libssl.so.10 --exclude libtinfo.so.5 --exclude libxml2.so.2 $CASAINSTALL/dist/casatools*.whl
```
#### Test casatools (Optional)

1. (Optional) Create a virtual environment for testing purposes. It is recommended, although not strictly necessary. If done this way, only the python sessions that activate the environment will have access to casatools. Otherwise casatools will be installed in `$HOME` and be available to all python sessions which is probably not what most of developers want. You can use the --system-site-packages option to venv, which will use the python packages from your environment as installed with the instructions mentioned above using your package manager (or by other method you have used). However, for some platforms, including RHEL 8 and macOS that won't work out of the box and therefore is not recommended.
```
    $ python3 -m venv $CASATESTDIR/test_env # You can use python3 -m venv --system-site-packages $CASATESTDIR/test_env in Ubuntu and Fedora  
    $ . $CASATESTDIR/test_env/bin/activate 
```
1. Install the casatools wheel. NOTE: The uninstall command ensures that a potential previous installation is uninstalled first. Otherwise pip install won't install the new one if it thinks it has the same version. The pip uninstall command is harmless if this is the first time casatools is installed.
```
    $ pip uninstall casatools
    $ pip install $CASAINSTALL/dist/casatools*whl
    $ pip install casatestutils
    $ pip install casadata
```
1. Run the tests. Note that until CAS-13968 is not fixed, there will be failures when running all tools tests with pytest in test_tool_table. You can avoid these failures by running all tests with "python casatools/tests/run.py"
```
    $ python -m pytest $CASASRC/casatools/tests/tools/
```
### Create CASA casatasks wheel

The casatasks wheel creation and installation has not changed int he modular build system and is not yet fully PEP-517 compliant. Hence the build creation needs to take place on-source.

1.  Go to the casatasks source:
```
    $ cd $CASASRC/casatasks
```
1.    Remove the output directory to avoid confusion with old created wheels:
```
    $ rm -rf dist
```
1. Create casatasks wheel. The resulting wheel file is stored under $CASASRC/casatasks/dist. Requires the Python wheel module (question)
```
    $ ./setup.py bdist_wheel
```

1. This creates a casatasks wheel under `$CASASRC/casatasks/dist`. The current status of the casatasks setup.py is not fully PEP-517 compatible, that's why it is not possible to get the wheel under `$CASAINSTALL/dist` like the casatools case.

#### Test casatasks (Optional)

It is assumed that the steps to test casatools (see above) have already been performed.

1.     (Optional). Activate the same virtual environment that was used to test casatools (see above). Omit this step if you prefer to install in your $HOME.
```
    $ . $CASATESTDIR/test_env/bin/activate 
```
1.    Install casatasks wheel. NOTE: The uninstall command ensures that a potential previous installation is uninstalled first. Otherwise pip install won't install the new one if it thinks it has the same version. The pip uninstall command is harmless if this is the first time casatasks is installed.
```
    $ cd $CASASRC/casatasks
    $ pip uninstall casatasks
    $ pip install ./dist/casatasks*whl
```
1.    Optionally install casadata or point to a given location of casadata in your $HOME/.casa/config.py
```
    $ pip install casadata
```
1.     Run the tests:
```
    $ python $CASASRC/casatasks/tests/run.py
```

### Incremental builds for developers

The modular builds system uses standard cmake/make for building the casacore and casacpp modules. That means that if a change in the source code is done in any of those libraries, to recompile it is enough to do:
```
    $ make install -j `getconf _NPROCESSORS_ONLN`
```     
However, if the changes involve adding new files or directories then the cmake step needs to be re-run, using the cmake commands explained in the casacore or casaccpp section.

If the change is in casacore, then it is recommended to also run `make` in casacpp and create a new casatools wheel. If the change is in casacpp it is enough to create a new casatools wheel. Strictly speaking, creating a new casatools is not necessary if the change involves _only_ modifications of the C++ .cc files. However it is required if interface with casatools is changed, either due to changes in .h header files or template .tcc files.

If the casatools is created in a _portable_ way using the option `-C="--build-option=--mod-closure"` then the wheel itself contains a copy of the casacpp libraries and therefore the wheel needs to be recreated _always_ that there is a change in casacore or casacpp, even if the change involves only .cc files.

A new casatools wheel needs to be obviously recreated if there is a change in the tools directories, namely those under `casatools/xml, casatools/src/tools`.

Similarly, a new casatasks wheel would need to be created if there is any change under casatasks. 

If a new casatools wheel is created and needs to be tested using casatasks then the newly created casatools wheel needs to be reinstalled in the testing environment. The same is true if a new casatasks wheel is created. Since the newly generated wheel has the same version, you need to uninstall before:
```
    $ pip uninstall -y casatools
    $ pip install $CASAINSTALL/dist/casatools*whl
```

### Installation using Makefile

To streamline the installation of the different CASA modules a Makefile has been created that does everything in one go. If you use this Makefile you don't need to follow the different steps mentioned above, since they are basically encoded in the Makefile. However you still need to install the software prerequisites steps before attempting this, following the instructions outlines above.

The Makefile allows to have a first installation of a branch of CASA with the following steps:

1. Choose a root directory where the Makefile will be downloaded. If you don't modify the Makefile (see below), this directory will also be the root directory for cloning the git code, building and installing everything.
```
    $ mkdir my_casa_build
    $ cd my_casa_build
```
1. Download the Makefile script. This script is actually of the CASA6 repository and lies at the top level. You can get the latest version with this command:
```
    $ wget "https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/raw/Makefile?at=refs%2Fheads%2Fmaster" -O Makefile
```
1. Edit the Makefile to suit your needs. In particular, the most important variable is the first one: `CASA_BRANCH`. The rest have reasonable defaults that can be however modified to suite your needs. Note that the Makefile variables CASASRC, CASAINSTALL, CASABUILD, CASATESTDIR play the same role as the environmental variables mentioned in the manual steps above. They do not need to be redefined as, just modify in the Makefile if so needed. Please refer to sections above for a description of those variables.

1. Clone the git repository and build all CASA components casacore, casacpp, casatools, casatasks and casashell:

```
   $ make firstcasa
```
1. Use the venv under CASAVENVDIR dir to test the installation _or_ install the wheels under `CASAINSTALL/dist` (Substitute CASAINSTALL with the value in your Makefile):
```
   $ pip uninstall -y casatools casatasks
   $ pip install CASAINSTALL/dist
```

 It is also possible to run individual steps with the `make` targets `libsakura, casacore, casacpp, casatools, casatasks, casashell`. Those steps will depend on each other, so running `make casashell` will actually run all the other ones, which might be suboptimal. There are also individual make targets that do not depend on each other: `libsakura, casacore-configure, casacore-build, casacpp-configure, casacpp-build, casatools-wheel, casatasks-wheel, casashell-wheel`, Note that all these targets are there for convenience and they basically run the different modular steps, so sometimes it might be more convenient to use the standard commands used by the build system like `cmake`, `make` and `python -m build` tools. Note that once you have build everything with the `Makefile` you can use the individual steps mentioned in previous sections to rebuild casacpp, casatools, etc, as long as you set the environmental variables to the same values in the `Makefile`.

 The target `make clean` will clean all the working and output directories and can be used to start clean from scratch.

## Testing with Viewer and PlotMS in a development environment

To test interactive clean, or Viewer and PlotMS in general, a specific variants of the wheels need to be used. This is needed so that the GRPC versions are compatible between Casatools and Viewer/PlotMS. The previous version of the build system built GRPC as part of the Casatools package but the current system uses prebuilt rpms/ports.

It is also important to note that protobuf version 3.20.x is the latest supported version. This can be installed with pip: "pip install protobuf==3.20.3".

Linux versions of PlotMS and Viewer

https://casa-pip.nrao.edu/repository/casa-test-wheel/packages/casaviewer/1.7.1/casaviewer-1.7.1-py3-none-manylinux_2_5_x86_64.manylinux1_x86_64.manylinux_2_17_x86_64.manylinux2014_x86_64.whl

https://casa-pip.nrao.edu/repository/casa-test-wheel/packages/casaplotms/1.9.1/casaplotms-1.9.1-py3-none-manylinux_2_5_x86_64.manylinux1_x86_64.manylinux_2_17_x86_64.manylinux2014_x86_64.whl


macOS version of PlotMS (Viewer is no longer built on Macs).

https://casa-pip.nrao.edu/repository/casa-test-wheel/packages/casaplotms/1.9.1/casaplotms-1.9.1-py3-none-macosx_12_0_x86_64.whl

When a feature requires both casacore and casa change
 
Making changes in Casacore 


1. Create a Casacore fork in GitHub

2. Create a Casa branch in BitBucket

3. Clone the repository and checkout your branch


git clone --recursive https://open-bitbucket.nrao.edu/scm/casa/casa6.git
cd casa6/casatools
git checkout CAS-1234


4. Create a casacore branch

cd casacore
git checkout master
git pull
git branch mycasacorefeature
git checkout mycasacorefeature


5. Make your changes in casacore

6. Make your changes in the rest of the branch

7. Test locally

8. Push the Casacore changes to your fork in GitHub

cd casa6/casacore
git remote add mycasacore https://github.com/vsuorant/casacore
git push mycasacore mycasacorefeature


9. Create a pull request in GitHub

10. Wait for the pull request to be applied (you must wait since the master submodule doesn't know about your fork, so you can't point the submodule there)

11. Update the submodule reference in your branch

cd casa6/casatools
git checkout CAS-1234
cd casacore
git checkout master
git pull
cd ..
git add casacore
git commit --amend (this will amend your latest commit. If you would rather have a separate commit, leave the --amend out)


12. Push your changes to BitBucket

git push origin CAS-1234


Switching Casacore remotes and branches

Sometimes you need or want to add more remotes for Casacore changes. To add a remote do:

git remote add mycasacore https://github.com/vsuorant/casacore
If you want to make your casacore master "track" the master in the new remote, do the following:

git fetch mycasacore
 
git checkout -B master central-casacore/master



## CASA's GUIs

The primary GUIs [casaviewer](https://open-bitbucket.nrao.edu/projects/CASA/repos/casaviewer/browse)
and [casaplotms](https://open-bitbucket.nrao.edu/projects/CASA/repos/casaplotms/browse) are available
as separate modules. These allow the GUIs to be loaded and used.

Some of the GUI tools available as part of CASA 5 are only available in the packaged tar-file based
distribution of CASA 6. The GUIs that are currently only available in the monolithic version of
CASA 6 currently includes:

  1. casabrowser
  2. casafeather
  3. casalogger
  4. casaplotserver
