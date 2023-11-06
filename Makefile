##
# CASA6 Modular Makefile
#
# @file
# @version 0.1

CASA_BRANCH		= master
CASA_REPO	      	= https://open-bitbucket.nrao.edu:/scm/casa/casa6.git
CASACORE_DATA_REPO 	= ftp://ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar

CASA_BUILD_TYPE     = RelWithDebInfo
CASACORE_BUILD_TYPE = RelWithDebInfo

LIBSAKURA_VERSION	= 5.1.3
CASASHELL_BRANCH	= master

# Number of cores used for compilation (default: all available in the machine)
NCORES = $(shell getconf _NPROCESSORS_ONLN)

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#--------------------------------------------------------------------------------------------------------
#
# Package-level dir structure
ROOT           	= $(shell pwd)

SRCDIR		= ${ROOT}/src
CASASRC        	= ${SRCDIR}/casa6
CASAINSTALL	= ${ROOT}/install
CASATESTDIR	= ${ROOT}/test
CASAVENVDIR     = ${ROOT}/venv
CASABUILD	= ${ROOT}/build

#INSTALLPREFIX  = ${CASAINSTALL}
#
# Common options to install artifacts of all packages in a single location
#
#INSTALLOPTS = -DCMAKE_INSTALL_PREFIX=${INSTALLPREFIX} \
#		  -DCMAKE_INSTALL_BINDIR=sbin \
#		  -DCMAKE_INSTALL_LIBDIR=lib


#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
.ONESHELL: # Execute all make commands in the same shell. No need to end every line with \


firstcasa: init casa-clone libsakura casacore casacpp venv-build casatools casatasks casashell
casa: libsakura casacore casacpp venv-build casatools casatasks casashell

clean:
	rm -rf ${SRCDIR} ${CASASRC} ${CASABUILD} ${CASAINSTALL} ${CASATESTDIR} ${CASAVENVDIR}

init:
	mkdir -p  ${SRCDIR} ${CASASRC} ${CASABUILD} ${CASAINSTALL} ${CASATESTDIR} ${CASAVENVDIR}

casa-clone:
	cd ${SRCDIR}; git clone -b ${CASA_BRANCH} --recursive ${CASA_REPO}

libsakura:
	cd ${SRCDIR};
	curl -L https://github.com/tnakazato/sakura/archive/refs/tags/libsakura-${LIBSAKURA_VERSION}.tar.gz | gunzip | tar -xvf -

	mkdir -p ${CASABUILD}/libsakura
	cd ${CASABUILD}/libsakura
	cmake  \
	-DCMAKE_INSTALL_PREFIX=${CASAINSTALL} \
	-DCMAKE_BUILD_TYPE=${CASA_BUILD_TYPE} \
	-DBUILD_DOC:BOOL=OFF \
	-DPYTHON_BINDING:BOOL=OFF \
	-DSIMD_ARCH=GENERIC \
	-DENABLE_TEST:BOOL=OFF \
	-DUseCcache=1 \
	${SRCDIR}/sakura-libsakura-${LIBSAKURA_VERSION}/libsakura/

	make install -j ${NCORES}


casacore: casacore-build casacore-configure

casacore-configure:
	cd ${CASAINSTALL}

	if [ ! -d ./data ]; then
		mkdir -p ./data
		cd ./data
		curl -L ${CASACORE_DATA_PATH} | gunzip | tar -xvf -
		cd ../
	fi

	mkdir -p ${CASABUILD}/casacore
	cd ${CASABUILD}/casacore
	cmake \
	-DCMAKE_INSTALL_PREFIX=${CASAINSTALL} \
	-DDATA_DIR=${CASAINSTALL}/data \
	-DCMAKE_BUILD_TYPE=${CASA_BUILD_TYPE} \
	-DCMAKE_BUILD_PREFIX=${CASABUILD} \
	-DUSE_OPENMP=ON \
	-DUSE_THREADS=ON \
	-DBUILD_FFTPACK_DEPRECATED=ON \
	-DBUILD_TESTING=ON \
	-DBUILD_PYTHON3=OFF \
	-DBUILD_DYSCO=ON \
	-DPORTABLE=ON\
	-DUSE_PCH=OFF \
	-DUseCcache=1 \
	${CASASRC}/casatools/casacore

casacore-build : casacore-configure
	cd ${CASABUILD}/casacore
	make install -j ${NCORES}


casacpp: libsakura casacore casacpp-build

casacpp-needs-configure: $(CASABUILD)/casacpp/Makefile

casacpp-configure: clean_casacpp_build $(CASABUILD)/casacpp/Makefile

clean-casacpp-build:
	rm -rf $(CASABUILD)/casacpp

$(CASABUILD)/casacpp/Makefile:
	if [ -d ${CASABUILD}/casacpp ]; then rm -rf ${CASABUILD}/casacpp; fi
	mkdir -p ${CASABUILD}/casacpp
	cd ${CASABUILD}/casacpp
	PATH=/usr/lib64/openmpi/bin/:${PATH}
	PKG_CONFIG_PATH=${CASAINSTALL}/lib/pkgconfig cmake \
		-DCMAKE_INSTALL_PREFIX=${CASAINSTALL} \
		-DPKG_CONFIG_USE_CMAKE_PREFIX_PATH=${CASAINSTALL} \
		-DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
		${CASASRC}/casatools/src/code

casacpp-build : casacpp-needs-configure
	cd ${CASABUILD}/casacpp
	make install -j ${NCORES}

venv-build: ${CASAVENVDIR}/bin/activate

${CASAVENVDIR}/bin/activate:
	python3 -m venv ${CASAVENVDIR}
	source ${CASAVENVDIR}/bin/activate

casatools: casacpp casatools-wheel

casatools-wheel: venv-build
	-deactivate # Disable any running virtual environments
	source ${CASAVENVDIR}/bin/activate

	if [ -d ${CASABUILD}/casatools ]; then rm -rf ${CASABUILD}/casatools; fi
	mkdir -p ${CASABUILD}/casatools
	mkdir -p ${CASAINSTALL}/dist
	cd ${CASABUILD}/casatools

	pip install build
	export CMAKE_BUILD_PARALLEL_LEVEL=${NCORES}

	PKG_CONFIG_PATH=${CASAINSTALL}/lib/pkgconfig python3 -m build -o ${CASAINSTALL}/dist ${CASASRC}/casatools

	pip uninstall -y casatools
	pip install ${CASAINSTALL}/dist/casatools*whl
	pip install casadata
	deactivate

casatasks: casatools casatasks-wheel

casatasks-wheel: venv-build
	-deactivate
	source ${CASAVENVDIR}/bin/activate
	pip install --upgrade setuptools
	pip install --upgrade wheel
	cd ${CASASRC}/casatasks

	if [ -d ./dist ]; then rm -rf ./dist; fi
	if [ -d ./build ]; then rm -rf ./build; fi

	./setup.py bdist_wheel

	pip uninstall -y casatasks
	pip install ./dist/casatasks*whl
	\cp -f ./dist/casatasks*whl ${CASAINSTALL}/dist

casashell: casatasks casashell-wheel

casashell-wheel: venv-build
	-deactivate
	source ${CASAVENVDIR}/bin/activate

	cd ${SRCDIR}
	if [ -d casashell ]; then rm -rf casashell; fi
	git clone -b ${CASASHELL_BRANCH} --recursive https://open-bitbucket.nrao.edu/scm/casa/casashell.git
	cd casashell/
	./setup.py bdist_wheel
	pip install ./dist/casashell*whl


# end
