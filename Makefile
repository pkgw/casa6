##
# CASA6 Modular Makefile
# The intention of this file is to serve as a quick way to put together
# the different pieces of the modular build system to have a full build
# up to casashell with a single "make" command. The target users is
# developers with an understanding of the CASA build procedure.
#
# While the Makefile should work out of the box it is meant to be
# customized fo the individual developer needs. For instance, it
# will checkout the code from git, which might not be neccesary
# for some developers.
#
# To customize for a given branch build, please change the
# CASA_BRANCH variable below to the branch of interest.
# All the steps, from checking out the code to the build directories
# uses directories under ROOT, which by default points to the current
# working directory. This can be changed modriyfying the ROOT variable
# below.
#
# To get a full build from scratch, type "make firstcasa" in a
# directory which contains only this Makefile. Afterwards, the individual
# targets can be used to run different steps. Note that this Makefile
# to perform incremental builds of casacpp one would type
# "make casacpp-build" after the first "make firstcasa"
#
# This Makefile does not run any test.

CASA_BRANCH         = master
CASA_REPO           = https://open-bitbucket.nrao.edu:/scm/casa/casa6.git
CASACORE_DATA_REPO  = ftp://ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar

CASA_BUILD_TYPE     = RelWithDebInfo
CASACORE_BUILD_TYPE = RelWithDebInfo

LIBSAKURA_VERSION   = 5.2.1
CASASHELL_BRANCH    = master

# Number of cores used for compilation (default: all available in the machine)
NCORES              = $(shell getconf _NPROCESSORS_ONLN)

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#--------------------------------------------------------------------------------------------------------
#
# Package-level dir structure
ROOT        = $(shell pwd)

SRCDIR      = $(ROOT)/src
CASASRC     = $(SRCDIR)/casa6
CASAINSTALL = $(ROOT)/install
CASATESTDIR = $(ROOT)/test
CASAVENVDIR = $(ROOT)/venv
CASABUILD   = $(ROOT)/build

#INSTALLPREFIX  = $(CASAINSTALL)
#
# Common options to install artifacts of all packages in a single location
#
#INSTALLOPTS = -DCMAKE_INSTALL_PREFIX=$(INSTALLPREFIX) \
#		  -DCMAKE_INSTALL_BINDIR=sbin \
#		  -DCMAKE_INSTALL_LIBDIR=lib


#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

firstcasa: init casa-clone libsakura casacore casacpp venv-build casatools casatasks casashell
	echo CASA has been built successfully.
	You can run it with:
	\$ . $(CASAVENVDIR)/bin/activate
	\$python
	>>> import casatasks

casa: libsakura casacore casacpp venv-build casatools casatasks casashell

clean:
	rm -rf $(SRCDIR) $(CASASRC) $(CASABUILD) $(CASAINSTALL) $(CASATESTDIR) $(CASAVENVDIR)

init:
	mkdir -p $(SRCDIR) $(CASASRC) $(CASABUILD) $(CASAINSTALL) $(CASATESTDIR) $(CASAVENVDIR)

casa-clone:
	git -C $(SRCDIR) clone -b $(CASA_BRANCH) --recursive $(CASA_REPO)

libsakura:
	curl -L https://github.com/tnakazato/sakura/archive/refs/tags/libsakura-$(LIBSAKURA_VERSION).tar.gz | gunzip | tar -xvf - -C $(SRCDIR)

	mkdir -p $(CASABUILD)/libsakura
	cmake  \
		-DCMAKE_INSTALL_PREFIX=$(CASAINSTALL) \
		-DCMAKE_BUILD_TYPE=$(CASA_BUILD_TYPE) \
		-DBUILD_DOC:BOOL=OFF \
		-DPYTHON_BINDING:BOOL=OFF \
		-DSIMD_ARCH=GENERIC \
		-DENABLE_TEST:BOOL=OFF \
		-DUseCcache=1 \
		$(SRCDIR)/sakura-libsakura-$(LIBSAKURA_VERSION)/libsakura/ \
		-B $(CASABUILD)/libsakura

	$(MAKE) -C $(CASABUILD)/libsakura  install -j $(NCORES)


casacore: casacore-build casacore-configure

casacore-configure:
	if [ ! -d $(CASAINSTALL)/data ]; then \
		mkdir -p $(CASAINSTALL)/data ; \
		curl -L $(CASACORE_DATA_REPO) | gunzip | tar -xvf - -C $(CASAINSTALL)/data ; \
	fi

	mkdir -p $(CASABUILD)/casacore
	cd $(CASABUILD)/casacore
	cmake \
		-DCMAKE_INSTALL_PREFIX=$(CASAINSTALL) \
		-DDATA_DIR=$(CASAINSTALL)/data \
		-DCMAKE_BUILD_TYPE=$(CASA_BUILD_TYPE) \
		-DCMAKE_BUILD_PREFIX=$(CASABUILD) \
		-DUSE_OPENMP=ON \
		-DUSE_THREADS=ON \
		-DBUILD_FFTPACK_DEPRECATED=ON \
		-DBUILD_TESTING=ON \
		-DBUILD_PYTHON3=OFF \
		-DBUILD_DYSCO=ON \
		-DPORTABLE=ON \
		-DUSE_PCH=OFF \
		-DUseCcache=1 \
		$(CASASRC)/casatools/casacore \
		-B $(CASABUILD)/casacore

casacore-build : casacore-configure
	$(MAKE) -C $(CASABUILD)/casacore install -j $(NCORES)


casacpp: libsakura casacore casacpp-build

casacpp-needs-configure: $(CASABUILD)/casacpp/Makefile

casacpp-configure: clean_casacpp_build $(CASABUILD)/casacpp/Makefile

clean-casacpp-build:
	rm -rf $(CASABUILD)/casacpp

$(CASABUILD)/casacpp/Makefile:
	if [ -d $(CASABUILD)/casacpp ]; then rm -rf $(CASABUILD)/casacpp; fi
	mkdir -p $(CASABUILD)/casacpp
	cd $(CASABUILD)/casacpp
	PATH=/usr/lib64/openmpi/bin/:$(PATH) \
		PKG_CONFIG_PATH=$(CASAINSTALL)/lib/pkgconfig \
		cmake \
		-DCMAKE_INSTALL_PREFIX=$(CASAINSTALL) \
		-DPKG_CONFIG_USE_CMAKE_PREFIX_PATH=$(CASAINSTALL) \
		-DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
		$(CASASRC)/casatools/src/code \
		-B $(CASABUILD)/casacpp


casacpp-build : casacpp-needs-configure
	$(MAKE) -C $(CASABUILD)/casacpp  install -j $(NCORES)

venv-build: $(CASAVENVDIR)/bin/activate

$(CASAVENVDIR)/bin/activate:
	-deactivate # Disable any running virtual environments
	python3 -m venv $(CASAVENVDIR)
	. $(CASAVENVDIR)/bin/activate

casatools: casacpp casatools-wheel

casatools-wheel: venv-build
	if [ -d $(CASABUILD)/casatools ]; then rm -rf $(CASABUILD)/casatools; fi
	mkdir -p $(CASABUILD)/casatools
	if [ -d $(CASAINSTALL)/dist ]; then rm -rf $(CASAINSTALL)/dist; fi
	mkdir -p $(CASAINSTALL)/dist

	# Run the build in a venv
	# Disable any potentially running virtual environments beforehand
	-deactivate ; \
		. $(CASAVENVDIR)/bin/activate ; \
		pip install build ; \
		export CMAKE_BUILD_PARALLEL_LEVEL=$(NCORES) ; \
		cd $(CASABUILD)/casatools; PKG_CONFIG_PATH=$(CASAINSTALL)/lib/pkgconfig python3 -m build -o $(CASAINSTALL)/dist $(CASASRC)/casatools ; \
		pip uninstall -y casatools ; \
		pip install $(CASAINSTALL)/dist/casatools*whl ; \
		pip install casadata ; \
		deactivate

casatasks: casatools casatasks-wheel

casatasks-wheel: venv-build
	if [ -d $(CASASRC)/casatasks/dist ]; then rm -rf $(CASASRC)/casatasks/dist; fi
	mkdir -p $(CASASRC)/casatasks/dist
	if [ -d $(CASASRC)/casatasks/build ]; then rm -rf $(CASASRC)/casatasks/build; fi

	# Run the build in a venv
	# Disable any potentially running virtual environments beforehand
	-deactivate ; \
		. $(CASAVENVDIR)/bin/activate ; \
		pip install --upgrade setuptools ; \
		pip install --upgrade wheel ; \
		cd $(CASASRC)/casatasks ; \
		./setup.py bdist_wheel ; \
		pip uninstall -y casatasks ; \
		pip install $(CASASRC)/casatasks/dist/casatasks*whl ; \
		\cp -f $(CASASRC)/casatasks/dist/casatasks*whl $(CASAINSTALL)/dist

casashell: casatasks casashell-wheel

casashell-wheel: venv-build
	if [ -d $(SRCDIR)/casashell ]; then rm -rf $(SRCDIR)/casashell; fi
	git -C $(SRCDIR) clone -b $(CASASHELL_BRANCH) --recursive https://open-bitbucket.nrao.edu/scm/casa/casashell.git

	# Run the build in a venv
	# Disable any potentially running virtual environments beforehand
	-deactivate ; \
		. $(CASAVENVDIR)/bin/activate ; \
		cd $(SRCDIR)/casashell ; \
		./setup.py bdist_wheel ; \
		pip uninstall -y casashell ; \
		pip install $(SRCDIR)/casashell/dist/casashell*whl ; \
		\cp -f $(SRCDIR)/casashell/dist/casashell*whl $(CASAINSTALL)/dist

# end
