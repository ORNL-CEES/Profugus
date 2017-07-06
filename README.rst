Description
***********

Profugus is an open-source mini-application (mini-app) for radiation
transport and reactor applications.  It contains the fundamental
computational kernels used in the Exnihilo code suite from Oak Ridge
National Laboratory. However, Exnihilo is production code with a
substantial user base. Furthermore, Exnihilo is export controlled.
This makes collaboration with computer scientists and computer
engineers difficult.  Profugus is designed to bridge that gap.  By
encapsulating the core numerical algorithms in an abbreviated code
base that is open-source, computer scientists can analyze the
algorithms and easily make code-architectural changes to test
performance without compromising the production code values of
Exnihilo.

Profugus is **not** meant to be production software with respect to
problem analysis.  The computational kernels in Profugus are designed
to analyze *performance*, not *correctness*.  Nonetheless, users of
Profugus can setup and run problems with enough *real-world* features
to be useful as proof-of-concept for actual production work.

Profugus is also used as a test-bed mini-app for the ASCR projects
XPRESS and MCREX.


Profugus Development Team
=========================

Profugus contains the computational kernels from the Exnihilo code
base. The core Exnihilo development team consists of the following
scientists (listed alphabetically)

* Kevin Clarno <clarnokt@ornl.gov>

* Greg Davidson <davidsongg@ornl.gov>

* Tom Evans <evanstm@ornl.gov>

* Steven Hamilton <hamiltonsp@ornl.gov>

* Seth Johnson <johnsonsr@ornl.gov>

* Tara Pandya <pandyatm@ornl.gov>

* Stuart Slattery <slatterysr@ornl.gov>

* Rachel Slaybaugh <slaybaugh@berkeley.edu>

Profugus is developed and maintained by the following team:

* Tom Evans <evanstm@ornl.gov>

* Steven Hamilton <hamiltonsp@ornl.gov>

* Stuart Slattery <slatterysr@ornl.gov>


Profugus Packages
=================

Profugus constains the following code mini-apps:

**SPn**
   Mini-app of Simplified Spherical Harmonics ("SPn") computational
   kernel.

**MC**
   Mini-app of Monte Carlo computational kernel.

Documentation on the mini-apps (methods and execution) are given in
*Running The Mini-Applications*. There are several support packages
that are used by the mini-apps.  These are:

**Utils**
   Utilities and testing framework used by other packages.

**Matprop**
   Data structures for storing material properties (cross-sections).

Profugus is designed to build and run with a minimum of dependencies.
However, there are some requirements.  The third-party software (TPLs)
necessary to build Profugus is all open-source and freely available.
The TPLs for Profugus are listed in the following table:

+-----------------------+---------------+---------------------------------------+
| TPL                   | Required      | Comments                              |
+=======================+===============+=======================================+
| TriBITS               | Yes           | TriBITS is the Trilinos build system  |
+-----------------------+---------------+---------------------------------------+
| Trilinos              | Yes           | Trilinos is an open-source solvers    |
+-----------------------+---------------+---------------------------------------+
| BLAS/LAPACK           | Yes           | Use vendor-specific implementation    |
+-----------------------+---------------+---------------------------------------+
| MPI                   | No            | OpenMPI and MPICH are suggested       |
+-----------------------+---------------+---------------------------------------+
| HDF5                  | No            | Parallel HDF5 will allow output of    |
+-----------------------+---------------+---------------------------------------+

Building The Code
*****************

The most straightforward method for building Profugus is to use the
scripts in `Profugus/install`.  Profugus uses the `TriBITS` build
system.  This system is a set of package-based extensions to standard
cmake.  So, first you need to obtain *Trilinos* and *TriBITS* and put
them in your top-level Profugus directory::

   > cd Profugus
   > git clone https://github.com/TriBITSPub/TriBITS.git
   > ln -s $PATH_TO_TRILINOS .

The preferred mechanism for using the build scripts is to make a
*target* directory where the build is to be performed::

   > pwd
     /home/me
   > mkdir debug
   > cd debug
   > mkdir target
   > cd target
   > pwd
     /home/me/debug/target

The `install` directory contains several example build scripts.
General options for all platforms (which can be overridden at
configure time) are specified in the `install/base.cmake`::

   ##---------------------------------------------------------------------------##
   ## CMAKE BASE FILE
   ##---------------------------------------------------------------------------##

   # Default build all packages
   SET(Profugus_ENABLE_Utils   ON CACHE BOOL "")
   SET(Profugus_ENABLE_Matprop ON CACHE BOOL "")
   SET(Profugus_ENABLE_SPn     ON CACHE BOOL "")
   SET(Profugus_ENABLE_MC      ON CACHE BOOL "")

   # Turn on tests
   SET(Profugus_ENABLE_TESTS ON CACHE BOOL "")
   SET(Profugus_TEST_CATEGORIES "BASIC" CACHE STRING "")

   # Turn on SS code and optional packages by default
   SET(Profugus_ENABLE_ALL_FORWARD_DEP_PACKAGES OFF CACHE BOOL "")
   SET(Profugus_ENABLE_ALL_OPTIONAL_PACKAGES    ON  CACHE BOOL "")
   SET(Profugus_ENABLE_SECONDARY_STABLE_CODE    ON  CACHE BOOL "")

   # Up the max num procs
   SET(MPI_EXEC_MAX_NUMPROCS 8 CACHE STRING "")

   # Turn off binutils
   SET(Teuchos_ENABLE_BinUtils OFF CACHE BOOL "")

   # Turn off Zoltan2
   SET(Profugus_ENABLE_Zoltan2 OFF CACHE BOOL "")

   # Compiler options
   SET(BUILD_SHARED_LIBS ON CACHE BOOL "")
   SET(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated-declarations" CACHE STRING "")
   SET(Profugus_ENABLE_CXX11:BOOL=ON)

   # TriBITS stuff
   SET(Profugus_ENABLE_INSTALL_CMAKE_CONFIG_FILES OFF CACHE BOOL "")
   SET(Profugus_DEPS_XML_OUTPUT_FILE "" CACHE FILEPATH "")

By default, all of the packages inside of Profugus are turned on.
Furthermore, *C++-11* is **required**.  The default options specify
the appropriate compiler flags for gcc.  The tests are also turned on
by default; to disable tests in any upstream package simply do not
explicitly *ENABLE* that package.  For example, to build the *SPn*
package and all of its tests but only include required *source* from
upstream packages, the user would specify::

   SET(Profugus_ENABLE_SPn ON CACHE BOOL "")

In this case, only the pieces of *Utils* needed to build *SPn* are
compiled. All tests can be turned off by setting
**Profugus_ENABLE_TESTS** to **OFF**.

The `install` directory contains several build scripts that are all
suffixed by the platform name.  For example, to build on a Linux
*x86_64* system the "install/cmake_x86_64.sh" script can be used::

   #!/bin/sh
   ##---------------------------------------------------------------------------##
   ## CMAKE FOR X86_64
   ##---------------------------------------------------------------------------##

   # CLEANUP
   rm -rf CMakeCache.txt
   rm -rf CMakeFiles

   # SOURCE AND INSTALL
   SOURCE=<SET_SOURCE_DIR>
   INSTALL=<SET_INSTALL_DIR>

   # BUILD OPTIONS
   BUILD="DEBUG"
   MPI="ON"

   # TPL PATHS
   HDF5_PATH="/vendors/hdf5_parallel"
   MPI_PATH="/opt/openmpi/gcc/current"

   ##---------------------------------------------------------------------------##

   cmake \
   -DCMAKE_BUILD_TYPE:STRING="$BUILD" \
   -DTPL_ENABLE_MPI:BOOL=$MPI \
   -DCMAKE_INSTALL_PREFIX:PATH=$INSTALL \
   \
   -DMPI_BASE_DIR:PATH=$MPI_PATH \
   \
   -DTPL_ENABLE_HDF5:BOOL=ON \
   -DHDF5_INCLUDE_DIRS:PATH=$HDF5_PATH/include \
   -DHDF5_LIBRARY_DIRS:PATH=$HDF5_PATH/lib \
   \
   -DBLAS_LIBRARY_DIRS:PATH=/vendors/gcc/atlas/lib \
   -DLAPACK_LIBRARY_DIRS:PATH=/vendors/gcc/atlas/lib \
   -DBLAS_LIBRARY_NAMES:STRING="f77blas;cblas;atlas" \
   -DLAPACK_LIBRARY_NAMES:STRING="lapack" \
   \
   -DProfugus_CONFIGURE_OPTIONS_FILE:FILEPATH="${SOURCE}/install/base.cmake" \
   -DProfugus_ASSERT_MISSING_PACKAGES:BOOL=OFF \
   \
   ${SOURCE}

   ##---------------------------------------------------------------------------##
   ## end of cmake_x86_64.sh
   ##---------------------------------------------------------------------------##

The source and install locations must be set. Also, to enable a
optimized build set **BUILD** to **RELEASE**.  Adjust the paths and
libraries for LAPACK to fit your platform.  The example assumes that
the ATLAS LAPACK is available.  Any standard LAPACK distribution will
work. HDF5 is **not** required, to build/run/test the applications;
however, problem output will be severely curtailed if a parallel HDF5
option is not provided.  If HDF5 is not available, setting::

   -DTPL_ENABLE_HDF5:BOOL=OFF \

will disable HDF5.

To complete the configuration, execute this script inside the *target*
directory and then make/test/install::

   > pwd
     /home/me/debug/target
   > sh /home/me/Profugus/install/cmake_x86_64.sh
   > make -j 8
   > ctest -j 8
   > make -j 8 install
