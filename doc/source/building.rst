*****************
Building The Code
*****************

.. highlight:: sh

The most straightforward method for building |Profugus| is to use the scripts
in :file:`Profugus/install`.  |Profugus| uses the :download:`TriBits build
system <../../Trilinos/TrilinosBuildQuickRef.pdf>`.  This system is a set of
package-based extensions to standard cmake_.

The preferred mechanism for using the build scripts is to make a *target*
directory where the build is to be performed::

 > pwd
   /home/me
 > mkdir debug
 > cd debug
 > mkdir target
 > cd target
 > pwd
   /home/me/debug/target

The :file:`install` directory contains several example build scripts. General
options for all platforms (which can be overridden at configure time) are
specified in the :file:`install/base.cmake`:

.. literalinclude:: install/base.cmake
   :language: cmake

By default, all of the packages inside of |Profugus| are turned on.
Furthermore, `C++-11` is **required**.  The default options specify the
appropriate compiler flags for gcc_.  The tests are also turned on by default;
to disable tests in any upstream package simply do not explicitly `ENABLE`
that package.  For example, to build the `SPn` package and all of its tests
but only include required *source* from upstream packages, the user would
specify:

.. code-block:: cmake

   SET(Profugus_ENABLE_SPn ON CACHE BOOL "")

In this case, only the pieces of `Utils` needed to build `SPn` are compiled.
All tests can be turned off by setting :makevar:`Profugus_ENABLE_TESTS` to
:makevar:`OFF`.

The :file:`install` directory contains several build scripts that are all
suffixed by the platform name.  For example, to build on a Linux `x86_64`
system the :file:`install/cmake_x86_64.sh` script can be used:

.. literalinclude:: install/cmake_x86_64.sh
   :language: sh

The source and install locations must be set. Also, to enable a optimized
build set :makevar:`BUILD` to :makevar:`RELEASE`.  Adjust the paths and
libraries for LAPACK_ to fit your platform.  The example assumes that the
ATLAS_ LAPACK_ is available.  Any standard LAPACK_ distribution will
work. HDF5_ is **not** required, to build/run/test the applications; however,
problem output will be severely curtailed if a parallel HDF5_ option is not
provided.  If HDF5_ is not available, setting::

 -DTPL_ENABLE_HDF5:BOOL=OFF \

will disable HDF5_.

To complete the configuration, execute this script inside the *target*
directory and then make/test/install::

 > pwd
   /home/me/debug/target
 > sh /home/me/Profugus/install/cmake_x86_64.sh
 > make -j 8
 > ctest -j 8
 > make -j 8 install

.. _cmake: http://www.cmake.org
.. _LAPACK: http://www.netlib.org
.. _ATLAS: http://math-atlas.sourceforge.net
.. _gcc: http://gcc.gnu.org
.. _HDF5: http://www.hdfgroup.org/HDF5
