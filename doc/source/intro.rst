***********
Description
***********

|Profugus| is an open-source mini-application (mini-app) for radiation
transport and reactor applications.  It contains the fundamental computational
kernels used in the Exnihilo_ code suite from Oak Ridge National Laboratory.
However, Exnihilo_ is production code with a substantial user base.
Furthermore, Exnihilo_ is export controlled.  This makes collaboration with
computer scientists and computer engineers difficult.  |Profugus| is designed
to bridge that gap.  By encapsulating the core numerical algorithms in an
abbreviated code base that is open-source, computer scientists can analyze the
algorithms and easily make code-architectural changes to test performance
without compromising the production code values of Exnihilo_.

|Profugus| is **not** meant to be production software with respect to problem
analysis.  The computational kernels in |Profugus| are designed to analyze
*performance*, not *correctness*.  Nonetheless, users of |Profugus| can setup
and run problems with enough *real-world* features to be useful as
proof-of-concept for actual production work.

|Profugus| is also used as a test-bed mini-app for the ASCR projects XPRESS
and MCREX.

Profugus Development Team
=========================

|Profugus| contains the computational kernels from the Exnihilo_ code
base. The core Exnihilo_ development team consists of the following scientists
(listed alphabetically)

* Kevin Clarno <clarnokt@ornl.gov>
* Greg Davidson <davidsongg@ornl.gov>
* Tom Evans <evanstm@ornl.gov>
* Steven Hamilton <hamiltonsr@ornl.gov>
* Seth Johnson <johnsonsr@ornl.gov>
* Tara Pandya <pandyatm@ornl.gov>
* Stuart Slattery <slatterysr@ornl.gov>
* Rachel Slaybaugh <slaybaugh@berkeley.edu>

|Profugus| is developed and maintained by the following team:

* Tom Evans <evanstm@ornl.gov>
* Steven Hamilton <hamiltonsr@ornl.gov>
* Stuart Slattery <slatterysr@ornl.gov>

Profugus Packages
=================

|Profugus| constains the following code mini-apps:

**SPn**
  Mini-app of Simplified Spherical Harmonics (SPn) computational kernel.

Documentation on the mini-apps (methods and execution) are given in
:ref:`running_mini_app`. There are several support packages that are used by
the mini-apps.  These are:

**Utils**
  Utilities and testing framework used by other packages.

**Matprop**
  Data structures for storing material properties (cross-sections).

|Profugus| is designed to build and run with a minimum of dependencies.
However, there are some requirements.  The third-party software (TPLs)
necessary to build |Profugus| is all open-source and freely available.  The
TPLs for |Profugus| are listed in the following table:

.. table:: Profugus TPLs

+---------------------+-------------+-------------------------------------+
|         TPL         |  Required   |         Comments                    |
+=====================+=============+=====================================+
|Trilinos_            |Yes          |Trilinos_ is an open-source solvers  |
|                     |             |library.                             |
+---------------------+-------------+-------------------------------------+
|BLAS/LAPACK          |Yes          |Use vendor-specific implementation   |
|                     |             |when possible.  The default          |
|                     |             |implementations can be found at      |
|                     |             |netlib_.                             |
+---------------------+-------------+-------------------------------------+
|MPI_                 |No           |OpenMPI_ and MPICH_ are suggested    |
|                     |             |options. Vendor-provided options will|
|                     |             |also work.                           |
+---------------------+-------------+-------------------------------------+
|HDF5                 |No           |Parallel HDF5_ will allow output of  |
|                     |             |solution fiilds.  It is not needed   |
|                     |             |for general problem execution.       |
+---------------------+-------------+-------------------------------------+

.. _Exnihilo: denovo@email.ornl.gov
.. _Trilinos: http://trilinos.sandia.gov
.. _MPI: http://www.mcs.anl.gov/mpi
.. _MPICH: http://www.mcs.anl.gov/mpi
.. _OpenMPI: http://www.openmpi.org
.. _netlib: http://www.netlib.org
.. _HDF5: http://www.hdfgroup.org/HDF5
