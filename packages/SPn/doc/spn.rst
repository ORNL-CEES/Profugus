Running the SPn Mini-Application
================================

.. highlight:: sh

The ``SPn`` mini-app is compiled into the executable :program:`xspn` that has
the following options:

.. program:: xspn

.. option:: -i

   Load the subsequent xml input file.

The ``SPn`` mini-app solves the :math:`SP_N` form of the neutron transport
equation.  The equations it solves are described in :download:`in this
technical note <spn_doc/spn.pdf>`.

Inputs are xml-driven.  To run the application do the following::

 > mpirun -np 1 ./xspn -i inf_med.xml

Infinite-Medium Problem Example
-------------------------------

In this case :file:`inf_med.xml` is the xml input file.  We can walk through
this case; it is one of the example inputs in :file:`packages/SPn/examples`:

.. literalinclude:: spn_examples/inf_med.xml
   :language: xml

This is an input of a 1-group infinite-medium problem that has an analytical
solution.  Integrating the Boltzmann transport equation over energy and angle
gives the following balance equation:

.. math::
   :label: balance_eq

   \nabla\cdot\mathbf{J} + (\sigma_{\text{t}}-\sigma_{\text{s}})\phi
    = Q

In an infinite medium (all reflecting boundaries with a uniform source and
material everywhere), the divergence of the current vanishes and the solution
for :math:`\phi` is simply

.. math::
   :label: inf_med_sol

   \phi = \frac{Q}{(\sigma_{\text{t}}-\sigma_{\text{s}})}

The geometric model for ``SPn`` is a simplified input for nuclear reactor
problems (although it can be used to describe other types of problems). In
this application area radial dimensions are *xy** and axial dimensions are
*z*.  A core is generally defined radially by 2D maps of assemblies and
pin-cells and axially by 1D dimensions.  In this model a **CORE** is the
outermost object.  A **CORE** contains **ASSEMBLIES**.  The **MATERIALS** are
referenced by the **ASSEMBLIES**.  For fixed-source problems a **SOURCE** may
be defined at the **CORE** level.  Finally, a set of problem parameters that
dictates how the code is run is defined in **PROBLEM**.  These are summarized
below:

**CORE**
  Describes the outer geometric configuration of the problem.  A core defines
  the axial (*z*-axis) definitions and a 2D assembly map.

**ASSEMBLIES**
  Individual assembly types are defined in this block.  Each assembly is an
  :math:`N\times M` array of *pin-cells*.  Each pin cell contains a single
  material, but can be descretized by a finer computational mesh grid.

**SOURCE**
  Describes fixed sources.  When this option is present, ``SPn`` is run in
  fixed-source mode; otherwise, an eigenvalue problem is run.

**MATERIAL**
  Specifies the xml-file containing the material cross sections and the list
  of materials referenced in the **ASSEMBLIES** block.

**PROBLEM**
  Determines all problem solution parameters including computational mesh
  discretization, SPN-order of the problem (1, 3, 5, 7), Pn-order of the
  scattering expansion, boundary conditions, and solver parameters.  Note that
  solver parameters can be passed to directly to Trilinos_ solvers by
  describing the appropriate parameters as shown below.

Details on how the problem parameters are used to set up a problem can be
ascertained by perusing the well-documented source code in :download:`the
Problem_Builder class <../../packages/SPn/driver/Problem_Builder.hh>`
(this is a mini-app after all, so we do not consider it unreasonable to look
at code!).

.. note::

   All of the 2D arrays in the xml-input are COLUMN-MAJOR (FORTRAN) ordered.
   In other words, the *i*-indices of an array cycle fastest. Thus, 2D arrays
   are stored internally as arrays with dimension ``[N_j][N_i]``, reflecting
   the fact that C arrays are stored ROW-MAJOR.

   Conversely, the scattering cross sections in the cross section xml files
   are ROW-MAJOR ordered.  Thus, the 2D arrays have dimensions ``[g][g']``.

   This will become clean in the examples.

.. highlight:: xml

Let's us now step through the input step-by-step. The **CORE** block for this
problem defines the following parameters::

 <Parameter name="axial list" type="Array(string)" value="{core}"/>
 <Parameter name="axial height" type="Array(double)" value="{ 1.00000e+01}"/>
 <Parameter name="core" type="TwoDArray(int)" value="1x1:{0}"/>

This first parameter defines a list of axial core maps.  The second parameter
defines the height of each axial level and should have the same number of
entries as the ``axial list`` parameter.  The final set of parameters are the
core maps listed in hte ``axial list`` array.  In this case there is a single
axial core map named ``core`` (the names can be unique).  The ``core`` map
specifies a core containing a single assembly.  This indices in the core map
are ordered :math:`[0,\ldots,N)` and refer to assemblies defined in the
**ASSEMBLIES** block ``assembly list`` parameter. This brings us to the
parameters defined in the **ASSEMBLIES** block::

 <Parameter name="pin pitch" type="double" value="10.0"/>
 <Parameter name="assembly list" type="Array(string)" value="{assembly}"/>
 <Parameter name="assembly" type="TwoDArray(int)" value="1x1:{0}"/>

The ``pin pitch`` parameter gives the width of each (square) pin-cell.  The
``assembly list`` parameter gives the list of assembly types.  As in the
**CORE** block, the final set of parameters are the 2D assembly maps for each
assembly defined in the ``assembly list``.  Every assembly must have the same
dimensions, but they can be composed of different materials.  Each assembly
map refers to a material defined in the **MATERIAL** block ``mat list``
parameter and is ordered :math:`[0,\ldots,N)`.  So, this block defines a
single assembly type that is :math:`10\times 10` cm in radial width that
contains a single pin-cell with material 0. Examining the 2 parameters in the
**MATERIAL** block::

 <Parameter name="xs library" type="string" value="xs_1G.xml"/>
 <Parameter name="mat list" type="Array(string)" value="{scatterer}"/>

We see from the ``mat list`` parameter that material 0 corresponds to the
material ``scatterer``.  This name refers to a material defined in the
cross-section file ``xs_1G.xml`` that is loaded based on the value of the ``xs
library`` parameter. The cross section file, :file:`xs_1G.xml`, containst the
following data:

.. literalinclude:: spn_examples/xs_1G.xml
   :language: xml

Thus, the geometric description of this problem is a :math:`10\times 10\times
10` cm box containing the material ``scatterer`` that has a total cross
section of 1.0 and a scattering cross section of 0.9 (resulting in an
aborption/removal cross section of 0.1).

To make this a fixed-source problem, the **SOURCE** block must be present::

 <Parameter name="source list" type="Array(string)" value="{uniform}"/>
 <Parameter name="source map" type="TwoDArray(int)" value="1x1:{0}"/>
 <Parameter name="axial source" type="Array(double)" value="{ 1.00000e+00}"/>
 <ParameterList name="uniform">
   <Parameter name="strength" type="TwoDArray(double)" value="1x1:{ 1.00000e+00}"/>
   <Parameter name="shape" type="Array(double)" value="{ 1.00000e+00}"/>
 </ParameterList>

Here, the ``source list`` parameter gives the list of sources. The ``source
map`` parameter has the same dimensions as the core maps and gives the source
index that resides in each core location.  A value of -1 indicates that there
is no source in a given location.  A value in the range :math:`[0,\ldots,N)`
refers to the source indicated by the ``source list``.  The ``axial source``
parameter has the same number of entries as the ``axial list`` and ``axial
height`` parameters in the **CORE** block.  The ``axial source`` is a
multiplier that is applied to the source at each axial level.  Each radial
source is defined in its own sublist.  In this case, there is a single source,
``uniform``.  Each source sublist contains a ``strength`` and ``shape``
parameter.  The ``strength`` parameter gives a :math:`N\times M` 2D axial
array of the source strength in each pin-cell.  The ``shape`` is dimensioned
by the number of groups in the ``xs library`` parameter ``num groups``.  It
gives the energy-spectral shape of the source.  Thus, in any pin-cell, the
group source is given by the product of the shape and strength:

.. math::
   :label: source

   q_{\text{ext}\,i,j}^g = \text{shape[g]} \times \text{strength[j][i]}

And, as described above, the ``strength`` map is COLUMN-MAJOR ordered in the
xml file.

Finally, the **PROBLEM** database provides entries that control the problem
setup and solver.  Many of these entries will be set by default, although all
can be overrided.  In this example, we only set a few parameters::

 <Parameter name="radial mesh" type="int" value="10"/>
 <Parameter name="axial mesh" type="Array(int)" value="{10}"/>
 <Parameter name="symmetry" type="string" value="full"/>
 <Parameter name="Pn_order" type="int" value="0"/>
 <Parameter name="SPn_order" type="int" value="1"/>
 <Parameter name="problem_name" type="string" value="inf_med"/>

Here the ``radial mesh`` indicates that each pin-cell is meshed
:math:`10\times 10`.  The ``axial mesh`` value is an array giving the number
of computational mesh cells in each axial level.  The ``symmetry`` parameter
tells that mesh generator to use the full problem description and not apply
any symmetry conditions.  The ``Pn_order`` parameter gives the scattering
order for the solver; it must be less than or equal to the ``pn order``
specified in the cross section library file.  The ``SPn_order`` gives the
order of the :math:`SP_N` approximation.  Finally, the ``problem_name`` gives
the base name that will be used for all output files. We have not specified
any specific solver options, so the defaults will be used.

The final specification of this problem is :math:`10\times 10\times 10` cm box
with computational mesh cells with dimension :math:`1\times 1\times 1` cm
resulting in 1000 total cells (:math:`10\times 10\times 10`).  There is a
uniform 1 particle/cc source throughout the box.  The box has a uniform
material of ``scatterer`` as defined in the :file:`xs_1G.xml` file.  The
solver will run a ``SP1`` calculation with ``P0`` scattering.

.. note::

   The only symmetry option currently supported is ``full`; however, ``qtr``
   symmetry will be added for 1/4 symmetry.

The outputs for this problem are contained in the :file:`examples` directory.
Automatically, the code will output a final problem xml file so that the user
can see what defaults were added.  In this case, the output xml file is stored
in :file:`inf_med_db.xml`:

.. literalinclude:: spn_examples/inf_med_db.xml
   :language: xml

Perusing this database, we see the default solver options that were set. These
can be overriden by adding a ``Stratimikos`` parameterlist to the
input xml. Other items of note are::

 <Parameter  name="num_blocks_i" type="int" value="1"/>
 <Parameter  name="num_blocks_j" type="int" value="1"/>
 <Parameter  name="g_first" type="int" value="0"/>
 <Parameter  name="g_last" type="int" value="0"/>

The ``num_blocks`` parameters allow the user to specify a parallel
decomposition.  The total number of processes is:

.. math::

   N_p = N_i \times N_j

.. note::

   The ``SPn`` mini-app currenly only supports *xy* decompositions. This is an
   historical requirement so that the ``SPn`` matches a similar decomposition
   in the |Exnihilo| production code for another physics model.  It is not a
   *true* restriction.

The ``g_first`` and ``g_last`` parameters allow the user to specify a range of
groups to run over.  For example, if a 23-group cross section set is loaded,
but the user is only interested in the solution in groups 0 and 1 set::

 <Parameter  name="g_first" type="int" value="0"/>
 <Parameter  name="g_last" type="int" value="1"/>

If HDF5 is available, the group-wise fluxes are output in the file
:file:`inf_med_output.h5`.  From :eq:`inf_med_sol` the solution to this
problem should be :math:`\phi = 10.0` everywhere.  Using :program:`h5dump` on
the output file yields:

.. code-block:: sh

 HDF5 "SPn_output.h5" {
 GROUP "/" {
    GROUP "fluxes" {
       DATASET "group_0" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 10, 10, 10 ) / ( 10, 10, 10 ) }
          DATA {
          (0,0,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (0,1,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (0,2,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (0,3,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (0,4,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (0,5,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (0,6,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (0,7,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (0,8,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (0,9,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (1,0,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (1,1,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (1,2,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (1,3,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (1,4,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (1,5,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (1,6,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (1,7,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (1,8,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (1,9,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (2,0,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (2,1,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (2,2,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (2,3,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (2,4,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (2,5,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (2,6,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (2,7,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (2,8,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (2,9,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (3,0,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (3,1,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (3,2,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (3,3,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (3,4,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (3,5,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (3,6,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (3,7,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (3,8,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (3,9,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (4,0,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (4,1,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (4,2,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (4,3,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (4,4,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (4,5,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (4,6,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (4,7,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (4,8,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (4,9,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (5,0,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (5,1,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (5,2,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (5,3,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (5,4,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (5,5,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (5,6,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (5,7,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (5,8,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (5,9,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (6,0,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (6,1,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (6,2,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (6,3,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (6,4,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (6,5,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (6,6,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (6,7,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (6,8,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (6,9,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (7,0,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (7,1,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (7,2,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (7,3,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (7,4,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (7,5,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (7,6,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (7,7,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (7,8,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (7,9,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (8,0,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (8,1,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (8,2,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (8,3,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (8,4,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (8,5,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (8,6,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (8,7,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (8,8,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (8,9,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (9,0,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (9,1,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (9,2,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (9,3,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (9,4,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (9,5,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (9,6,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (9,7,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (9,8,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
          (9,9,0): 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
          }
          ATTRIBUTE "data_order" {
             DATATYPE  H5T_STD_I32LE
             DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
             DATA {
             (0): 0
             }
          }
       }
    }
 }
 }

Thus, the correct solution is attained for this problem.

.. _Trilinos: http://trilinos.sandia.gov
