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

 > mpirun -np 4 ./xspn -i inf_med.xml

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

.. _Trilinos: http://trilinos.sandia.gov
