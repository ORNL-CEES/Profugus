//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/rtk/RTK_Array.hh
 * \author Thomas M. Evans
 * \date   Tue Dec 21 12:46:26 2010
 * \brief  RTK_Array class definition.
 * \note   Copyright (C) 2010 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_rtk_RTK_Array_hh
#define geometry_rtk_RTK_Array_hh

#include <vector>
#include <algorithm>
#include <cmath>
#include <ostream>
#include <utility>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "utils/Definitions.hh"
#include "utils/Vector_Lite.hh"
#include "utils/SP.hh"
#include "RTK_State.hh"

// Include Pin_Cell for specializations.
#include "RTK_Cell.hh"

namespace denovo
{

//===========================================================================//
/*!
 * \class RTK_Array
 * \brief Defines an array of Reactor ToolKit objects to assemble a reactor
 * core for MC tracking.
 *
 * Boundary conditions (reflecting faces) should only be set in the highest
 * level object.
 */
/*!
 * \example geometry/rtk/test/tstRTK_Array.cc
 *
 * Test of RTK_Array.
 */
//===========================================================================//

template<class T>
class RTK_Array
{
  public:
    //@{
    //! Typedefs.
    typedef T                      Object_t;
    typedef nemesis::SP<Object_t>  SP_Object;
    typedef std::vector<SP_Object> Object_Array;
    typedef def::Vec_Int           Vec_Int;
    typedef def::Vec_Dbl           Vec_Dbl;
    typedef RTK_State              Geo_State_t;
    typedef def::Space_Vector      Space_Vector;
    //@}

  private:
    // >>> DATA

    // Array sizes
    nemesis::Vector_Lite<int, 3> d_N;

    // Layout of objects in array.
    Vec_Int d_layout;

    // Array of core objects.
    Object_Array d_objects;

    // Array dimensions.
    Vec_Dbl d_x, d_y, d_z;

    // Coordinates in real space of lower-left coordinate of array.
    Space_Vector d_corner;

    // Length in each dimension.
    Space_Vector d_length;

    // Reflecting faces.
    Vec_Int d_reflect;

    // Number of cells in each array element.
    Vec_Int d_num_cells;
    Vec_Int d_Nc_offset;
    int     d_total_cells;

    // Outer vessel.
    nemesis::Vector_Lite<double, 2> d_r;
    nemesis::Vector_Lite<double, 2> d_origin;    // origin of vessel
    bool                   d_vessel;    // true if outer vessel defined
    int                    d_vessel_id; // material id for outer vessel

  public:
    // Constructor.
    RTK_Array(int Nx, int Ny, int Nz, int num_objects);

    // >>> SETUP FUNCTIONALITY

    //@{
    //! Return object id.
    int id(int i, int j, int k) const { return d_layout[index(i, j, k)]; }
    int& id(int i, int j, int k) { return d_layout[index(i, j, k)]; }
    //@}

    // Assign an object to the array.
    void assign_object(SP_Object object, int vid);

    // Finalize the construction of the array.
    void complete(double low_x, double low_y, double low_z);

    // Set reflecting faces.
    void set_reflecting(const Vec_Int &reflecting_faces);

    // Set outer vessel.
    void set_vessel(double r0, double r1, int vid);

    // >>> TRACKING FUNCTIONALITY

    // Initialize a state.
    void initialize(const Space_Vector &r, Geo_State_t &state) const;

    // Track to next boundary.
    void distance_to_boundary(const Space_Vector &r,
                              const Space_Vector &omega,
                              Geo_State_t &state) const;

    // Update a state at collision sites.
    void update_state(Geo_State_t &state) const;

    // Cross a surface.
    void cross_surface(const Space_Vector &r, Geo_State_t &state);

    // Find the object a point is in.
    int find_object(const Space_Vector &r, Geo_State_t &state) const;
    int find_object_on_boundary(const Space_Vector &r, int face, int face_type,
                                Geo_State_t &state) const;

    // >>> ACCESSORS

    // Return object.
    inline const Object_t& object(int i, int j, int k) const;

    // Return object.
    inline Object_t& object(int index) const;

    //! Number of objects.
    int num_objects() const { return d_objects.size(); }

    // Return the current material id.
    inline int matid(const Geo_State_t &state) const;

    //! Total number of cells.
    int num_cells() const { return d_total_cells; }

    // Return the current cell id.
    inline int cellid(const Geo_State_t &state) const;

    // Get extents of this geometry element in the parent reference frame
    inline void get_extents(Space_Vector &lower, Space_Vector &upper) const;

    //@{
    //! Array size.
    int size(int dimension) const { return d_N[dimension]; }
    int size() const { return d_N[0] * d_N[1] * d_N[2]; }
    //@}

    //! Array indexing.
    inline int index(int i, int j, int k) const;

    //! Pitch.
    double pitch(int dim) const { return d_length[dim]; }

    //! Height.
    double height() const { return d_length[2]; }

    //! Level of array.
    int level() const { return d_level; }

    // Calculate the level of this array.
    static int calc_level();

    //! Construction of array completed?
    bool completed() const { return d_completed; }

    //! Total number of regions in the array
    int num_regions() const;

    //! Bottom corner of array
    const Space_Vector& corner(){ return d_corner; }

    // set mapped cells for symmetry
    void set_mapped_cells(
        int this_offset, int map_offset,
        Vec_Int::iterator t_it, Vec_Int::iterator m_it);

    // Diagnostic output.
    void output(std::ostream &out) const;

    //! Query if vessel is defined.
    bool has_vessel() const { return d_vessel; }

  private:
    // >>> IMPLEMENTATION

    // Required to access private members of same class.
    template<class X> friend class RTK_Array;

    // Logical array coordinates.
    typedef nemesis::Vector_Lite<int, 3> Logical_Array;

    // Vector of integer pairs.
    typedef std::vector< std::pair<int, int> > Vec_Int_Pair;
    typedef Vec_Int_Pair::const_iterator       Vec_Int_Pair_Itr;

    // Get object.
    inline SP_Object object(const Geo_State_t &state) const;

    // Diagnostic output.
    void output(std::ostream &out, int level, int obj_id) const;

    // Count the cells in each level.
    int count_cells();
    inline int cell_count_dispatch(int i, int j, int k);

    // Add vessel to all objects in the array.
    void add_vessel(double R0, double R1, double xoff, double yoff, int id);

    // Add the vessel to an object in the array.
    void add_vessel_to_object(int i, int j, int k, double R0, double R1,
                              double xoff, double yoff, int vid);

    // Return the widths by dimension.
    double dx(int i) const { return d_x[i+1] - d_x[i]; }
    double dy(int j) const { return d_y[j+1] - d_y[j]; }
    double dz(int k) const { return d_z[k+1] - d_z[k]; }

    // Transform to object coordinate system.
    inline Space_Vector transform(const Space_Vector &r,
                                  const Geo_State_t &state) const;

    // Cross surface into next array element.
    inline void calc_high_face(Geo_State_t &state, int face_type,
                               int exiting_face);
    inline void calc_low_face(Geo_State_t &state, int face_type,
                              int exiting_face);

    // Determine boundary crossings at each level in the array.
    void determine_boundary_crossings(Geo_State_t &state);

    // Update the coordinates of each level in the array.
    void update_coordinates(const Space_Vector &r, Geo_State_t &state);

    // Calculate level.
    int d_level;

    // Completed flag for checks.
    bool d_completed;
};

} // end namespace denovo

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "RTK_Array.i.hh"

#endif // geometry_rtk_RTK_Array_hh

//---------------------------------------------------------------------------//
//              end of geometry/RTK_Array.hh
//---------------------------------------------------------------------------//
