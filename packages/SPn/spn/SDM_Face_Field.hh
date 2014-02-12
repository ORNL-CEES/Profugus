//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/SDM_Face_Field.hh
 * \author Thomas M. Evans
 * \date   Tue Oct 30 17:05:49 2012
 * \brief  SDM_Face_Field class definition.
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_SDM_Face_Field_hh
#define spn_SDM_Face_Field_hh

#include <vector>

#include <SPn/config.h>
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_BLAS.hpp"

#include "harness/DBC.hh"
#include "Mesh.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class SDM_Face_Field
 * \brief Face-field that continously stores SerialDenseMatrices.
 *
 * Data is stored in (abscissa, ordinate) pairs. To access the beginning of
 * the continuous data in the field (for communication purposes), use
 * data_pointer() and data_size().  To get an iterator that walks member-wise
 * through the continuous data, use begin_data() and end_data().
 *
 * The use-case for this class is:
 * \code

   // defines a field on an X face for 3x3 block matrices
   SDM_Face_Field field(mesh, X, 3);

   SerialDenseMatrix<int, double> D1;
   // fill D1

   // add D1 to the field
   field.insert(0, 0, D1);

   // communicate the field
   send(field.data_pointer(), field.data_size(), 1);

   // ...

   // receive the field
   recv(field.data_pointer(), field.data_size(), 0);

   // get views to the matrices
   SerialDenseMatrix<int, double> D1 = field.view(0, 0);

   // use D1 ...

   \endcode
 *
 * The key here is that the matrices that are accessed using view() are \b
 * views only.  Do not let the field go out of scope while using the view;
 * RESULTS WILL BE UNDEFINED.
 */
/*!
 * \example spn/test/tstSDM_Face_Field.cc
 *
 * Test of SDM_Face_Field.
 */
//===========================================================================//

class SDM_Face_Field
{
  public:
    //@{
    //! Typedefs.
    typedef Teuchos::SerialDenseMatrix<int, double> Serial_Matrix;
    typedef std::vector<double>                     Vec_Matrix;
    typedef Mesh                                    Mesh_t;
    typedef const double *                          const_pointer;
    typedef double *                                pointer;
    //@}

  private:
    // >>> DATA

    // Field of SerialDenseMatrix data.
    Vec_Matrix d_field;

    // Face type.
    int d_face;

    // Number of cells in (abscissa, ordinate) and total.
    int d_N_abscissa;
    int d_N_ordinate;
    int d_N;

    // Rank and size of each MxM block matrix in the vector.
    int d_M;
    int d_size;

  public:
    // Constructor.
    SDM_Face_Field(const Mesh_t &mesh, int face, int M);

    // Get a view of the matrix.
    inline Serial_Matrix view(int abscissa, int ordinate);

    // Insert a matrix into the field.
    inline void insert(int abscissa, int ordinate, const Serial_Matrix &m);

    // Fast copy of data from another equivalent face-field.
    inline void fast_copy(const SDM_Face_Field &b);

    // Fast copy from a data stream of values.
    inline void fast_copy(const_pointer begin, const_pointer end);

    //! Number of elements (block-matrices) in field.
    int size() const { return d_N; }

    //! Size of abscissa.
    int abscissa() const { return d_N_abscissa; }

    //! Size of ordinate.
    int ordinate() const { return d_N_ordinate; }

    // Swap with another SDM_Face_Field.
    inline void swap(SDM_Face_Field &b);

    // Return the face for this field.
    int face() const { return d_face; }

    // Get a constant pointer to the data for continuous accessing.
    inline const_pointer data_pointer() const;

    // Get a mutable pointer to the data for continuous accessing (careful).
    inline pointer data_pointer();

    // Get size of continuous data in the field.
    int data_size() const { return d_size * d_N; }

    //@{
    //! Data pointer iterators.
    const_pointer begin_data() const { return data_pointer(); }
    const_pointer end_data() const { return data_pointer() + data_size(); }
    pointer begin_data() { return data_pointer(); }
    pointer end_data() { return data_pointer() + data_size(); }
    //@}

  private:
    // >>> IMPLEMENTATION

    // Convert abcissa and ordinate to a pointer location in the field.
    inline const_pointer convert(int abscissa, int ordinate) const;
    inline pointer convert(int abscissa, int ordinate);

    // Teuchos BLAS interface.
    Teuchos::BLAS<int, double> d_blas;
};

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "SDM_Face_Field.i.hh"

#endif // spn_SDM_Face_Field_hh

//---------------------------------------------------------------------------//
//              end of spn/SDM_Face_Field.hh
//---------------------------------------------------------------------------//
