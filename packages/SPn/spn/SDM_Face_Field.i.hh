//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/SDM_Face_Field.i.hh
 * \author Thomas M. Evans
 * \date   Tue Oct 30 22:38:08 2012
 * \brief  Member definitions of class SDM_Face_Field.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_SDM_Face_Field_i_hh
#define spn_SDM_Face_Field_i_hh

namespace profugus
{

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Get a view of the matrix from the field.
 *
 * This returns a view of a SerialDenseMatrix from the specified location.
 */
SDM_Face_Field::Serial_Matrix SDM_Face_Field::view(int abscissa,
                                                   int ordinate)
{
    Serial_Matrix m(
        Teuchos::View, convert(abscissa, ordinate), d_M, d_M, d_M);
    ENSURE(m.numRows() == d_M);
    ENSURE(m.numCols() == d_M);
    return m;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Insert a \f$n\times n\f$ matrix into the field.
 *
 * This returns a view of a SerialDenseMatrix from the specified location.
 */
void SDM_Face_Field::insert(int                  abscissa,
                            int                  ordinate,
                            const Serial_Matrix &m)
{
    REQUIRE(m.numRows() == d_M);
    REQUIRE(m.numCols() == d_M);
    REQUIRE(m.stride() == d_M);

    d_blas.COPY(d_size, m.values(), 1, convert(abscissa, ordinate), 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief \f$ O(1)\f$ swap with another SDM_Face_Field.
 *
 * Given SDM_Face_Field a,b calling swap as follows:
 * \code
     a.swap(b);
 * \endcode
 * will result in a having b's values and b having a's values.  Both fields
 * must be on the same mesh (block) and face type (XYZ).
 *
 * \param b face field which will take a's values with a taking b's values
 */
void SDM_Face_Field::swap(SDM_Face_Field &b)
{
    REQUIRE(b.d_face == d_face);
    REQUIRE(b.d_N_abscissa == d_N_abscissa);
    REQUIRE(b.d_N_ordinate == d_N_ordinate);
    REQUIRE(b.d_N == d_N);
    REQUIRE(b.d_M == d_M);
    REQUIRE(b.d_size == d_size);

    d_field.swap(b.d_field);

    ENSURE(b.d_field.size() == d_field.size());
    ENSURE(static_cast<int>(d_field.size()) == d_N * d_size);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Fast copy of a block face field.
 *
 * An example:
 * \code
     SDM_Face_Field x(mesh, X, 4);
     SDM_Face_Field y(mesh, X, 4);

     // copy x data into y
     y.fast_copy(x);
 * \endcode
 *
 * The face types must be the same, ie. X-face to X-face, etc.
 *
 * \param b SDM_Face_Field that is copied into the calling object
 */
void SDM_Face_Field::fast_copy(const SDM_Face_Field &b)
{
    REQUIRE(b.d_face == d_face);
    REQUIRE(b.d_N_abscissa == d_N_abscissa);
    REQUIRE(b.d_N_ordinate == d_N_ordinate);
    REQUIRE(b.d_N == d_N);
    REQUIRE(b.d_M == d_M);
    REQUIRE(b.d_size == d_size);

    d_blas.COPY(d_N * d_size, &b.d_field[0], 1, &d_field[0], 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Fast copy of a block face field from a data stream.
 *
 * An example:
 * \code
     SDM_Face_Field x(mesh, X, 4);

     vector<double> values;

     // fast copy values into x
     x.fast_copy(&values[0], &values[0] + values.size());
 * \endcode
 *
 * The face types must be the same, ie. X-face to X-face, etc.  The values
 * must be correcly sized.
 */
void SDM_Face_Field::fast_copy(const_pointer begin,
                               const_pointer end)
{
    REQUIRE(end - begin == d_size * d_N);

    d_blas.COPY(d_N * d_size, begin, 1, &d_field[0], 1);
}

//---------------------------------------------------------------------------//
//@{
/*!
 * \brief Get a pointer to the continuous field data.
 */
SDM_Face_Field::const_pointer SDM_Face_Field::data_pointer() const
{
    REQUIRE(!d_field.empty());
    return &d_field[0];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Mutable version of data_pointer().
 */
SDM_Face_Field::pointer SDM_Face_Field::data_pointer()
{
    REQUIRE(!d_field.empty());
    return &d_field[0];
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
//@{
/*!
 * \brief Convert (abscissa, ordinate) to a starting pointer location of
 * block-matrix data.
 */
SDM_Face_Field::const_pointer SDM_Face_Field::convert(int abscissa,
                                                      int ordinate) const
{
    REQUIRE(abscissa < d_N_abscissa);
    REQUIRE(ordinate < d_N_ordinate);

    ENSURE(abscissa + ordinate * d_N_abscissa < d_N);
    ENSURE((abscissa + ordinate * d_N_abscissa) * d_size
            <= d_size * (d_N - 1));
    return &d_field[0] + (abscissa + ordinate * d_N_abscissa) * d_size;
}

SDM_Face_Field::pointer SDM_Face_Field::convert(int abscissa,
                                                int ordinate)
{
    REQUIRE(abscissa < d_N_abscissa);
    REQUIRE(ordinate < d_N_ordinate);

    ENSURE(abscissa + ordinate * d_N_abscissa < d_N);
    ENSURE((abscissa + ordinate * d_N_abscissa) * d_size
            <= d_size * (d_N - 1));
    return &d_field[0] + (abscissa + ordinate * d_N_abscissa) * d_size;
}
//@}

} // end namespace profugus

#endif // spn_SDM_Face_Field_i_hh

//---------------------------------------------------------------------------//
//              end of spn/SDM_Face_Field.i.hh
//---------------------------------------------------------------------------//
