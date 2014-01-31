//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/XS.hh
 * \author Thomas M. Evans
 * \date   Wed Jan 29 15:27:36 2014
 * \brief  XS class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef xs_XS_hh
#define xs_XS_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TwoDArray.hpp"

#include <set>
#include <vector>
#include "harness/DBC.hh"
#include "utils/Static_Map.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class XS
 * \brief Cross-section container class.
 */
/*!
 * \example xs/test/tstXS.cc
 *
 * Test of XS.
 */
//===========================================================================//

class XS
{
  public:
    //@{
    // Typedefs.
    typedef Teuchos::SerialDenseVector<int, double> Vector;
    typedef Teuchos::SerialDenseMatrix<int, double> Matrix;
    typedef Teuchos::RCP<Vector>                    RCP_Vector;
    typedef Teuchos::RCP<Matrix>                    RCP_Matrix;
    typedef Static_Map<int, RCP_Vector>             Hash_Vector;
    typedef Static_Map<int, RCP_Matrix>             Hash_Matrix;
    typedef Teuchos::Array<double>                  OneDArray;
    typedef Teuchos::TwoDArray<double>              TwoDArray;
    //@}

    //! Cross section types.
    enum XS_Types
    {
        TOTAL = 0,
        SIG_F,
        NU_SIG_F,
        CHI,
        END_XS_TYPES
    };

  private:
    // >>> DATA

    // Pn order of scattering.
    int d_pn;

    // Number of groups.
    int d_Ng;

    // Final number of materials.
    int d_Nm;

    // Hash table of totals.
    std::vector<Hash_Vector> d_totals;

    // Hash table of scattering (vector dimensioned over number of moments).
    std::vector<Hash_Matrix> d_scatter;

  public:
    XS();

    // >>> SETTING XS DATA

    // Set number of groups and Pn order that is stored.
    void set(int Pn_order, int num_groups);

    // Add 1-D cross sections to the database.
    void add(int matid, int type, const OneDArray &data);

    // Add scattering cross sections to the database.
    void add(int matid, int pn, const TwoDArray &data);

    // Complete assignment.
    void complete();

    // >>> ACCESSORS

    //! Pn order of data.
    int pn_order() const { return d_pn; }

    //! Number of groups in data.
    int num_groups() const { return d_Ng; }

    //! Number of materials in database (invalid until complete() called).
    int num_mat() const { return d_Nm; }

  private:
    // >>> IMPLEMENTATION

    typedef std::set<int>           Set_Int;
    typedef std::vector<Set_Int>    Vec_Set;
    typedef Set_Int::const_iterator set_iter;

    Vec_Set d_inst_totals;
    Vec_Set d_inst_scat;
};

} // end namespace profugus

//---------------------------------------------------------------------------//

#include "XS.i.hh"

#endif // xs_XS_hh

//---------------------------------------------------------------------------//
//                 end of XS.hh
//---------------------------------------------------------------------------//
