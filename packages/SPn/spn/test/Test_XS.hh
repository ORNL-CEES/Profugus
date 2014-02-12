//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/Test_XS.hh
 * \author Thomas M. Evans
 * \date   Fri Oct 19 12:41:28 2012
 * \brief  Cross-sections for testing.
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_test_Test_XS_hh
#define spn_test_Test_XS_hh

#include <vector>

#include "Teuchos_RCP.hpp"

#include "xs/Mat_DB.hh"

//---------------------------------------------------------------------------//
// 12 Group Cross Sections
//---------------------------------------------------------------------------//

namespace twelve_grp
{

extern double T[12];
extern double S0[12][12];
extern double S1[12][12];
extern double S2[12][12];
extern double S3[12][12];

Teuchos::RCP<profugus::Mat_DB> make_mat(int Pn, int Nc);

}

//---------------------------------------------------------------------------//
// 3 Group Cross Sections
//---------------------------------------------------------------------------//

namespace three_grp
{

extern double T[3];
extern double S0[3][3];
extern double S1[3][3];
extern double S2[3][3];
extern double S3[3][3];

Teuchos::RCP<profugus::Mat_DB> make_mat(int Pn, int Nc);

Teuchos::RCP<profugus::Mat_DB> make_mat(int Pn, const std::vector<int> &matids,
                                        const std::vector<double> &f,
                                        const std::vector<int> &cell2mid);
}

//---------------------------------------------------------------------------//
// 1 Group Cross Sections
//---------------------------------------------------------------------------//

namespace one_grp
{

extern double T[1];
extern double S0[1][1];
extern double S1[1][1];
extern double S2[1][1];
extern double S3[1][1];

Teuchos::RCP<profugus::Mat_DB> make_mat(int Pn, int Nc);

Teuchos::RCP<profugus::Mat_DB> make_mat(int Pn, const std::vector<int> &matids,
                                        const std::vector<double> &f,
                                        const std::vector<int> &cell2mid);
}

//---------------------------------------------------------------------------//
// 2 Group Cross Sections
//---------------------------------------------------------------------------//

namespace two_grp
{

extern double T[2];
extern double S0[2][2];
extern double S1[2][2];
extern double S2[2][2];
extern double S3[2][2];

Teuchos::RCP<profugus::Mat_DB> make_mat(int Pn, int Nc);

Teuchos::RCP<profugus::Mat_DB> make_mat(int Pn, const std::vector<int> &matids,
                                        const std::vector<double> &f,
                                        const std::vector<int> &cell2mid);
}

#endif // spn_test_Test_XS_hh

//---------------------------------------------------------------------------//
//              end of spn/Test_XS.hh
//---------------------------------------------------------------------------//
