//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/LinAlgTypedefs.hh
 * \author Steven Hamilton
 * \date   Mon Jul 01 11:48:09 2013
 * \brief  Collection of typedefs for Tpetra classes.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_LinAlgTypedefs_hh
#define solvers_LinAlgTypedefs_hh

#include "Teuchos_RCP.hpp"

#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"
#include "Tpetra_Operator.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Kokkos_DefaultNode.hpp"

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_Matrix.hpp"

namespace profugus
{

enum LinAlgType {EPETRA, TPETRA, XPETRA};

template <LinAlgType T>
struct LinAlgTypedefs
{
};

template <>
struct LinAlgTypedefs<EPETRA>
{
    typedef double                    ST;
    typedef int                       LO;
    typedef int                       GO;
    typedef KokkosClassic::SerialNode NODE;
    typedef Epetra_MultiVector        MV;
    typedef Epetra_Operator           OP;
    typedef Epetra_Map                MAP;
    typedef Epetra_CrsMatrix          MATRIX;
    typedef Epetra_CrsGraph           GRAPH;
};

template <>
struct LinAlgTypedefs<TPETRA>
{
    typedef double                                      ST;
    typedef int                                         LO;
    typedef int                                         GO;
    typedef KokkosClassic::DefaultNode::DefaultNodeType NODE;
    typedef Tpetra::MultiVector<ST,LO,GO,NODE>          MV;
    typedef Tpetra::Operator<ST,LO,GO,NODE>             OP;
    typedef Tpetra::Map<LO,GO,NODE>                     MAP;
    typedef Tpetra::CrsMatrix<ST,LO,GO,NODE>            MATRIX;
    typedef Tpetra::RowGraph<LO,GO,NODE>                GRAPH;
};

template <>
struct LinAlgTypedefs<XPETRA>
{
    typedef double                                      ST;
    typedef int                                         LO;
    typedef int                                         GO;
    typedef KokkosClassic::DefaultNode::DefaultNodeType NODE;
    typedef Xpetra::MultiVector<ST,LO,GO,NODE>          MV;
    typedef Xpetra::CrsMatrix<ST,LO,GO,NODE>            CRS_MATRIX;
    typedef Xpetra::Matrix<ST,LO,GO,NODE>               MATRIX;
};

} // end namespace profugus

#endif // solvers_LinAlgTypedefs_hh

//---------------------------------------------------------------------------//
//              end of solvers/LinAlgTypedefs.hh
//---------------------------------------------------------------------------//
