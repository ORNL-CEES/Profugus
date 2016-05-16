//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/LinAlgTypedefs.hh
 * \author Steven Hamilton
 * \date   Mon Jul 01 11:48:09 2013
 * \brief  Collection of typedefs for Tpetra classes.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_solvers_LinAlgTypedefs_hh
#define SPn_solvers_LinAlgTypedefs_hh

#include <SPn/config.h>

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

#ifdef USE_MUELU
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_Matrix.hpp"
#endif

namespace profugus
{

struct EpetraTypes
{
    typedef double                                      ST;
    typedef int                                         LO;
    typedef int                                         GO;
    typedef KokkosClassic::DefaultNode::DefaultNodeType NODE;
    typedef Epetra_MultiVector                          MV;
    typedef Epetra_Vector                               VECTOR;
    typedef Epetra_Operator                             OP;
    typedef Epetra_Map                                  MAP;
    typedef Epetra_CrsMatrix                            MATRIX;
    typedef Epetra_CrsGraph                             GRAPH;
};

struct TpetraTypes
{
    typedef double                                      ST;
    typedef int                                         LO;
    typedef int                                         GO;
    typedef KokkosClassic::DefaultNode::DefaultNodeType NODE;
    typedef Tpetra::MultiVector<ST,LO,GO,NODE>          MV;
    typedef Tpetra::Vector<ST,LO,GO,NODE>               VECTOR;
    typedef Tpetra::Operator<ST,LO,GO,NODE>             OP;
    typedef Tpetra::Map<LO,GO,NODE>                     MAP;
    typedef Tpetra::CrsMatrix<ST,LO,GO,NODE>            MATRIX;
    typedef Tpetra::CrsGraph<LO,GO,NODE>                GRAPH;
};

#ifdef USE_MUELU
struct XpetraTypes
{
    typedef double                                      ST;
    typedef int                                         LO;
    typedef int                                         GO;
    typedef KokkosClassic::DefaultNode::DefaultNodeType NODE;
    typedef Xpetra::MultiVector<ST,LO,GO,NODE>          MV;
    typedef Xpetra::CrsMatrix<ST,LO,GO,NODE>            CRS_MATRIX;
    typedef Xpetra::Matrix<ST,LO,GO,NODE>               MATRIX;
};
#endif

} // end namespace profugus

#endif // SPn_solvers_LinAlgTypedefs_hh

//---------------------------------------------------------------------------//
//              end of solvers/LinAlgTypedefs.hh
//---------------------------------------------------------------------------//
