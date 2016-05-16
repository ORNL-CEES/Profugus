//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/EigenvalueSolverBuilder.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Mon Feb 24 13:49:22 2014
 * \brief  EigenvalueSolverBuilder class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_solvers_EigenvalueSolverBuilder_hh
#define SPn_solvers_EigenvalueSolverBuilder_hh

#include "EigenvalueSolver.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class EigenvalueSolverBuilder
 * \brief Factory class for creating EigenvalueSolver.
 *
 * This class provides a generic interface for creating an eigensolver
 * involving an Epetra/Tpetra operator (or operators)
 *
 * \sa EigenvalueSolverBuilder.cc for detailed descriptions.
 */
/*!
 * \example solvers/test/tstEigenvalueSolverBuilder.cc
 *
 * Test of EigenvalueSolverBuilder.
 */
//===========================================================================//

template <class T>
class EigenvalueSolverBuilder
{
  public:
    //@{
    //! Typedefs.
    typedef typename T::MV                        MV;
    typedef typename T::OP                        OP;
    typedef EigenvalueSolver<T>                   EigenvalueSolver_t;
    typedef Teuchos::RCP<EigenvalueSolver_t>      RCP_EigenvalueSolver;
    typedef Teuchos::RCP<Teuchos::ParameterList>  RCP_ParameterList;
    //@}

    // Build standard eigenvalue solver
    static RCP_EigenvalueSolver build_solver(RCP_ParameterList db,
                                             Teuchos::RCP<OP> A);

    // Build generalized eigenvalue solver with preconditioner
    static RCP_EigenvalueSolver build_solver(
            RCP_ParameterList db,
            Teuchos::RCP<OP> A,
            Teuchos::RCP<OP> B,
            Teuchos::RCP<OP> P = Teuchos::null);
};

} // end namespace profugus

#endif // SPn_solvers_EigenvalueSolverBuilder_hh

//---------------------------------------------------------------------------//
//                 end of EigenvalueSolverBuilder.hh
//---------------------------------------------------------------------------//
