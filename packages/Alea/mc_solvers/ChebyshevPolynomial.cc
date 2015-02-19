//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ChebyshevPolynomial.cc
 * \author Steven Hamilton
 * \brief  ChebyshevPolynomial class definitions.
 */
//---------------------------------------------------------------------------//

#include <chrono>

#include "ChebyshevPolynomial.hh"

#include "AleaTypedefs.hh"
#include "PolynomialUtils.hh"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziMultiVecTraits.hpp"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Construct Chebyshev polynomial.
 */
//---------------------------------------------------------------------------//
ChebyshevPolynomial::ChebyshevPolynomial(Teuchos::RCP<const MATRIX> A,
        Teuchos::RCP<Teuchos::ParameterList> pl)
  : Polynomial(A,pl)
{
    TEUCHOS_ASSERT( b_A != Teuchos::null );

    if( b_verbosity >= LOW )
        std::cout << "Creating Chebyshev polynomial coefficients"
            << " of order " << b_m << std::endl;

    // If we're computing eigenvalues, do that now
    bool compute_eigenvalues =
        b_poly_pl->get<bool>("compute_eigenvalues",true);
    if( compute_eigenvalues )
    {
        computeEigenvalues();
    }
    // Otherwise read values from input pl
    else
    {
        TEUCHOS_ASSERT( b_poly_pl->isType<SCALAR>("lambda_min") );
        TEUCHOS_ASSERT( b_poly_pl->isType<SCALAR>("lambda_max") );
        d_lambda_min = b_poly_pl->get<SCALAR>("lambda_min");
        d_lambda_max = b_poly_pl->get<SCALAR>("lambda_max");
    }

    // Compute Chebyshev parameters alpha,beta such that
    //  the eigenvalues of (alpha*I + beta*A) are contained in [-1,1]
    SCALAR alpha = -(d_lambda_max+d_lambda_min)/(d_lambda_max-d_lambda_min);
    SCALAR beta  = 2/(d_lambda_max - d_lambda_min);

    // If a contraction is specified, modify alpha,beta accordingly
    // This modifies the mapping such that [lambda_min,lambda_max]
    //  gets mapped to [-c,c] rather than [-1,1]
    if( b_poly_pl->isType<SCALAR>("chebyshev_contraction") )
    {
        SCALAR contraction = b_poly_pl->get<SCALAR>("chebyshev_contraction");
        alpha *= contraction;
        beta *= contraction;
    }

    // We will create the polynomial in the "power" basis
    // The Chebyshev polynomials have a natural basis as alpha + beta*x,
    //  but we actually want the coefficients of p(x) s.t.
    //  I - x*p(x) = Tn(x).  In order to do this, we compute the Chebyshev
    //  coefficients in the natural basis, convert to the power basis,
    //  then convert to coefficients of p(x)
    b_native_basis = Teuchos::rcp( new PolynomialBasis("power") );

    // Get "raw" Chebyshev coefficients Tn(x)
    // We can view this as the coeffcients of the Chebyshev polynomial
    // w.r.t. the basis (alpha + beta*x)
    Teuchos::ArrayRCP<const SCALAR> cheby_coeffs =
        PolynomialUtils::getChebyshevCoefficients(b_m+1);
    TEUCHOS_ASSERT( cheby_coeffs.size() == (b_m+2) );

    // Compute the coefficients of Tn(alpha+beta*x)
    // This is accomplished by performing a change of basis from
    // (alpha + beta*x) to the standard basis of "x"
    PolynomialBasis cheby_basis("arbitrary");
    cheby_basis.setBasisCoefficients(alpha,beta);
    Teuchos::ArrayRCP<const SCALAR> tmp_coeffs =
        b_native_basis->transformBasis(cheby_coeffs,cheby_basis);
    TEUCHOS_ASSERT( tmp_coeffs.size() == (b_m+2) );

    // Coefficients of preconditioner are c[i] = -tmp[i+1]/tmp[0]
    b_coeffs.resize(b_m+1);
    for( int i=0; i<b_m+1; ++i )
    {
        b_coeffs[i] = -tmp_coeffs[i+1]/tmp_coeffs[0];
    }

    // Set values of alpha, beta in target basis
    // Default to alpha, beta used in Chebyshev mapping but allow override
    if( !b_poly_pl->isType<SCALAR>("polynomial_basis_alpha") )
    {
        b_poly_pl->set("polynomial_basis_alpha",alpha);
    }

    if( !b_poly_pl->isType<SCALAR>("polynomial_basis_alpha") )
    {
        b_poly_pl->set("polynomial_basis_beta",beta);
    }
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Compute smallest and largest magnitude eigenvalues
//---------------------------------------------------------------------------//
void ChebyshevPolynomial::computeEigenvalues()
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::milliseconds;
    high_resolution_clock::time_point time_start = high_resolution_clock::now();

    if( b_verbosity >= HIGH )
        std::cout << "Computing matrix eigenvalues" << std::endl;

    // Create eigenvector
    bool zero_out = true;
    Teuchos::RCP<MV> x( new MV(b_A->getDomainMap(),1,zero_out) );

    // Anasazi eigenproblem
    typedef Anasazi::BasicEigenproblem<SCALAR,MV,OP> Eigenproblem;
    Teuchos::RCP<Eigenproblem> problem( new Eigenproblem(b_A,x) );
    problem->setNEV(1);
    problem->setProblem();

    // Set basic parameters
    Teuchos::ParameterList anasazi_pl;
    int global_length = x->getGlobalLength();
    int num_blocks = std::min(20,global_length-2);
    anasazi_pl.set<int>("Num Blocks",num_blocks);
    anasazi_pl.set<SCALAR>("Convergence Tolerance",1.0e-8);

    // Use BlockKrylovSchur for now, can try out different solvers later
    typedef Anasazi::BlockKrylovSchurSolMgr<SCALAR,MV,OP> KrylovSchur;

    // Compute largest magnitude
    anasazi_pl.set<std::string>("Which","LM");
    Teuchos::RCP<KrylovSchur> solver( new KrylovSchur(problem,anasazi_pl) );

    solver->solve();

    SCALAR realpart = solver->getRitzValues()[0].realpart;
    SCALAR imagpart = solver->getRitzValues()[0].imagpart;
    d_lambda_max = realpart;
    if( b_verbosity >= HIGH )
    {
        if( SCALAR_TRAITS::magnitude(imagpart) < 1.0e-8 )
        {
            std::cout << "Largest magnitude eigenvalue: "
                      << realpart << std::endl;
        }
        else
        {
            std::cout << "Largest magnitude eigenvalue: "
                      << realpart << " +/- " << imagpart << "i" << std::endl;
        }
    }

    // Compute largest magnitude
    anasazi_pl.set<std::string>("Which","SM");
    solver = Teuchos::rcp( new KrylovSchur(problem,anasazi_pl) );

    solver->solve();

    realpart = solver->getRitzValues()[0].realpart;
    imagpart = solver->getRitzValues()[0].imagpart;
    d_lambda_min = realpart;
    if( b_verbosity >= HIGH )
    {
        if( SCALAR_TRAITS::magnitude(imagpart) < 1.0e-8 )
        {
            std::cout << "Smallest magnitude eigenvalue: "
                      << realpart << std::endl;
        }
        else
        {
            std::cout << "Smallest magnitude eigenvalue: "
                      << realpart << " +/- " << imagpart << "i" << std::endl;
        }
    }

    high_resolution_clock::time_point time_end = high_resolution_clock::now();
    if( b_verbosity >= HIGH )
    {
        std::cout << "Eigenvalue calculation took "
            << duration_cast<milliseconds>(time_end-time_start).count()
            << " milliseconds" << std::endl;
    }

}

} // namespace alea

