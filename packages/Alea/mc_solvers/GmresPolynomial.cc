//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GmresPolynomial.cc
 * \author Steven Hamilton
 * \brief  GmresPolynomial class definitions.
 */
//---------------------------------------------------------------------------//

#include <chrono>
#include <random>

#include "GmresPolynomial.hh"
#include "AleaTypedefs.hh"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"

#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziTpetraAdapter.hpp"

namespace alea
{

// Use Anasazi traits classes for easy manipulations
typedef Anasazi::MultiVecTraits<SCALAR,MV>       MVT;
typedef Anasazi::OperatorTraits<SCALAR,MV,OP>    OPT;
typedef Teuchos::SerialDenseMatrix<LO,SCALAR>    SDM;

//---------------------------------------------------------------------------//
/*!
 * \brief Construct GMRES polynomial.
 */
//---------------------------------------------------------------------------//
GmresPolynomial::GmresPolynomial(Teuchos::RCP<const MATRIX> A,
                 Teuchos::RCP<Teuchos::ParameterList> pl)
  : Polynomial(A,pl)
{
    TEUCHOS_ASSERT( b_A != Teuchos::null );

    std::string gmres_type =
        b_poly_pl->get<std::string>("gmres_type","qr");

    TEUCHOS_ASSERT( gmres_type == "arnoldi" ||
                    gmres_type == "fom"     ||
                    gmres_type == "normal"  ||
                    gmres_type == "qr" );

    // All GMRES polynomials are constructed in power basis
    b_native_basis = Teuchos::rcp( new PolynomialBasis("power") );

    if( gmres_type == "arnoldi" )
        buildGmresPolyFromArnoldi();
    else if( gmres_type == "fom" )
        buildFomPolyFromArnoldi();
    else if( gmres_type == "normal" )
        buildGmresPolyFromNormalEqns();
    else if( gmres_type == "qr" )
        buildGmresPolyFromQrDecomp();
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Build GMRES coefficients from Arnoldi decomposition
//---------------------------------------------------------------------------//
void GmresPolynomial::buildGmresPolyFromArnoldi()
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::milliseconds;
    high_resolution_clock::time_point time_start = high_resolution_clock::now();

    if( b_verbosity >= LOW )
        std::cout << "Creating GMRES polynomial coefficients"
            << " of order " << b_m << " from Arnoldi process" << std::endl;

    Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SCALAR> > arnoldi_mat =
        computeArnoldiMatrix(b_m+1);

    TEUCHOS_ASSERT( arnoldi_mat != Teuchos::null );

    Teuchos::ArrayRCP<const Anasazi::Value<SCALAR> > gmres_roots =
        computeGmresRoots(arnoldi_mat);

    TEUCHOS_ASSERT( gmres_roots.size() == b_m+1 );

    b_coeffs = computeCoeffsFromRoots(gmres_roots);

    high_resolution_clock::time_point time_end = high_resolution_clock::now();
    if( b_verbosity >= HIGH )
    {
        std::cout << "GMRES root calculation took "
            << duration_cast<milliseconds>(time_end-time_start).count()
            << " milliseconds" << std::endl;
    }

}

//---------------------------------------------------------------------------//
// Build FOM coefficients from Arnoldi decomposition
//---------------------------------------------------------------------------//
void GmresPolynomial::buildFomPolyFromArnoldi()
{
    if( b_verbosity >= LOW )
        std::cout << "Creating FOM polynomial coefficients"
            << " of order " << b_m << " from Arnoldi process" << std::endl;

    Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SCALAR> > arnoldi_mat =
        computeArnoldiMatrix(b_m+1);

    TEUCHOS_ASSERT( arnoldi_mat != Teuchos::null );

    Teuchos::ArrayRCP<const Anasazi::Value<SCALAR> > fom_roots =
        computeFomRoots(arnoldi_mat);

    TEUCHOS_ASSERT( fom_roots.size() == b_m+1 );

    b_coeffs = computeCoeffsFromRoots(fom_roots);
}

//---------------------------------------------------------------------------//
// Build GMRES coefficients from normal equations solution of (A*K)c=b
//---------------------------------------------------------------------------//
void GmresPolynomial::buildGmresPolyFromNormalEqns()
{

    Teuchos::RCP<const MV> b = computeInitialKrylovVector();

    // Get standard basis for Krylov subspace
    Teuchos::RCP<MV> K = computeKrylovStandardBasis(b_m+1, b);

    // Discard first vector (b)
    std::vector<LO> inds(b_m+1);
    for( LO i=0; i<b_m+1; ++i )
        inds[i]=i+1;
    Teuchos::RCP<const MV> AK = MVT::CloneView(*K,inds);

    Teuchos::RCP<SDM> AK_trans_AK( new SDM(b_m+1, b_m+1) );
    Teuchos::RCP<SDM> coeffs(      new SDM(b_m+1, 1) );

    MVT::MvTransMv(1.0,*AK,*AK,*AK_trans_AK);
    MVT::MvTransMv(1.0,*AK,*b,*coeffs);

    // Add regularization to diagonal
    SCALAR reg = b_poly_pl->get<SCALAR>("gmres_regularization",0.0);
    for( LO i=0; i<b_m+1; ++i )
    {
        (*AK_trans_AK)(i,i) += reg;
    }

    // Solve using LAPACK Cholesky decomp
    Teuchos::LAPACK<LO,SCALAR> lapack;
    int info;
    lapack.POSV('L',b_m+1,1,AK_trans_AK->values(),AK_trans_AK->stride(),
        coeffs->values(),coeffs->stride(),&info);
    TEUCHOS_ASSERT( info == 0 );

    // Copy coefficients into array in base class
    b_coeffs.resize(b_m+1);
    std::copy( coeffs->values(), coeffs->values()+b_m+1, b_coeffs.begin() );
}

//---------------------------------------------------------------------------//
// Build GMRES coefficients from QR solution of (K*A)c=b
//---------------------------------------------------------------------------//
void GmresPolynomial::buildGmresPolyFromQrDecomp()
{
    Teuchos::RCP<const MV> b = computeInitialKrylovVector();

    // Get standard basis for Krylov subspace
    Teuchos::RCP<MV> K = computeKrylovStandardBasis(b_m+1, b);

    // Discard first vector (b)
    std::vector<int> inds(b_m+1);
    for( int i=0; i<b_m+1; ++i )
        inds[i]=i+1;
    Teuchos::RCP<MV> AK = MVT::CloneViewNonConst(*K,inds);

    // Perform QR factorization of AK
    Anasazi::BasicOrthoManager<SCALAR,MV,OP> orthoman;
    Teuchos::RCP<SDM> R( new SDM(b_m+1,b_m+1) );
    int rank = orthoman.normalizeMat(*AK,R);
    TEUCHOS_ASSERT( rank == b_m+1 );

    // Compute Q^T * b
    SDM Q_trans_b(b_m+1,1);
    MVT::MvTransMv(1.0,*AK,*b,Q_trans_b);

    // R^{-1} (Q^T * b)
    // Solve using LAPACK triangular solve
    Teuchos::LAPACK<LO,SCALAR> lapack;
    int info;
    lapack.TRTRS('U','N','N',b_m+1,1,R->values(),R->stride(),
        Q_trans_b.values(),Q_trans_b.stride(),&info);
    TEUCHOS_ASSERT( info == 0 );

    // Copy coefficients into array in base class
    b_coeffs.resize(b_m+1);
    std::copy( Q_trans_b.values(), Q_trans_b.values()+b_m+1,
               b_coeffs.begin() );
}

//---------------------------------------------------------------------------//
// Compute random unit vector to seed Krylov subspace.
//---------------------------------------------------------------------------//
Teuchos::RCP<const MV> GmresPolynomial::computeInitialKrylovVector() const
{
    // Multivector to hold basis vectors
    Teuchos::RCP<MV> v( new MV(b_A->getDomainMap(),1) );

    // Start with random vector rather than constant, produces far better
    // convergence behavior, consistent with observations in
    // Joubert, SISC Vol. 15 (1994).
    std::vector<SCALAR> vec_norm(1);
    {
        std::mt19937_64 engine;
        bool reproducible = b_poly_pl->get("reproducible_random",false);
        if( reproducible )
        {
            // Ensure reproducibility in parallel by using same seed on
            // every process and discarding as many numbers as the global
            // index of the first local element
            // This is expensive and should only be done for small problems
            int seed = b_poly_pl->get("random_seed",12345);
            engine.seed(seed);

            int discard = v->getMap()->getGlobalElement(0);
            engine.discard(discard);
        }
        else
        {
            std::random_device rd;
            engine.seed(rd());
        }

        Teuchos::ArrayRCP<double> v_data = v->getDataNonConst(0);
        std::uniform_real_distribution<SCALAR> dist(0.0,1.0);
        for( int i=0; i<v_data.size(); ++i )
        {
            v_data[i] = dist(engine);
        }
    }
    MVT::MvNorm(*v,vec_norm);
    MVT::MvScale(*v,1.0/vec_norm[0]);

    /*
    // Use this for getting initial vector for MATLAB calculation of test data
    Teuchos::ArrayRCP<const SCALAR> data = v->getData(0);
    std::cout << "Initial Krylov vector" << std::endl;
    for( int i=0; i<data.size(); ++i )
    {
        std::cout << i << " " << std::setprecision(14) << data[i] << std::endl;
    }
    */

    return v;
}

//---------------------------------------------------------------------------//
// Convert roots of polynomial to coefficients
//---------------------------------------------------------------------------//
Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SCALAR> >
GmresPolynomial::computeArnoldiMatrix(int order) const
{
    TEUCHOS_ASSERT( b_A  != Teuchos::null );
    TEUCHOS_ASSERT( order > 0 );

    if( b_verbosity >= HIGH )
        std::cout << "Computing GMRES roots" << std::endl;

    // Do a quick Arnoldi iteration to compute Hessenberg matrix H

    // Multivector to hold basis vectors
    MV V(b_A->getDomainMap(),order+1);

    // Dense upper Hessenberg matrix
    bool zero_out = true;
    Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SCALAR> > H(
        new Teuchos::SerialDenseMatrix<LO,SCALAR>(order+1,order,zero_out) );

    // Set initial vector
    std::vector<SCALAR> dots(1), vec_norm(1);
    std::vector<int> v_ind(1), q_ind(1);
    v_ind[0] = 0;
    Teuchos::RCP<MV> v = MVT::CloneViewNonConst(V,v_ind);
    Teuchos::RCP<const MV> q;

    // Get initial vector
    Teuchos::RCP<const MV> b = computeInitialKrylovVector();
    MVT::SetBlock(*b,v_ind,*v);

    // Start Arnoldi iteration: Trefethen & Bau Algorithm 33.1
    for( int n=0; n<order; ++n )
    {
        // v = A*q_n
        v_ind[0] = n+1;
        v = MVT::CloneViewNonConst(V,v_ind);
        q_ind[0] = n;
        q = MVT::CloneView(V,q_ind);
        OPT::Apply(*b_A,*q,*v);


        for( int j=0; j<n+1; ++j )
        {
            // H_jn = q_j^* v
            q_ind[0] = j;
            q = MVT::CloneView(V,q_ind);
            MVT::MvDot(*v,*q,dots);
            (*H)(j,n) = dots[0];

            // v = v - H_jn * q_j
            MVT::MvAddMv(SCALAR_TRAITS::one(),*v,-(*H)(j,n),*q,*v);
        }

        // H_n+1,n = ||v||
        MVT::MvNorm(*v,vec_norm);
        (*H)(n+1,n) = vec_norm[0];

        // q_n+1 = v / ||v||
        MVT::MvScale(*v,1.0/vec_norm[0]);
    }

    return H;
}

//---------------------------------------------------------------------------//
// Compute GMRES polynomial roots from Arnoldi matrix
//---------------------------------------------------------------------------//
Teuchos::ArrayRCP<const Anasazi::Value<SCALAR> >
GmresPolynomial::computeGmresRoots(
    Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SCALAR> > H ) const
{
    LO order = H->numCols();
    TEUCHOS_ASSERT( H->numRows() == order+1 );

    // Compute eigenvalues of leading NxN block of H
    Teuchos::SerialDenseMatrix<LO,SCALAR> H_square(
        Teuchos::View,*H,order,order);

    // Matrix to hold (H^T * H)
    Teuchos::SerialDenseMatrix<LO,SCALAR> HtransH(order,order);
    int err = HtransH.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*H,*H,0.0);
    TEUCHOS_ASSERT( err == 0 );

    Teuchos::ArrayRCP<SCALAR> alpha_r(order), alpha_i(order), beta_r(order);
    LO worksize = 8*order;
    Teuchos::ArrayRCP<SCALAR> work(worksize);
    int info;
    Teuchos::LAPACK<LO,SCALAR> lapack;
    lapack.GGEV('N','N',order,HtransH.values(),HtransH.stride(),
        H_square.values(),H_square.stride(), alpha_r.get(),alpha_i.get(),
        beta_r.get(),0,order,0,order,work.get(),worksize,&info);
    TEUCHOS_ASSERT( info == 0 );

    Teuchos::ArrayRCP<Anasazi::Value<SCALAR> > roots(order);
    for( int i=0; i<order; ++i )
    {
        roots[i].realpart = alpha_r[i] / beta_r[i];
        roots[i].imagpart = alpha_i[i] / beta_r[i];
    }

    if( b_verbosity >= HIGH )
    {
        std::cout << "GMRES roots" << std::endl;
        for( int i=0; i<order; ++i )
        {
            std::cout << i << " " << roots[i].realpart << " "
                << roots[i].imagpart << std::endl;
        }
    }

    return roots;

}

//---------------------------------------------------------------------------//
// Compute FOM polynomial roots from Arnoldi matrix
//---------------------------------------------------------------------------//
Teuchos::ArrayRCP<const Anasazi::Value<SCALAR> >
GmresPolynomial::computeFomRoots(
    Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SCALAR> > H ) const
{
    LO order = H->numCols();
    TEUCHOS_ASSERT( H->numRows() == order+1 );

    // Compute eigenvalues of leading NxN block of H
    Teuchos::SerialDenseMatrix<LO,SCALAR> H_square(Teuchos::View,*H,order,order);

    Teuchos::ArrayRCP<SCALAR> real_parts(order), imag_parts(order);
    LO worksize = 6*order;
    Teuchos::ArrayRCP<SCALAR> work(worksize);
    int info;
    Teuchos::LAPACK<LO,SCALAR> lapack;
    lapack.HSEQR('E','N',order,1,order,H_square.values(),H_square.stride(),
        real_parts.get(),imag_parts.get(),0,order,work.get(),worksize,&info);

    TEUCHOS_ASSERT( info == 0 );

    Teuchos::ArrayRCP<Anasazi::Value<SCALAR> > roots(order);
    for( int i=0; i<order; ++i )
    {
        roots[i].realpart = real_parts[i];
        roots[i].imagpart = imag_parts[i];
    }

    if( b_verbosity >= HIGH )
    {
        std::cout << "FOM roots" << std::endl;
        for( int i=0; i<order; ++i )
        {
            std::cout << i << " " << real_parts[i] << " "
                << imag_parts[i] << std::endl;
        }
    }

    return roots;

}

//---------------------------------------------------------------------------//
// Convert roots of polynomial to coefficients
//---------------------------------------------------------------------------//
Teuchos::ArrayRCP<SCALAR> GmresPolynomial::computeCoeffsFromRoots(
    Teuchos::ArrayRCP<const Anasazi::Value<SCALAR> > roots) const
{
    Teuchos::ArrayRCP<SCALAR> coeffs(roots.size());

    // Expand polynomial from factored form p(x) = (x-l1)(x-l2)...
    //  into monomial basis p(x) = c0 + c1*x + c2*x^2 + ...
    Teuchos::ArrayRCP<SCALAR> c(1,1.0);
    int i=0;
    while( i<roots.size() )
    {
        SCALAR a = roots[i].realpart;
        SCALAR b = roots[i].imagpart;

        // Real roots: p_{n+1}(x) = (a-x)*p_{n}(x)
        if( b == 0.0 )
        {
            Teuchos::ArrayRCP<SCALAR> c1(c.size()+1,0.0);
            Teuchos::ArrayRCP<SCALAR> c2(c.size()+1,0.0);
            for( int j=0; j<i+1; ++j )
            {
                c1[j+1] = -c[j];
                c2[j] = a * c[j];
            }
            c.resize(c.size()+1,0.0);
            for( int j=0; j<c.size(); ++j )
            {
                c[j] = c1[j] + c2[j];
            }
            i++;
        }
        // Conjugate pair: p_{n+2}(x) = (x^2 - 2*a*x + (a^2+b^2))*p_{n}(x)
        else
        {
            Teuchos::ArrayRCP<SCALAR> c1(c.size()+2,0.0);
            Teuchos::ArrayRCP<SCALAR> c2(c.size()+2,0.0);
            Teuchos::ArrayRCP<SCALAR> c3(c.size()+2,0.0);
            for( int j=0; j<i+1; ++j )
            {
                c1[j+2] = c[j];
                c2[j+1] = -2.0*a*c[j];
                c3[j]   = (a*a + b*b) * c[j];
            }
            c.resize(c.size()+2,0.0);
            for( int j=0; j<c.size(); ++j )
            {
                c[j] = c1[j] + c2[j] + c3[j];
            }
            i+=2;
        }

    }
    TEUCHOS_ASSERT( c.size() == coeffs.size()+1 );

    for( int i=0; i<coeffs.size(); ++i )
    {
        coeffs[i] = -c[i+1] / c[0];
    }

    return coeffs;
}

//---------------------------------------------------------------------------//
// Compute "standard" basis for Krylov subspace: K = [b, A*b, ..., A^m*b]
//---------------------------------------------------------------------------//
Teuchos::RCP<MV> GmresPolynomial::computeKrylovStandardBasis(int order,
    Teuchos::RCP<const MV> b) const
{
    Teuchos::RCP<MV> K(new MV(b_A->getDomainMap(),order+1));

    std::vector<SCALAR> vec_norm(1);
    std::vector<int> v1_ind(1), v2_ind(1);
    Teuchos::RCP<const MV> v1;
    Teuchos::RCP<MV> v2;

    // Set first vector to b
    v2_ind[0] = 0;
    v2 = MVT::CloneViewNonConst(*K,v2_ind);
    MVT::SetBlock(*b,v2_ind,*v2);

    for( int i=0; i<order; ++i )
    {
        v1_ind[0] = i;
        v2_ind[0] = i+1;
        v1 = MVT::CloneView(*K,v1_ind);
        v2 = MVT::CloneViewNonConst(*K,v2_ind);
        OPT::Apply(*b_A,*v1,*v2);
    }

    return K;
}

} // namespace alea

