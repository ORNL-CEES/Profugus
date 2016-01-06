//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/test/tstTeuchos.cc
 * \author Thomas M. Evans
 * \date   Wed Jan 29 15:27:12 2014
 * \brief  Test teuchos capabilities.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <set>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TwoDArray.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

//---------------------------------------------------------------------------//
// HELPERS
//---------------------------------------------------------------------------//

class Foo
{
  private:
    int d_x;

  public:
    Foo() {/*...*/}
    int  x() const { return d_x; }
    void set_x(int x) { d_x = x; }
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class Teuchos_Test : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture
    typedef Teuchos::RCP<Foo>                       RCP_Foo;
    typedef Teuchos::ParameterList                  ParameterList;
    typedef Teuchos::RCP<ParameterList>             RCP_ParameterList;
    typedef Teuchos::Comm<int>                      Comm;
    typedef Teuchos::RCP<const Comm>                RCP_Comm;
    typedef std::vector<double>                     Vec_Dbl;
    typedef std::vector<int>                        Vec_Int;
    typedef Teuchos::Array<double>                  Array_Dbl;
    typedef Teuchos::Array<int>                     Array_Int;
    typedef Teuchos::TwoDArray<double>              TwoDArray;
    typedef Teuchos::SerialDenseVector<int, double> Vector;
    typedef Teuchos::SerialDenseMatrix<int, double> Matrix;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();
        comm  = Teuchos::DefaultComm<int>::getComm();
    }

  protected:
    // >>> Data that get re-initialized between tests

    int node, nodes;
    RCP_Comm comm;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Teuchos_Test, RCP)
{
    RCP_Foo f1 = Teuchos::rcp(new Foo);
    RCP_Foo f2;
    RCP_Foo f3(new Foo);

    EXPECT_FALSE(f1.is_null());
    EXPECT_TRUE(f2.is_null());
    EXPECT_FALSE(f3.is_null());
}

//---------------------------------------------------------------------------//

TEST_F(Teuchos_Test, pl_write)
{
    ParameterList p("xs");
    ParameterList m0("mat 0"), m5("mat 10");

    double tot0[] = {1.0, 2.0};
    double tot5[] = {10.0, 11.0};

    double sct0[][2] = {{0.7, 0.0}, {0.2, 1.8}};
    double sct5[][2] = {{0.4, 0.1}, {0.5, 1.2}};

    Array_Dbl sig0(tot0, tot0 + 2);
    Array_Dbl sig5(tot5, tot5 + 2);

    TwoDArray sigs0(2, 2, 0.0);
    TwoDArray sigs5(2, 2, 0.0);

    for (int g = 0; g < 2; ++g)
    {
        for (int gp = 0; gp < 2; ++gp)
        {
            sigs0(g, gp) = sct0[g][gp];
            sigs5(g, gp) = sct5[g][gp];
        }
    }

    m0.set("total", sig0);
    m5.set("total", sig5);

    m0.set("scat", sigs0);
    m5.set("scat", sigs5);

    p.set("num groups", 2);
    p.set("mat 0", m0);
    p.set("mat 5", m5);

    Teuchos::writeParameterListToXmlOStream(p, std::cout);
}

//---------------------------------------------------------------------------//

TEST_F(Teuchos_Test, pl_read)
{
    RCP_ParameterList p = Teuchos::rcp(new ParameterList("cross sections"));
    Teuchos::updateParametersFromXmlFileAndBroadcast(
        "xs3GP0.xml", p.ptr(), *comm);

    EXPECT_EQ(4, p->numParams());
    EXPECT_TRUE(p->isParameter("num groups"));
    EXPECT_EQ(3, p->get<int>("num groups"));

    std::set<std::string> xsnames;
    for (ParameterList::ConstIterator itr = p->begin(); itr != p->end(); ++itr)
    {
        if (p->isSublist(itr->first))
        {
            xsnames.insert(itr->first);
        }
    }
    EXPECT_EQ(2, xsnames.size());

    for (std::set<std::string>::const_iterator itr = xsnames.begin();
         itr != xsnames.end(); ++itr)
    {
        const ParameterList &m = p->sublist(*itr);
        EXPECT_TRUE(m.isParameter("sigma_t"));
        EXPECT_TRUE(m.isParameter("sigma_s0"));

        const Array_Dbl &sigma_t = m.get<Array_Dbl>("sigma_t");

        if (*itr == "mat 0")
        {
            EXPECT_EQ(1.0, sigma_t[0]);
            EXPECT_EQ(2.0, sigma_t[1]);
            EXPECT_EQ(3.0, sigma_t[2]);
        }
        else
        {
            EXPECT_EQ(10.0, sigma_t[0]);
            EXPECT_EQ(11.0, sigma_t[1]);
            EXPECT_EQ(12.0, sigma_t[2]);
        }
    }
}

//---------------------------------------------------------------------------//

TEST_F(Teuchos_Test, get)
{
    ParameterList p;
    EXPECT_FALSE(p.isParameter("A"));
    p.get("A", 1);
    EXPECT_TRUE(p.isParameter("A"));
    EXPECT_EQ(1, p.get<int>("A"));
}

//---------------------------------------------------------------------------//

TEST_F(Teuchos_Test, sublists)
{
    RCP_ParameterList p = Teuchos::rcp(new ParameterList("stuff"));
    p->set("A", 1);
    p->set("B", 2.5);

    EXPECT_EQ(2, p->numParams());

    {
        RCP_ParameterList sub = Teuchos::sublist(p, "MySubList");
        sub->set("subA", 10);
        sub->set("subB", 4.1);
    }

    EXPECT_EQ(3, p->numParams());
    EXPECT_TRUE(p->isSublist("MySubList"));

    const ParameterList &s = p->sublist("MySubList");
    EXPECT_EQ(2, s.numParams());

    {
        ParameterList &s2 = p->sublist("MySubList2");
        EXPECT_EQ(4, p->numParams());
        s2.set("C", 12);
    }

    RCP_ParameterList sub2 = Teuchos::sublist(p, "MySubList2");
    EXPECT_EQ(4, p->numParams());
    EXPECT_EQ(1, sub2->numParams());
}

//---------------------------------------------------------------------------//

TEST_F(Teuchos_Test, vector)
{
    Array_Dbl x(4, 0.0);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    x[3] = 4.0;

    Vector v(Teuchos::Copy, x.getRawPtr(), 4);
    v[3] = 5.0;

    EXPECT_EQ(x[0], v[0]);
    EXPECT_EQ(x[1], v[1]);
    EXPECT_EQ(x[2], v[2]);
    EXPECT_EQ(5.0, v[3]);
    EXPECT_EQ(4.0, x[3]);

    EXPECT_DOUBLE_EQ(11.0, v.normOne());
}

//---------------------------------------------------------------------------//

TEST_F(Teuchos_Test, validation)
{
    RCP_ParameterList vA = Teuchos::rcp(new ParameterList("valid"));
    RCP_ParameterList A  = Teuchos::rcp(new ParameterList("A"));

    // valid parameters
    vA->set("a", 0);
    vA->set("b", 0);
    vA->set("c", 0);
    vA->set("d", 0);

    // actual list
    A->get("a", 10);
    A->get("b", 11);
    A->get("c", 12);
    A->get("d", 13);

    A->validateParameters(*vA);
}

//---------------------------------------------------------------------------//

TEST_F(Teuchos_Test, matrix)
{
    TwoDArray x(3, 3, 0.0);
    x(0, 0) = 1.0;
    x(0, 1) = 2.0;
    x(1, 0) = 10.0;
    x(1, 1) = 20.0;
    x(1, 2) = 30.0;
    x(2, 0) = 100.0;
    x(2, 1) = 200.0;
    x(2, 2) = 300.0;

    Matrix m(3, 3);
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                m(i, j) = x(i, j);
            }
        }
    }

    m(2, 1) = 205.0;

    EXPECT_EQ(x(0, 0), m(0, 0));
    EXPECT_EQ(x(0, 1), m(0, 1));
    EXPECT_EQ(x(0, 2), m(0, 2));
    EXPECT_EQ(x(1, 0), m(1, 0));
    EXPECT_EQ(x(1, 1), m(1, 1));
    EXPECT_EQ(x(1, 2), m(1, 2));
    EXPECT_EQ(x(2, 0), m(2, 0));
    EXPECT_EQ(205.0, m(2, 1));
    EXPECT_EQ(200.0, x(2, 1));
    EXPECT_EQ(x(2, 2), m(2, 2));

    EXPECT_EQ(1.0,   *(m.values()));
    EXPECT_EQ(10.0,  *(m.values()+1));
    EXPECT_EQ(100.0, *(m.values()+2));
    EXPECT_EQ(2.0,   *(m.values()+3));
    EXPECT_EQ(20.0,  *(m.values()+4));
    EXPECT_EQ(205.0, *(m.values()+5));
    EXPECT_EQ(0.0,   *(m.values()+6));
    EXPECT_EQ(30.0,  *(m.values()+7));
    EXPECT_EQ(300.0, *(m.values()+8));
}

//---------------------------------------------------------------------------//
//                 end of tstTeuchos.cc
//---------------------------------------------------------------------------//
