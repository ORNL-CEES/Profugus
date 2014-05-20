//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cxx11/tstPointers.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 26 13:03:29 2014
 * \brief  CXX-11 test of pointers (including smart pointers).
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <cstring>
#include <typeinfo>
#include <memory>

//---------------------------------------------------------------------------//
// NullPtr tests
//---------------------------------------------------------------------------//

template<class T>
T doit(const T *data)
{
    if (data)
    {
        return *data;
    }

    return T();
}

//---------------------------------------------------------------------------//

TEST(NullPtr, definitions)
{
    int *x = nullptr;
    EXPECT_FALSE(x);

    int y = 12;
    x     = &y;
    EXPECT_TRUE(x);
}

//---------------------------------------------------------------------------//

TEST(NullPtr, arguments)
{
    double x = 12.5;
    int    y = 2;

    double *p = &x;
    int    *q = &y;

    EXPECT_EQ(12.5, doit(p));
    EXPECT_EQ(2,    doit(q));

    p = nullptr;
    q = nullptr;

    EXPECT_EQ(0.0, doit(p));
    EXPECT_EQ(0,   doit(q));
}

//---------------------------------------------------------------------------//
// SHARED_PTR
//---------------------------------------------------------------------------//

int nfoos = 0;
int nbars = 0;

class Foo
{
  protected:
    int v;

  public:
    Foo()
        : v(0)
    {
        nfoos++;
    }

    explicit Foo(int i)
        : v(i)
    {
        nfoos++;
    }

    Foo(const Foo &f)
        : v(f.v)
    {
        nfoos++;
    }

    virtual ~Foo()
    {
        nfoos--;
    }

    virtual int vf() { return v; }

    int f() { return v+1; }

    int cf() const { return v+1; }
};

//---------------------------------------------------------------------------//

class Bar : public Foo
{
  private:
    Bar(const Bar &);

  public:
    explicit Bar(int i)
        : Foo(i)
    {
        nbars++;
    }

    virtual ~Bar()
    {
        nbars--;
    }

    virtual int vf() { return Foo::f() + 1; }

    int f() { return Foo::f() + 2; }

    void set(int j) { v = j; }
};

//---------------------------------------------------------------------------//

class Test_shared_ptr : public testing::Test
{
  protected:
    typedef std::shared_ptr<Foo> SP_Foo;
    typedef std::shared_ptr<Bar> SP_Bar;

  protected:

    void SetUp()
    {
        expect(0, 0);
    }

    void TearDown()
    {
        expect(0, 0);
    }

    void expect(int nf, int nb)
    {
        EXPECT_EQ(nf, nfoos);
        EXPECT_EQ(nb, nbars);
    }
};

//---------------------------------------------------------------------------//

TEST_F(Test_shared_ptr, construction)
{
    // make a Foo and Bar
    {
        SP_Foo spfoo(new Foo(1));
        SP_Bar spbar(new Bar(2));

        // there should be 2 Foos, 1 Bars
        expect(2, 1);

        EXPECT_EQ(1, spfoo.use_count());
        EXPECT_EQ(1, spbar.use_count());

        EXPECT_TRUE(spfoo.unique());
        EXPECT_TRUE(spbar.unique());
    }
}

//---------------------------------------------------------------------------//

TEST_F(Test_shared_ptr, bool)
{
    {
        SP_Foo spfoo(new Foo(1));
        SP_Bar spbar(new Bar(2));

        if (!spfoo)
        {
            FAIL() << "foo ptr not assigned.";
        }

        if (!spbar)
        {
            FAIL() << "foo ptr not assigned.";
        }
    }
}

//---------------------------------------------------------------------------//

TEST_F(Test_shared_ptr, types)
{
    {
        SP_Foo spfoo = std::make_shared<Foo>(10);
        SP_Bar spbar = std::make_shared<Bar>(20);
        SP_Foo spfb  = std::make_shared<Bar>(30);

        expect(3, 2);

        EXPECT_EQ(11, spfoo->f());
        EXPECT_EQ(10, spfoo->vf());
        EXPECT_EQ(11, spfoo->cf());

        EXPECT_EQ(23, spbar->f());
        EXPECT_EQ(22, spbar->vf());
        EXPECT_EQ(21, spbar->cf());

        EXPECT_EQ(31, spfb->f());
        EXPECT_EQ(32, spfb->vf());
        EXPECT_EQ(31, spfb->cf());
    }
}

//---------------------------------------------------------------------------//

TEST_F(Test_shared_ptr, assignment)
{
    {
        // make some objects
        SP_Foo spfoo;
        SP_Bar spbar;
        SP_Foo spfoo2;
        {
            spbar = std::make_shared<Bar>(50);
            expect(1, 1);

            EXPECT_EQ(53, spbar->f());
            EXPECT_EQ(52, spbar->vf());

            // now assign to base class SP
            spfoo = spbar;
            expect(1, 1);

            EXPECT_FALSE(spfoo.unique());
            EXPECT_FALSE(spbar.unique());

            EXPECT_EQ(51, spfoo->f());
            EXPECT_EQ(52, spfoo->vf());

            EXPECT_EQ(typeid(Foo *), typeid(spfoo.get()));
            EXPECT_EQ(typeid(Bar), typeid(*spfoo.get()));
            EXPECT_EQ(typeid(Bar *), typeid(spbar.get()));

            // now do copy construction
            SP_Foo rspfoo(spbar);
            expect(1, 1);

            EXPECT_FALSE(rspfoo.unique());

            EXPECT_EQ(51, rspfoo->f());
            EXPECT_EQ(52, rspfoo->vf());

            EXPECT_EQ(typeid(Foo *), typeid(rspfoo.get()));
            EXPECT_EQ(typeid(Bar), typeid(*rspfoo.get()));

            // now check assignment with X *
            rspfoo = std::make_shared<Bar>(12);
            expect(2, 2);

            EXPECT_TRUE(rspfoo.unique());
            EXPECT_FALSE(spfoo.unique());
            EXPECT_FALSE(spbar.unique());

            EXPECT_EQ(13, rspfoo->f());
            EXPECT_EQ(14, rspfoo->vf());

            EXPECT_EQ(typeid(Foo *), typeid(rspfoo.get()));
            EXPECT_EQ(typeid(Bar), typeid(*rspfoo.get()));

            // assign SPfoo2 to a bar
            spfoo2 = std::make_shared<Bar>(20);
            expect(3, 3);

            EXPECT_TRUE(rspfoo.unique());
            EXPECT_TRUE(spfoo2.unique());
            EXPECT_FALSE(spfoo.unique());
            EXPECT_FALSE(spbar.unique());
        }

        // still have 2 objects
        expect(2, 2);

        spfoo.reset(new Foo(10));
        expect(3, 2);
        EXPECT_TRUE(spfoo2.unique());
        EXPECT_TRUE(spfoo.unique());
        EXPECT_TRUE(spbar.unique());

        EXPECT_EQ(11, spfoo->f());
        EXPECT_EQ(10, spfoo->vf());

        spbar.reset();
        expect(2, 1);

        EXPECT_TRUE(!static_cast<bool>(spbar));

        std::swap(spfoo, spfoo2);
        expect(2, 1);

        EXPECT_EQ(11, spfoo2->f());
        EXPECT_EQ(10, spfoo2->vf());

        EXPECT_EQ(21, spfoo->f());
        EXPECT_EQ(22, spfoo->vf());
    }
}

//---------------------------------------------------------------------------//

TEST_F(Test_shared_ptr, equivalence)
{
    SP_Foo a(std::make_shared<Foo>());
    SP_Foo b = a;

    EXPECT_EQ(a, b);

    SP_Bar c(std::make_shared<Bar>(12));
    SP_Bar d(c);

    d->set(11);

    EXPECT_EQ(c, d);

    SP_Foo e = d;

    EXPECT_EQ(c, e);

    EXPECT_EQ(12, e->f());
    EXPECT_EQ(13, e->vf());
    EXPECT_EQ(12, e->cf());

    EXPECT_EQ(14, d->f());
    EXPECT_EQ(13, d->vf());
    EXPECT_EQ(12, d->cf());
}

//---------------------------------------------------------------------------//
//                 end of tstPointers.cc
//---------------------------------------------------------------------------//
