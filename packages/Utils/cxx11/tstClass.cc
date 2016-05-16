//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/cxx11/tstClass.cc
 * \author Thomas M. Evans
 * \date   Wed May 14 09:36:01 2014
 * \brief  C++-11 Class-features.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include <iostream>

//---------------------------------------------------------------------------//
// Simple Inheritance
//---------------------------------------------------------------------------//

class Foo
{
  public:
    Foo() {/*...*/}
    virtual ~Foo() = 0;
    virtual void doit() { std::cout << "Foo knows nothing." << std::endl; }
    virtual void mine() = 0;
};

Foo::~Foo()
{
}

//---------------------------------------------------------------------------//

class Bar : public Foo
{
  private:
    int d_y;

  public:
    Bar(int y) : d_y(y) {/*...*/}

    virtual void doit() override
    {
        std::cout << "Here is Bar's y: " << d_y << std::endl;
    }

    virtual void mine() final
    {
        std::cout << "Bar owns the world!" << std::endl;
    }

  protected:

    int y() const { return d_y; }
};

//---------------------------------------------------------------------------//

class Baz : public Bar
{
  public:
    Baz() : Bar(11) {/*...*/}

    virtual void doit() override
    {
        std::cout << "Baz says y = " << Bar::y() << std::endl;
    }

    // Cannot override mine; it is final.
};

//---------------------------------------------------------------------------//

class Bat final : public Foo
{
  private:
    int d_y;

  public:
    Bat(int y) : d_y(y) {/*...*/}

    void doit()
    {
        std::cout << "Here is Bat's y: " << d_y << std::endl;
    }

    void mine()
    {
        std::cout << "Bat doesn't mind!" << std::endl;
    }
};

//---------------------------------------------------------------------------//
/*
  Cannot override Bat it is final.

  class Biz : public Bat
  {
  };

 */

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Inheritance, override)
{
    Bar b(10);
    b.mine();
    b.doit();

    Foo &f = b;
    f.mine();
    f.doit();
}

//---------------------------------------------------------------------------//

TEST(Inheritance, final)
{
    Baz b;
    b.mine();
    b.doit();

    Foo &f = b;
    f.mine();
    f.doit();

    Bar &r = b;
    r.mine();
    r.doit();
}

//---------------------------------------------------------------------------//

TEST(Inheritance, class_final)
{
    Bat b(100);
    b.doit();
    b.mine();
}

//---------------------------------------------------------------------------//
//                 end of tstClass.cc
//---------------------------------------------------------------------------//
