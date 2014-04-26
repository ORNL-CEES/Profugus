//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Bank.hh
 * \author Seth R Johnson and Thomas M. Evans
 * \date   Friday April 25 16:46:18 2014
 * \brief  Bank class definition
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Bank_hh
#define mc_Bank_hh

#include <vector>
#include <memory>
#include "harness/DBC.hh"
#include "Particle.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Bank
 * \brief Base class for storage for particles
 *
 * The bank is essentially a stack of particles. To reduce storage, we use a
 * delayed-copy approach.
 */
/*!
 * \example mc/test/tstBank.cc
 *
 * Test of Bank.
 */

//===========================================================================//

class Bank
{
  public:
    //@{
    //! Useful typedefs.
    typedef Particle                    Particle_t;
    typedef std::shared_ptr<Particle_t> SP_Particle;
    typedef SP_Particle                 value_type;
    typedef size_t                      size_type;
    typedef SP_Particle&                reference;
    typedef const SP_Particle&          const_reference;
    //@}

  private:
    // Container type for particles
    typedef std::vector<SP_Particle> Stack_Particle;

    // Container type for number of copies per particle
    typedef std::vector<size_type> Stack_Count;

  private:
    // >>> DATA

    // Stored particles
    Stack_Particle d_particles;

    // Number of copies per particle
    Stack_Count d_count;

    // Number of total particles (sum of d_count)
    size_type d_total;

  public:
    Bank() : d_particles(), d_count(), d_total(0) { /* * */ }

    //! \name std::stack-like operatinos
    //@{

    //! Is the bank empty?
    bool empty() const
    {
        Check(d_particles.size() == d_count.size());
        Check(d_particles.empty() ? d_total == 0 : true);
        return d_total == 0;
    }

    //! Return the number of particles left to emit
    size_type size() const { return d_total; }
    size_type num_particles() const { return size(); }

    //! View the particle on the top of the stack
    const SP_Particle& top() const
    {
        Require(!empty());
        Require(!d_particles.empty());
        return d_particles.back();
    }
    const SP_Particle& back() const { return top(); }

    //! Just pushing a particle
    void basic_push(const SP_Particle& p, size_type count)
    {
        Require(p);
        Require(count > 0);

        // Add a copy of the particle to the stack
        d_particles.push_back(std::make_shared<Particle_t>(*p));
        d_count.push_back(count);

        d_total += count;
    }

    //! Push a particle (default to basic push)
    void push(const SP_Particle& p, size_type count = 1)
    {
        basic_push(p, count);
    }

    void push_back(const SP_Particle& p) { push(p, 1); }

    //! Emit the topmost particle from the stack
    SP_Particle pop() { return basic_pop(); }
    SP_Particle pop_back() { return pop(); }

    //@}

  public:
    //! \name Implementation-specific operatinos
    //@{

    //! Return the number of unique particles being stored
    size_type num_unique() const { return d_particles.size(); }

    //! Return the number of copies of the next particle
    size_type next_count() const
    {
        return (empty() ? 0 : d_count.back());
    }

    //@}

  private:
    // >>> IMPLEMENTATION

    //! Just emitting the topmost particle from the stack
    inline SP_Particle basic_pop();
};

} // end namespace profugus

//---------------------------------------------------------------------------//

#include "Bank.i.hh"

//---------------------------------------------------------------------------//

#endif // mc_Bank_hh

//---------------------------------------------------------------------------//
//              end of mc/Bank.hh
//---------------------------------------------------------------------------//
