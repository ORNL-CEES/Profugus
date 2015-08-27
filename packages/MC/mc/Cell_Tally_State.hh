//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Cell_Tally_State.hh
 * \author Thomas M. Evans
 * \date   Thu Aug 27 11:24:25 2015
 * \brief  Cell_Tally_State class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Cell_Tally_State_hh
#define MC_mc_Cell_Tally_State_hh

#include <unordered_map>

#include "Utils/harness/DBC.hh"
#include "Utils/utils/Member_Manager.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Cell_Tally_State
 * \brief Particle-private statue for cell tallies.
 */
//===========================================================================//

class Cell_Tally_State
{
  public:
    typedef std::unordered_map<int, double> History_Tally;

  private:
    // >>> DATA

    // Tally for a history.
    History_Tally d_hist;

  public:
    Cell_Tally_State() { /* * */ }

    //! Get the history of the tally.
    const History_Tally& state() const { return d_hist; }

    //@{
    //! Add contribution to a cell.
    double& operator()(int cell) { return d_hist[cell]; }
    //@}
};

//---------------------------------------------------------------------------//
// METADATA FOR CELL_TALLY_STATE
//---------------------------------------------------------------------------//

namespace metaclass
{

//===========================================================================//
/*!
 * \class Cell_Tally_State
 * \brief Particle-private statue for cell tallies.
 */
//===========================================================================//

class Member_Manager_Cell_Tally_State : public Member_Manager
{
  public:
    typedef Cell_Tally_State State_t;

  public:

    //! Call the constructor on already allocated memory
    inline void construct(void* data)
    {
        new (data) State_t();
    }

    //! Use the copy constructor to initialize already allocated memory
    inline  void copy_construct(void* data, const void* rhs_data)
    {
        const State_t* rhs = reinterpret_cast<const State_t*>(rhs_data);
        // Calling this constructor then tells the STL vector to allocate
        // space for all the entries, and to copy those.
        new (data) State_t(*rhs);
    }

    inline const char* unpack_construct(void* data, const char* buffer)
    {
        NOT_IMPLEMENTED("Cannot unpack Cell_Tally_State.");
    }

    //! Use operator= for assignment
    inline void assign(void* lhs_data, const void* rhs_data)
    {
        State_t* lhs       = reinterpret_cast<State_t*>(lhs_data);
        const State_t* rhs = reinterpret_cast<const State_t*>(rhs_data);

        *lhs = *rhs;
    }

    //! Call the destructor
    inline void destroy(void* data)
    {
        State_t* state = reinterpret_cast<State_t*>(data);
        state->~State_t();
    }

    //! Pack all items from the vector into the buffer
    inline char* pack(const void* data, char* buffer)
    {
        NOT_IMPLEMENTED("Cannot pack Cell_Tally_State.");
    }

    inline size_type packed_size(const void* data) const
    {
        NOT_IMPLEMENTED("Cannot do size of Cell_Tally_State.");
    }
};

} // end namespace metaclass

//---------------------------------------------------------------------------//

} // end namespace profugus

#endif // MC_mc_Cell_Tally_State_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Cell_Tally_State.hh
//---------------------------------------------------------------------------//
