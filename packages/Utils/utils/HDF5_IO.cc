//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/HDF5_IO.cc
 * \author Thomas M. Evans
 * \date   Thu Nov 07 21:53:44 2013
 * \brief  HDF5_IO member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "HDF5_IO.hh"

#include <vector>
#include <sstream>

#include "harness/Warnings.hh"
#include "comm/global.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// HDF5_IO BASE CLASS DEFINITIONS
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
HDF5_IO::HDF5_IO()
    : b_file(0)
    , b_node(profugus::node())
    , b_nodes(profugus::nodes())
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Virtual destructor.
 */
HDF5_IO::~HDF5_IO()
{
    try
    {
        HDF5_IO::close();
    }
    catch (const profugus::assertion& ass)
    {
        ADD_WARNING("While closing HDF5 file at " << b_filename
                << ", caught error during destructor: "
                << ass.what());
    }
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Close the HDF5 file.
 *
 * This is a null-op if the file is already closed.
 */
void HDF5_IO::close()
{
    if (b_file == 0)
        return;

    // Check file is at end of group
    Validate(b_loc_stack.size() <= 1, "HDF5 group was not ended. "
            "Last group name was '" << b_loc_stack.back().group_name << "'.");

    // close file on node 0
    herr_t status = H5Fclose(b_file);
    Insist (status >= 0, "HDF5 failed to close file.");

    b_mode = END_FILE_MODE;
    b_loc_stack.clear();

    b_file = 0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Query to see if keyed value is in the current group.
 */
bool HDF5_IO::query(const std_string &name) const
{
    return H5LTfind_dataset(current_loc(), name.c_str());
}

//---------------------------------------------------------------------------//
// HIERARCHY INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Create a new group inside the current group
 */
void HDF5_IO::begin_group(const std_string &dirname)
{
    Require(dirname.find('/') == std_string::npos); // disallow slashes

    Remember(size_t ng_orig = b_loc_stack.size());

    // group
    hid_t new_group = 0;

    // Create the group if writing
    if (!is_readonly())
    {
        new_group = H5Gcreate(
            current_loc(),   // file/group identifier
            dirname.c_str(), // relative name of the link to the new group
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT // property list identifiers
            );
    }
    else
    {
        new_group = H5Gopen(current_loc(), dirname.c_str(), H5P_DEFAULT);
    }

    Validate(new_group >= 0, "HDF5 failed to create/open group '"
             << b_loc_stack.back().group_name << "'");

    // Add it to our stack
    HDF5_Loc loc = {new_group, dirname};
    b_loc_stack.push_back(loc);

    Ensure(b_loc_stack.size() == ng_orig + 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Close the current group, returning to the parent group
 */
void HDF5_IO::end_group()
{
    Insist(b_loc_stack.size() > 1, "Can't close root group.");
    Remember(size_t ng_orig = b_loc_stack.size());

    herr_t status = H5Gclose(b_loc_stack.back().handle);
    if (status < 0)
        ADD_WARNING("HDF5 failed to close group '"
                << b_loc_stack.back().group_name << "'");

    b_loc_stack.pop_back();

    Ensure(b_loc_stack.size() == ng_orig - 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the absolute path to the current group.
 *
 * The result will contain the correct leading slash and a trailing slash.
 */
HDF5_IO::std_string HDF5_IO::abspath() const
{
    std::ostringstream result;
    for (std::vector<HDF5_Loc>::const_iterator it = b_loc_stack.begin(),
            end_it = b_loc_stack.end();
            it != end_it;
            ++it)
    {
        result << it->group_name << "/";
    }

    return result.str();
}

//---------------------------------------------------------------------------//
// PARALLEL_HDF5_WRITER::DECOMP DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
HDF5_IO::Decomp::Decomp(int        n,
                        Data_Order o)
    : ndims(n)
    , global(n, 0)
    , local(n, 0)
    , offset(n, 0)
    , order(COLUMN_MAJOR)
{
}

//---------------------------------------------------------------------------//
// HDF5_IO_NODE CLASS DEFINITIONS
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
HDF5_IO_Node::HDF5_IO_Node()
    : Base()
    , b_master(0)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Virtual destructor.
 */
HDF5_IO_Node::~HDF5_IO_Node()
{
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Close the file.
 */
void HDF5_IO_Node::close()
{
    // only close the file on the node that it was opened on
    if (b_node == b_master)
    {
        Base::close();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Query to see if keyed value is in the current group.
 */
bool HDF5_IO_Node::query(const std_string &name) const
{
    if (b_node != b_master)
        return false;
    return Base::query(name);
}

//---------------------------------------------------------------------------//
// HIERARCHY INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Create a new group inside the current group.
 *
 * Only runs on the designated host node.
 */
void HDF5_IO_Node::begin_group(const std_string &dirname)
{
    if (b_node != b_master)
        return;
    Base::begin_group(dirname);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Close the current group, returning to the parent group
 *
 * Only runs on the designated host node.
 */
void HDF5_IO_Node::end_group()
{
    if (b_node != b_master)
        return;
    Base::end_group();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the absolute path to the current group.
 *
 * The result will contain the correct leading slash and a trailing slash.
 *
 * Only runs on the designated host node.
 */
HDF5_IO_Node::std_string HDF5_IO_Node::abspath() const
{
    if (b_node != b_master)
        return "";
    return Base::abspath();
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of HDF5_IO.cc
//---------------------------------------------------------------------------//
