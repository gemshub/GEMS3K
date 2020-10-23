//-------------------------------------------------------------------
// $Id$
//
/// \file io_template.cpp
/// Implementation of service functions for writing/reading arrays in files
//
// Copyright (c) 2020 S.Dmytriyeva
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------

#include "io_template.h"

namespace  io_formats {

long int TRWArrays::findFld(  const std::string& name  ) const
{
    for( long int ii=0; ii < num_flds; ii++ )
        if( flds[ii].name == name )
            return ii;
    return -1;
}

std::string TRWArrays::testRead() const
{
    std::string ret = "";
    for( long int ii=0; ii < num_flds; ii++ )
        if( flds[ii].alws==1 && flds[ii].readed != 1 )
        {
            if( !ret.empty() )
                ret += ", ";
            ret += flds[ii].name;
        }
    return ret;
}

}  // io_formats
