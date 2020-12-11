//-------------------------------------------------------------------
// $Id$
/// \file verror.h
/// Declarations of classes TError and TFatalError for error handling.
//
// Copyright (C) 1996-2012 A.Rysin, S.Dmytriyeva
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
#ifndef VERROR_H
#define VERROR_H

#include <string>

struct TError
{
    std::string mess;
    std::string title;
    TError()
    {}

    TError( const std::string& titl, const std::string& msg):
        mess(msg),
        title(titl)
    {}

//    TError( const TError& other ):
//        mess(other.mess),
//        title(other.title)
//    {}

    virtual ~TError();
};


struct TFatalError:
        public TError
{
    TFatalError()
    {}

    TFatalError(const TError& err):
        TError(err)
    {}

    TFatalError(const std::string& titl, const std::string& msg):
        TError( titl, msg )
    {}

};


[[ noreturn ]]  void Error ( const std::string& title, const std::string& message );
void ErrorIf ( bool error, const std::string& title, const std::string& message );


#endif  // VERROR_H

