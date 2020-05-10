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

    TError(const std::string& titl, const std::string& msg):
            mess(msg),
            title(titl)
    {}

    TError(const TError& other ):
            mess(other.mess),
            title(other.title)
    {}

    virtual ~TError() = default;


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


inline
void Error (const std::string& title, const std::string& message)
{
    throw TError(title, message);
}

inline
void ErrorIf (bool error, const std::string& title, const std::string& message)
{
    if(error)
        throw TError(title, message);
}

//#ifndef IPMGEMPLUGIN

//#include "gstring.h"

//inline
//void Error (const gstring& title, const gstring& message)
//{
//    throw TError(title.c_str(), message.c_str());
//}

//inline
//void ErrorIf (bool error, const gstring& title, const gstring& message)
//{
//    if(error)
//        throw TError(title.c_str(), message.c_str());
//}

//#endif

inline
void Error (const char* title, const char* message)
{
    throw TError(title, message);
}

inline
void ErrorIf (bool error, const char* title, const char* message)
{
    if(error)
        throw TError(title, message);
}


#endif  // VERROR_H

