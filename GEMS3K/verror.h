//-------------------------------------------------------------------
// $Id$
//
// Error handling: classes TError and TFatalError
//
// Copyright (C) 1996-2001 A.Rysin, S.Dmytriyeva
// Uses  gstring class (C) A.Rysin 1999
//
// This file is part of the GEM-Selektor GUI library and GEMS3K
// code package
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail gems2.support@psi.ch
//-------------------------------------------------------------------
#ifndef _verror_h_
#define _verror_h_

#ifdef IPMGEMPLUGIN

#include <string>

using namespace std;
typedef string gstring;
static const size_t npos = string::npos;
//   static const size_t npos = static_cast<size_t>(-1);
//   static  const size_t npos=32767;   /wp sergey 2004 from below assignment

void strip(string& str);

#else

#include "gstring.h"

#endif

struct TError
{
    gstring mess;
    gstring title;
    TError()
    {}

    TError(const gstring& titl, const gstring& msg):
            mess(msg),
            title(titl)
    {}

    virtual ~TError()
    {}

};


struct TFatalError:
            public TError
{
    TFatalError()
    {}

    TFatalError(TError& err):
            TError(err)
    {}

    TFatalError(const gstring& titl, const gstring& msg):
            TError( titl, msg )
    {}

};


inline
void Error(const gstring& title, const gstring& message)
{
    throw TError(title, message);
}

inline
void ErrorIf(bool error, const gstring& title, const gstring& message)
{
    if(error)
        throw TError(title, message);
}


#endif
// _verror_h

