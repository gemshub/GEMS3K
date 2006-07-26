//-------------------------------------------------------------------
// $Id: verror.h 774 2006-07-26 08:45:45Z gems $
//
// Error handling classes TError and TFatalError
//
// Copyright (C) 1996-2001 A.Rysin, S.Dmytriyeva
// Uses  gstring class (C) A.Rysin 1999
//
// This file is part of the GEM-Vizor library which uses the
// Qt v.2.x GUI Toolkit (Troll Tech AS, http://www.trolltech.com)
// according to the Qt Duo Commercial license
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI/ for more information
// E-mail gems2.support@psi.ch
//-------------------------------------------------------------------
#ifndef _verror_h_
#define _verror_h_

#include "gstring.h"

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


#endif    // _verror_h

