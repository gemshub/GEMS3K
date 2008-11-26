// $Id: gstring.h 968 2007-12-13 13:23:32Z gems $
/***************************************************************************
	gstring class
	Version 1.02
        04/08/2000
 
	represents some (vital) part
	of stl::string functionality
	without copy-on-write technology
	(there's a lot of improvements could be done here :-)

// This file is part of the GEM-Vizor library and GEMIPM2K code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
//
	Copyright (c) 2000
	Author: Andy Rysin
	E-mail: arysin@yahoo.com
****************************************************************************/

#ifndef _gstring_h_
#define _gstring_h_

#include <ctype.h>
#include <string.h>
#include "array.h"

#ifdef IPMGEMPLUGIN
   static const size_t npos = static_cast<size_t>(-1);
//   static  const size_t npos=32767;   /wp sergey 2004 from below assignment

#endif

class gstring
{
//    internal class for keeping string values
struct str_val:
                public TOArray<char>
    {
        friend class gstring;
        str_val(size_t sz=30, size_t granul=30)
                :TOArray<char>(sz, granul)
        { }

        ~str_val()
        {}


        bool equals(const str_val& s) const;
        str_val* clone(size_t add_size=0);
        const char* c_str()
        {
            if( elem(GetCount()) )
            {
                Add('\0');
                count--;
            }
            return p;
        }
    };

    str_val* ps;

public:
#ifndef IPMGEMPLUGIN
    static const size_t npos = static_cast<size_t>(-1);
#endif
    // various constructors and destructor
    gstring()
    {
        ps = new str_val();
    }
    gstring(const gstring& s)
    {
        ps = s.ps->clone();
    }
    gstring(const gstring& s, size_t pos, size_t len=npos);
    gstring(const char* s, size_t pos=0, size_t len=npos);
    gstring(int num, char ch);

    ~gstring()
    {
        delete ps;
    }

// Selectors

    // returns length of the string
    size_t length() const
    {
        return ps->GetCount();
    }
    // returns true if string is empty - just for convinience
    size_t empty() const
    {
        return length() == 0;
    }

    // find something from position 'pos' or beginning if omited
    size_t find(const char* s, size_t pos = 0) const;
    size_t find(const gstring& s, size_t pos = 0) const;
    size_t find(const char ch, size_t pos = 0) const;
    size_t find_first_of(const char* s) const;

    // reverse find
    size_t rfind(const char* s) const;

    // returns substring
    gstring substr(size_t pos, size_t len = npos) const;



// Modificators
    
    // replaces the last occurance of the old_part with new_part
    gstring replace(const char* old_part,
                      const char* new_part);
    // erases from position 'p1' and size 'sz'
    void erase(size_t p1, size_t p2);

    // various gstring assigning
    const gstring& operator=(const gstring& s)
    {
        delete ps;
        ps = s.ps->clone();
        return *this;
    }
    const gstring& operator+=(const gstring& s);

// comparisons
    bool operator==(const gstring& s) const
    {
        return ps->equals(*s.ps);
    }
    
    bool operator!=(const gstring& s) const
    {
        return !(*this == s);
    }
    
    bool equals(const char* s) const;

    // various C-string assigning
    const gstring& operator=(const char* s);
    const gstring& operator+=(const char* s);
    const gstring& operator+=(const char ch)
    {
        ps->Add(ch);
        return *this;
    }

    // self-explanatory
    char& operator[](size_t ii)
    {
        return ps->elem(ii);
    }

    // returns c-style string pointer
    const char* c_str() const
    {
        return ps->c_str();
    }

    // strip spaces from both ends
    void strip();
};

inline
bool operator==(const gstring& str, const char* sp)
{
    return str.equals(sp);
}

inline
bool operator!=(const gstring& str, const char* sp)
{
    return !(str == sp);
}

inline
bool operator==(const char* sp, const gstring& str)
{
    return (str == sp);
}

inline
bool operator!=(const char* sp, const gstring& str)
{
    return !(str == sp);
}

inline
const gstring operator+(const gstring& str1, const gstring& str2)
{
    gstring res(str1);
    return res += str2;
}

// added for convenience because of frequent use
typedef TArrayF<gstring> TCStringArray;

#endif
//_gstring_h_
