//-------------------------------------------------------------------
// $Id: array.h 968 2007-12-13 13:23:32Z gems $
//
// Declaration/implementation of array classes collection
//
// Copyright (C) 1996-2001 A.Rysin
// Uses  gstring class (C) A.Rysin 1999
//
// This file is part of the GEM-Vizor library and GEMIPM2K
// code package
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail gems2.support@psi.ch
//-------------------------------------------------------------------

#ifndef _gems_array_h_
#define _gems_array_h_

#include <ctype.h>

#include <memory>
using namespace std;

#ifdef _INDEX_DEBUG_ON
# define IF_THROW(cond,type) if(cond) throw type();
#else
# define IF_THROW(cond,type)
#endif

// default granularity for memory allocation
const size_t uGRAN = 20;

// Default class for basic memory allocating procedures
template <class T>
class TAllocator
{
    TAllocator (const TAllocator<T>&)
    {}

    const TAllocator<T>& operator= (const TAllocator<T>&)
    {
        return *this;
    }
    T* realloc(size_t sz, size_t old_cnt);

protected:
    size_t gran;	// semi-constant
    size_t count;	// counts elements
    size_t size;	// keeps memory size
    T* p;		// points to array's memory

    size_t resize(size_t ind);

public:
    TAllocator(size_t sz=0, size_t gr=20);
    virtual ~TAllocator()
    {
        delete[] p;
    }
    size_t Remove(size_t ind);
    size_t GetCount() const
    {
        return count;
    }
    void Clear(bool free=true);

#ifdef _INDEX_DEBUG_ON
    struct EIndex
        {};
#endif
};

// Reallocates memory of the array
template <class T>
T*
TAllocator<T>::realloc(size_t sz, size_t old_cnt)
{
    T* pt = new T[sz];
    for( size_t ii=0; ii<old_cnt; ii++ )
        pt[ii] = this->p[ii];

    delete[] this->p;
    this->size = sz;
    return this->p = pt;
}

// Changes the size of the array and reallocates the memory
template <class T>
size_t
TAllocator<T>::resize(size_t ind)
{
    //  IF_THROW( ind<0, EIndex );
    size_t old_cnt = this->count;
    if( ind >= this->count )
        this->count = ind+1;
    if( ind < this->size )
        return this->size;

    size_t sz = ( ind/uGRAN + 1) * uGRAN;
    realloc(sz, old_cnt);

    return this->size;
}

// Constructor. Takes the initial size of the array and granularity value
template <class T>
TAllocator<T>::TAllocator(size_t sz, size_t gr):
        gran(gr),
        count(0),
        size( ( sz <= 0 ) ? uGRAN : ((sz - 1)/uGRAN + 1) * uGRAN ),
        p( new T[size] )
{
    //  IF_THROW( !p, EAlloc );
}

// removes element of the array but
// does not change the size of it's memory
template <class T>
size_t
TAllocator<T>::Remove(size_t ind)
{
    IF_THROW( (ind >= count), EIndex );
    for( size_t ii=ind; ii<count-1; ii++ )
        this->p[ii] = this->p[ii+1];
    this->count--;
    return ind;
}

// clears the array and resets memory allocation
template <class T>
void
TAllocator<T>::Clear(bool freeMem)
{
    delete[] this->p;
    this->count = 0;
    if( freeMem )
        this->size = uGRAN;
    //  else 'size' leaves the same

    p = new T[this->size];
    //  IF_THROW( !p, EAlloc );
}


// Class TOArray store objects as they are and destroys them on cleanup
//
template <class T>
class TOArray:
            public TAllocator<T>
{
public:
    TOArray(size_t sz=0, size_t gr=20):
            TAllocator<T>(sz, gr)
    {}

    const TOArray<T>& operator =(const TOArray<T>& a); // should be not allowed
    TOArray(const TOArray& a);
    virtual ~TOArray()
    {}

// selectors
    T& elem(size_t ind) const;
    T& operator[](size_t ind) const
    {
        return this->elem(ind);
    }
    T* _ptr()
    {
        return this->p;
    }

// modificators
    size_t Add(const T&);
    size_t AddAt(const T&, size_t ind);
    void SetToArray(T* p, size_t cnt);
};

// TOArray copy constructor
template <class T>
TOArray<T>::TOArray(const TOArray& a):
        TAllocator<T>(a.size, a.gran)
{
    this->count = a.count;
    for( size_t ii=0; ii<this->count; ii++ )
        this->p[ii] = a.p[ii];
}

// returns element at position 'index'
template <class T>
inline
T& TOArray<T>::elem(size_t index) const
{
    IF_THROW( (index >= this->count), EIndex );
    // maybe try this:
    // resize( ind );

    return this->p[index];
}

// TOArray assignment operator
template <class T>
const TOArray<T>&
TOArray<T>::operator=(const TOArray<T>& a)
{
    if( this->p == a.p )
        return *this;

    delete[] this->p;
    this->gran = a.gran;
    this->size = a.size;
    this->count = a.count;

    this->p = new T[this->size];
    for( size_t ii=0; ii<this->count; ii++ )
        this->p[ii] = a.p[ii];

    return *this;
}

// Adds element to the array
template <class T>
size_t
TOArray<T>::Add(const T& t)
{
    resize(this->count);
    this->p[this->count - 1] = t;
    return this->count - 1;
}

// Adds element to the array at position 'index'
template <class T>
size_t
TOArray<T>::AddAt(const T& t, size_t index)
{
    if( index < this->count )
        this->resize(this->count);
    else
        this->resize(index);

    for( size_t ii=this->count-1; ii>index; ii-- )
        this->p[ii] = this->p[ii-1];
    this->p[index] = t;
    return index;
}

// Resets the array and fills it with values from C-type array
//
template <class T>
void
TOArray<T>::SetToArray(T* ptr_n, size_t cnt)
{
    delete[] this->p;
    this->p = new T[cnt];
    //  IF_THROW( !p, EAlloc );

    this->count = this->size = cnt;
    for( size_t ii=0; ii<cnt; ii++ )
        this->p[ii] = ptr_n[ii];
}

// TArrayF - extends TOArray with Find() and Sort() functions
//    This allows TOArray to store types wich does not have
//    assignment and '<', '>' operators
//
template <class T>
class TArrayF:
            public TOArray< T >
{
public:
    TArrayF(size_t sz=0, size_t gr=20):
            TOArray<T>(sz, gr)
    {}

    virtual ~TArrayF()
    {}

// 
    int Find(const T& object) const;
    void Sort();
};

// searches the array for element which equals to 'object'
//    and returns index of the element in the array
//    pretty unefficient, but we can't be sure array is sorted
//
template <class T>
int
TArrayF<T>::Find(const T& object) const
{
    for( size_t ii=0; ii<this->GetCount(); ii++ )
        if( object == this->elem(ii) )
            return ii;

    return -1;
}

// sorts the elements of the array
template <class T>
void
TArrayF<T>::Sort()
{
    for( int ii=this->GetCount()-1; ii>0; ii-- )
    {
        int imax = ii;
        for( int jj=ii-1; jj>=0; jj-- )
            if( this->elem(imax) > this->elem(ii) )
            {
                T tmp = this->elem(imax);
                this->elem(imax) = this->elem(jj);
                this->elem(jj) = tmp;
                imax = jj;
            }
    }
}

// Added for convenience
typedef TArrayF<int> TCIntArray;


// TIArray stores pointers to object
//    returns and accepts references on pointed objects
//    and cleans the objects pointed by its elements
//   If you don't want objects to be cleaned use TOArray<type*>
//
template <class T>
class TIArray:
            public TAllocator< auto_ptr<T> >
{
    typedef auto_ptr<T> ptrT;

    // forbidden !!!
    const TIArray<T>& operator =(const TIArray<T>& a); // should be not allowed
    TIArray(const TIArray& a);

public:
    TIArray(size_t sz=0, size_t gr=20):
            TAllocator<ptrT>(sz, gr)
    {}

    virtual ~TIArray()
    {}

// selectors
    size_t Add(T* p);
    size_t AddAt(T* p, size_t ind);
    T& elem(size_t ind) const;
    T& operator[](size_t ind) const
    {
        return this->elem(ind);
    }
    int Find(const T&) const;
    //  T& GetLast() { return elem( GetCount()-1 ); }
    //  void SetToArray(T* p, size_t cnt);
    //  const T** _ptr() { return arr._ptr(); }
};


// Adds pointer to the element to the array
template <class T>
size_t
TIArray<T>::Add(T* t)
{
    resize(this->count);
    this->p[this->count - 1].reset(t);
    return this->count - 1;
}

// Adds pointer to the element to the array at position 'index'
template <class T>
size_t
TIArray<T>::AddAt(T* t, size_t ind)
{
    if( ind < this->count )
        this->resize(this->count);
    else
        this->resize(ind);
    for( size_t ii=this->count-1; ii>ind; ii-- )
        this->p[ii] = this->p[ii-1];
    this->p[ind].reset(t);
    return ind;
}

// returns element pointer at index points to
template <class T>
inline
T&
TIArray<T>::elem(size_t index) const
{
    IF_THROW( (index >= this->count), EAlloc );
    //  resize( ind );

    return *(this->p[index]);
}


// TIArrayF - extends TIArray with Find() function
//    This allows TIArray to store pointers to types wich does not have
//   assignment and '<', '>' operators
template <class T>
class TIArrayF:
            public TIArray< T >
{
public:
    TIArrayF(size_t size=0, size_t gran=20):
            TIArray<T>(size, gran)
    {}

    virtual ~TIArrayF()
    {}

    int Find(const T&) const;
};

template <class T>
int
TIArrayF<T>::Find(const T& object) const
{
    for( size_t ii=0; ii<this->GetCount(); ii++ )
        if( object == this->elem(ii) )
            return ii;
    return -1;
}

#endif
// End of array.h

