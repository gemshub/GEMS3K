//--------------------------------------------------------------------
// $Id$
//
/// \file nodearray.cpp
/// Implementation of TNodeArray class functionality - advanced
/// interface between GEM IPM and FMT node array
/// working with one DATACH structure and arrays of DATABR structures
//
// Copyright (c) 2004-2012 S.Dmytriyeva, D.Kulik
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
//

#include <cmath>
#include <algorithm>
#include "v_detail.h"

#ifdef NODEARRAYLEVEL

#include "nodearray.h"
#ifndef NOPARTICLEARRAY
#include "particlearray.h"
#endif

TNodeArray* TNodeArray::na;

// ------------------------------------------------------------------

bool TNodeArray::NeedGEMS( TNode* wrkNode, const TestModeGEMParam& modeParam, DATABR* C0, DATABR* C1  )
{
    bool NeedGEM = false;
    DATACH* CH = wrkNode->pCSD();  // DataCH structure
    double dc;

    if( modeParam.useSIA == S_OFF )
        NeedGEM = true;
    else
    {   C1->IterDone = 0;
        NeedGEM = false;
    }

    // Here we compare this node for current time and for previous time - works for AIA and PIA
    for( long int ic=0; ic < CH->nICb; ic++)    // do we check charge here?
    {
        // It has to be checked on minimal allowed c0 value
        if( C1->bIC[ic] < modeParam.cez )
            C1->bIC[ic] = modeParam.cez; // to prevent loss of Independent Component

        dc = C0->bIC[ic] - C1->bIC[ic];
        if( fabs( dc ) > std::min( modeParam.cdv, (C1->bIC[ic] * 1e-3 ) ))
        {
            NeedGEM = true;  // we still need to recalculate equilibrium
            // in this node because its vector b has changed
        }
    }
    C1->bIC[CH->nICb-1] = 0.;   // zeroing charge off in bulk composition

    return NeedGEM;
}

long int TNodeArray::SmartMode( const TestModeGEMParam& modeParam, long int ii, bool* iaN  )
{
    long int Mode = NEED_GEM_AIA;
    //bool* iaN = piaNode();     // indicators for IA in the nodes

    if( modeParam.useSIA == S_OFF )
        iaN[ii] = true;
    else
        iaN[ii] = false;

    if( modeParam.mode == NEED_GEM_SIA )
    {
        // smart algorithm
        if( iaN[ii] == true )
        {
            Mode = NEED_GEM_AIA;
        }
        else {
            Mode = NEED_GEM_SIA;
            if( modeParam.useSIA == S_ON )   // force loading of primal solution into GEMIPM
                Mode *= -1;            // othervise use internal (old) primal solution
        }
    }

    return Mode;
}

std::string TNodeArray::ErrorGEMsMessage( long int RetCode,  long int ii, long int step  )
{
    std::string err_msg;
    char buf[200];

    sprintf( buf, " Node= %-8ld  Step= %-8ld\n", ii, step );
    err_msg = buf;

    switch( RetCode )
    {
    case BAD_GEM_AIA:
        err_msg += "Bad GEM result using LPP AIA";
        break;
    case  ERR_GEM_AIA:
        err_msg += "GEM calculation error using LPP AIA";
        break;
    case  BAD_GEM_SIA:
        err_msg += "Bad GEM result using SIA";
        break;
    case  ERR_GEM_SIA:
        err_msg += "GEM calculation error using SIA";
        break;
    case  T_ERROR_GEM:  err_msg +=  "Terminal error in GEMS3K module";
    }

    return  err_msg;
}


//-------------------------------------------------------------------
// setNodeArray()
// Copying data from work DATABR structure into the node array
// (as specified in nodeTypes array, ndx index of dataBR files in
//    the ipmfiles_lst_name list).
//
//-------------------------------------------------------------------
void  TNodeArray::setNodeArray( long int ndx, long int* nodeTypes  )
{

    for( long int ii=0; ii<anNodes; ii++)
        if(  (!nodeTypes && ndx==0) ||
             ( nodeTypes && (nodeTypes[ii] == ndx/*i+1*/ )) )
        {
            calcNode->pCNode()->NodeHandle = ndx/*(i+1)*/;
            NodT0[ii] = allocNewDBR( calcNode);
            NodT1[ii] = allocNewDBR( calcNode);

            MoveWorkNodeToArray( calcNode, ii, anNodes, NodT0);
            //  CopyWorkNodeFromArray( calcNode, ii, anNodes,NodT0);
            MoveWorkNodeToArray( calcNode, ii, anNodes, NodT1);
            //  CopyWorkNodeFromArray( calcNode, ii, anNodes,NodT1);
        }
}

//---------------------------------------------------------//


// Copying data for node ii from node array into work DATABR structure
//
void TNodeArray::CopyWorkNodeFromArray( TNode* wrkNode, long int ii, long int nNodes, DATABRPTR* arr_BR )
{
    // from arr_BR[ii] to pCNode() structure
    if( ii < 0 || ii>= nNodes )
        return;
    // memory must be allocated before

    // mem_cpy( &wrkNode->pCNode()->NodeHandle, &arr_BR[ii]->NodeHandle, 6*sizeof(short));
    wrkNode->pCNode()->NodeHandle = arr_BR[ii]->NodeHandle;
    wrkNode->pCNode()->NodeTypeHY = arr_BR[ii]->NodeTypeHY;
    wrkNode->pCNode()->NodeTypeMT = arr_BR[ii]->NodeTypeMT;
    wrkNode->pCNode()->NodeStatusFMT = arr_BR[ii]->NodeStatusFMT;
    wrkNode->pCNode()->NodeStatusCH = arr_BR[ii]->NodeStatusCH;
    wrkNode->pCNode()->IterDone = arr_BR[ii]->IterDone;      //6
    // mem_cpy( &wrkNode->pCNode()->TK, &arr_BR[ii]->TK, 32*sizeof(double));
    wrkNode->pCNode()->TK = arr_BR[ii]->TK;
    wrkNode->pCNode()->P = arr_BR[ii]->P;
    wrkNode->pCNode()->Vs = arr_BR[ii]->Vs;
    wrkNode->pCNode()->Vi = arr_BR[ii]->Vi;
    wrkNode->pCNode()->Ms = arr_BR[ii]->Ms;
    wrkNode->pCNode()->Mi = arr_BR[ii]->Mi;
    wrkNode->pCNode()->Gs = arr_BR[ii]->Gs;
    wrkNode->pCNode()->Hs = arr_BR[ii]->Hs;
    wrkNode->pCNode()->Hi = arr_BR[ii]->Hi;
    wrkNode->pCNode()->IC = arr_BR[ii]->IC;
    wrkNode->pCNode()->pH = arr_BR[ii]->pH;
    wrkNode->pCNode()->pe = arr_BR[ii]->pe;
    wrkNode->pCNode()->Eh = arr_BR[ii]->Eh; //13

    wrkNode->pCNode()->Tm = arr_BR[ii]->Tm;
    wrkNode->pCNode()->dt = arr_BR[ii]->dt;
#ifdef NODEARRAYLEVEL
    wrkNode->pCNode()->Dif = arr_BR[ii]->Dif;
    wrkNode->pCNode()->Vt = arr_BR[ii]->Vt;
    wrkNode->pCNode()->vp = arr_BR[ii]->vp;
    wrkNode->pCNode()->eps = arr_BR[ii]->eps;
    wrkNode->pCNode()->Km = arr_BR[ii]->Km;
    wrkNode->pCNode()->Kf = arr_BR[ii]->Kf;
    wrkNode->pCNode()->S = arr_BR[ii]->S;
    wrkNode->pCNode()->Tr = arr_BR[ii]->Tr;
    wrkNode->pCNode()->h = arr_BR[ii]->h;
    wrkNode->pCNode()->rho = arr_BR[ii]->rho;
    wrkNode->pCNode()->al = arr_BR[ii]->al;
    wrkNode->pCNode()->at = arr_BR[ii]->at;
    wrkNode->pCNode()->av = arr_BR[ii]->av;
    wrkNode->pCNode()->hDl = arr_BR[ii]->hDl;
    wrkNode->pCNode()->hDt = arr_BR[ii]->hDt;
    wrkNode->pCNode()->hDv = arr_BR[ii]->hDv;
    wrkNode->pCNode()->nto = arr_BR[ii]->nto; //19
#endif
    // Dynamic data - dimensions see in DATACH.H and DATAMT.H structures
    // exchange of values occurs through lists of indices, e.g. xDC, xPH
    copyValues( wrkNode->pCNode()->xDC, arr_BR[ii]->xDC, wrkNode->pCSD()->nDCb );
    copyValues( wrkNode->pCNode()->gam, arr_BR[ii]->gam, wrkNode->pCSD()->nDCb );
    if( pCSD()->nAalp >0 )
        copyValues( wrkNode->pCNode()->aPH, arr_BR[ii]->aPH, wrkNode->pCSD()->nPHb );
    else  wrkNode->pCNode()->aPH = nullptr;
    copyValues( wrkNode->pCNode()->xPH, arr_BR[ii]->xPH, wrkNode->pCSD()->nPHb );
    copyValues( wrkNode->pCNode()->omPH, arr_BR[ii]->omPH, wrkNode->pCSD()->nPHb );
    copyValues( wrkNode->pCNode()->vPS, arr_BR[ii]->vPS, wrkNode->pCSD()->nPSb );
    copyValues( wrkNode->pCNode()->mPS, arr_BR[ii]->mPS, wrkNode->pCSD()->nPSb );

    copyValues( wrkNode->pCNode()->bPS, arr_BR[ii]->bPS,
                wrkNode->pCSD()->nPSb*wrkNode->pCSD()->nICb );
    copyValues( wrkNode->pCNode()->xPA, arr_BR[ii]->xPA, wrkNode->pCSD()->nPSb );
    copyValues( wrkNode->pCNode()->dul, arr_BR[ii]->dul, wrkNode->pCSD()->nDCb );
    copyValues( wrkNode->pCNode()->dll, arr_BR[ii]->dll, wrkNode->pCSD()->nDCb );
    copyValues( wrkNode->pCNode()->bIC, arr_BR[ii]->bIC, wrkNode->pCSD()->nICb );
    copyValues( wrkNode->pCNode()->rMB, arr_BR[ii]->rMB, wrkNode->pCSD()->nICb );
    copyValues( wrkNode->pCNode()->uIC, arr_BR[ii]->uIC, wrkNode->pCSD()->nICb );
    copyValues( wrkNode->pCNode()->bSP, arr_BR[ii]->bSP, wrkNode->pCSD()->nICb );
    copyValues( wrkNode->pCNode()->amru, arr_BR[ii]->amru, wrkNode->pCSD()->nPSb );
    copyValues( wrkNode->pCNode()->amrl, arr_BR[ii]->amrl, wrkNode->pCSD()->nPSb );
}

// new Copying data for node iNode back from work DATABR structure into the node array
void TNodeArray::MoveWorkNodeToArray( TNode* wrkNode, long int ii, long int nNodes, DATABRPTR* arr_BR )
{
    // from arr_BR[ii] to pCNode() structure
    if( ii < 0 || ii>= nNodes )
        return;
    // memory must be allocated before

    // mem_cpy( &wrkNode->pCNode()->NodeHandle, &arr_BR[ii]->NodeHandle, 6*sizeof(short));
    arr_BR[ii]->NodeHandle = wrkNode->pCNode()->NodeHandle;
    arr_BR[ii]->NodeTypeHY = wrkNode->pCNode()->NodeTypeHY;
    arr_BR[ii]->NodeTypeMT = wrkNode->pCNode()->NodeTypeMT;
    arr_BR[ii]->NodeStatusFMT = wrkNode->pCNode()->NodeStatusFMT;
    arr_BR[ii]->NodeStatusCH = wrkNode->pCNode()->NodeStatusCH;
    arr_BR[ii]->IterDone = wrkNode->pCNode()->IterDone;      //6
    // mem_cpy( &wrkNode->pCNode()->TK, &arr_BR[ii]->TK, 32*sizeof(double));
    arr_BR[ii]->TK = wrkNode->pCNode()->TK;
    arr_BR[ii]->P = wrkNode->pCNode()->P;
    arr_BR[ii]->Vs = wrkNode->pCNode()->Vs;
    arr_BR[ii]->Vi = wrkNode->pCNode()->Vi;
    arr_BR[ii]->Ms = wrkNode->pCNode()->Ms;
    arr_BR[ii]->Mi = wrkNode->pCNode()->Mi;
    arr_BR[ii]->Gs = wrkNode->pCNode()->Gs;
    arr_BR[ii]->Hs = wrkNode->pCNode()->Hs;
    arr_BR[ii]->Hi = wrkNode->pCNode()->Hi;
    arr_BR[ii]->IC = wrkNode->pCNode()->IC;
    arr_BR[ii]->pH = wrkNode->pCNode()->pH;
    arr_BR[ii]->pe = wrkNode->pCNode()->pe;
    arr_BR[ii]->Eh = wrkNode->pCNode()->Eh; //13

    arr_BR[ii]->Tm = wrkNode->pCNode()->Tm;
    arr_BR[ii]->dt = wrkNode->pCNode()->dt;
#ifdef NODEARRAYLEVEL
    arr_BR[ii]->Dif = wrkNode->pCNode()->Dif;
    arr_BR[ii]->Vt = wrkNode->pCNode()->Vt;
    arr_BR[ii]->vp = wrkNode->pCNode()->vp;
    arr_BR[ii]->eps = wrkNode->pCNode()->eps;
    arr_BR[ii]->Km = wrkNode->pCNode()->Km;
    arr_BR[ii]->Kf = wrkNode->pCNode()->Kf;
    arr_BR[ii]->S = wrkNode->pCNode()->S;
    arr_BR[ii]->Tr = wrkNode->pCNode()->Tr;
    arr_BR[ii]->h = wrkNode->pCNode()->h;
    arr_BR[ii]->rho = wrkNode->pCNode()->rho;
    arr_BR[ii]->al = wrkNode->pCNode()->al;
    arr_BR[ii]->at = wrkNode->pCNode()->at;
    arr_BR[ii]->av = wrkNode->pCNode()->av;
    arr_BR[ii]->hDl = wrkNode->pCNode()->hDl;
    arr_BR[ii]->hDt = wrkNode->pCNode()->hDt;
    arr_BR[ii]->hDv = wrkNode->pCNode()->hDv;
    arr_BR[ii]->nto = wrkNode->pCNode()->nto; //19
#endif
    // Dynamic data - dimensions see in DATACH.H and DATAMT.H structures
    // exchange of values occurs through lists of indices, e.g. xDC, xPH
    copyValues( arr_BR[ii]->xDC, wrkNode->pCNode()->xDC, wrkNode->pCSD()->nDCb );
    copyValues( arr_BR[ii]->gam, wrkNode->pCNode()->gam, wrkNode->pCSD()->nDCb );
    if( wrkNode->pCSD()->nAalp >0 )
        copyValues( arr_BR[ii]->aPH, wrkNode->pCNode()->aPH, wrkNode->pCSD()->nPHb );
    else  arr_BR[ii]->aPH = nullptr;
    copyValues( arr_BR[ii]->xPH, wrkNode->pCNode()->xPH, wrkNode->pCSD()->nPHb );
    copyValues( arr_BR[ii]->omPH, wrkNode->pCNode()->omPH, wrkNode->pCSD()->nPHb );
    copyValues( arr_BR[ii]->vPS, wrkNode->pCNode()->vPS, wrkNode->pCSD()->nPSb );
    copyValues( arr_BR[ii]->mPS, wrkNode->pCNode()->mPS, wrkNode->pCSD()->nPSb );

    copyValues( arr_BR[ii]->bPS, wrkNode->pCNode()->bPS,
                wrkNode->pCSD()->nPSb*wrkNode->pCSD()->nICb );
    copyValues( arr_BR[ii]->xPA, wrkNode->pCNode()->xPA, wrkNode->pCSD()->nPSb );
    copyValues( arr_BR[ii]->dul, wrkNode->pCNode()->dul, wrkNode->pCSD()->nDCb );
    copyValues( arr_BR[ii]->dll, wrkNode->pCNode()->dll, wrkNode->pCSD()->nDCb );
    copyValues( arr_BR[ii]->bIC, wrkNode->pCNode()->bIC, wrkNode->pCSD()->nICb );
    copyValues( arr_BR[ii]->rMB, wrkNode->pCNode()->rMB, wrkNode->pCSD()->nICb );
    copyValues( arr_BR[ii]->uIC, wrkNode->pCNode()->uIC, wrkNode->pCSD()->nICb );
    copyValues( arr_BR[ii]->bSP, wrkNode->pCNode()->bSP, wrkNode->pCSD()->nICb );
    copyValues( arr_BR[ii]->amru, wrkNode->pCNode()->amru, wrkNode->pCSD()->nPSb );
    copyValues( arr_BR[ii]->amrl, wrkNode->pCNode()->amrl, wrkNode->pCSD()->nPSb );
}



/** old Copying data for node iNode back from work DATABR structure into the node array
void TNodeArray::MoveWorkNodeToArray( TNode& wrkNode, long int ii, long int nNodes, DATABRPTR* arr_BR )
{
  if( ii < 0 || ii>= nNodes )
    return;
  if( arr_BR[ii] )
  {
       arr_BR[ii] = wrkNode.databr_free( arr_BR[ii] );
       // delete[] arr_BR[ii];
  }
  arr_BR[ii] = wrkNode.pCNode();
// alloc new memory
  wrkNode.allocNewDBR();
}*/

void TNodeArray::CopyNodeFromTo( long int ndx, long int nNod,
                                 DATABRPTR* arr_From, DATABRPTR* arr_To )
{
    if( !arr_From || !arr_To )
        return;
    CopyWorkNodeFromArray( calcNode, ndx, nNod, arr_From );
    MoveWorkNodeToArray( calcNode, ndx,  nNod, arr_To );
}

//---------------------------------------------------------
// Methods for working with node arrays (access to data from DBR)

// Calculate phase (carrier) mass, kg  of single component phase
double TNodeArray::get_mPH( long int ia, long int nodex, long int PHx )
{
    long int DCx = calcNode->Phx_to_DCx( Ph_xDB_to_xCH(PHx) );
    double val=0.;

    if( DCx >= pCSD()->nDCs && DCx < pCSD()->nDC )
    {
        val = pCSD()->DCmm[DCx];
        if( ia == 0)
            val *= pNodT0()[nodex]->xDC[calcNode->DC_xCH_to_xDB(DCx)];
        else
            val *= pNodT1()[nodex]->xDC[calcNode->DC_xCH_to_xDB(DCx)];
    }

    return val;
}

// Calculate phase volume (in cm3) of single - component phase
double TNodeArray::get_vPH( long int ia, long int nodex, long int PHx )
{
    long int DCx = calcNode->Phx_to_DCx( Ph_xDB_to_xCH(PHx) );
    double val=0.;

    if( DCx >= pCSD()->nDCs && DCx < pCSD()->nDC )
    {
        double T, P;
        if( ia == 0 )
        {
            T = pNodT0()[(nodex)]->TK;
            P = pNodT0()[(nodex)]->P;
            val = pNodT0()[nodex]->xDC[calcNode->DC_xCH_to_xDB(DCx)]; // number of moles
        }
        else
        {
            T = pNodT1()[(nodex)]->TK;
            P = pNodT1()[(nodex)]->P;
            val = pNodT1()[nodex]->xDC[calcNode->DC_xCH_to_xDB(DCx)];
        }
        val *= calcNode->DC_V0( DCx, P, T );
    }
    return val;
}


// Calculate bulk compositions  of single component phase
double TNodeArray::get_bPH( long int ia, long int nodex, long int PHx, long int ICx )
{
    long int DCx = calcNode->Phx_to_DCx( Ph_xDB_to_xCH(PHx) );
    double val=0.;

    if( DCx >= pCSD()->nDCs && DCx < pCSD()->nDC )
    {
        val = pCSD()->A[ pCSD()->xic[ICx] + DCx * pCSD()->nIC];
        if( ia == 0)
            val *= pNodT0()[nodex]->xDC[calcNode->DC_xCH_to_xDB(DCx)];
        else
            val *= pNodT1()[nodex]->xDC[calcNode->DC_xCH_to_xDB(DCx)];
    }

    return val;
}

//---------------------------------------------------------
// working with grid

// Set grid coordinate array use predefined array aGrid
// or set up regular scale
void TNodeArray::SetGrid( double aSize[3], double (*aGrid)[3] )
{
    long int i, j, k, ndx;
    size.x =  aSize[0];
    size.y =  aSize[1];
    size.z =  aSize[2];

    LOCATION delta( size.x/sizeN, size.y/sizeM, size.z/sizeK );
    if( !grid )
        grid = new LOCATION[ anNodes ];

    for( i = 0; i < sizeN; i++ )
        for( j = 0; j < sizeM; j++ )
            for( k = 0; k < sizeK; k++ )
            {
                ndx = iNode( i, j, k );
                if( aGrid )
                {
                    grid[ndx].x = aGrid[ndx][0];
                    grid[ndx].y = aGrid[ndx][1];
                    grid[ndx].z = aGrid[ndx][2];
                }
                else
                {
                    grid[ndx].x =delta.x*i;
                    grid[ndx].y =delta.y*j;
                    grid[ndx].z =delta.z*k;
                }
            }
    // outside limits settted in size
}

// test location in node
bool TNodeArray::isLocationInNode( long int iNode, LOCATION cxyz ) const
{
    long int i1, j1, k1;
    i1 = indN( iNode );
    j1 = indM( iNode );
    k1 = indK( iNode );
    return isLocationInNode( i1,j1, k1, cxyz ) ;
}

// test location in node
bool TNodeArray::isLocationInNode( long int ii, long int jj, long int kk, LOCATION cxyz ) const
{
    LOCATION maxl;

    if( ii<0 || ii>= sizeN || jj<0 ||
            jj >= sizeM || kk <0 || kk >= sizeK )
        return false;

    long int ndx = iNode( ii, jj, kk );
    maxl = getGrid( ii+1, jj+1, kk+1 ); // only for rectangular
    // must be changed
    //  x = const, find new y,z srez i
    // analiz pryamougol`nika pri y1 == const, poisk z21 i z22
    // analiz otrezka po z2
    if( grid[ndx].x <= cxyz.x &&
            ( cxyz.x < maxl.x || ( cxyz.x <= maxl.x && essentiallyEqual( size.x,maxl.x ) ) )&&
            grid[ndx].y <= cxyz.y &&  cxyz.y <= maxl.y &&
            grid[ndx].z <= cxyz.z &&  cxyz.z <= maxl.z )
        return true;

    return false; // location behind the node
}

// Finds a node absolute index for the current
// point location (uses grid coordinate array grid[])
// performance-important functions to be used e.g. in particle tracking methods
long int TNodeArray::FindNodeFromLocation( LOCATION cxyz, long int old_node ) const
{
    //  LOCATION maxl;
    long int i, j, k/*, ndx*/;

    if( old_node == -1 )
    { // check all nodes
        for( i = 0; i < sizeN; i++ )
            for( j = 0; j < sizeM; j++ )
                for( k = 0; k < sizeK; k++ )
                {
                    if(  isLocationInNode( i, j, k, cxyz ) )
                        return iNode( i, j, k );
                }
    }
    else // check only nearest nodes
    {
        long int i1, j1, k1;
        i1 = indN( old_node );
        j1 = indM( old_node );
        k1 = indK( old_node );

        for( i = i1-1; i <= i1+1; i++ )
            for( j = j1-1; j <= j1+1; j++ )
                for( k = k1-1; k <= k1+1; k++ )
                {
                    if(  isLocationInNode( i, j, k, cxyz ) )
                        return iNode( i, j, k );
                }
    }
    return -1; // behind region
}

// get current node location
// if iN, jN or kN more then corresponding sizeN, sizeM, sizeK
// return size of system
LOCATION TNodeArray::getGrid( long int iN, long int jN, long int kN ) const
{
    LOCATION loc;
    long int i1, j1, k1;

    // only for test
    if( iN < 0 || iN > sizeN ||
            jN < 0 || jN > sizeM ||
            kN < 0 || kN > sizeK  )
        Error( "", "getGrid - programm error");

    if( iN == sizeN )
        i1 = iN-1;
    else i1 = iN;
    if( jN == sizeM )
        j1 = jN-1;
    else j1 = jN;
    if( kN == sizeK )
        k1 = kN-1;
    else k1 = kN;

    loc = grid[ iNode( i1, j1, k1)];
    if( i1 != iN )
        loc.x = size.x;
    if( j1 != jN ) loc.y = size.y;
    if( k1 != kN ) loc.z = size.z;

    return loc;
}

// get 3D sizes for node (  from cxyz[0] - to cxyz[1] )
// only for rectangular -  must be changed
// for any must be LOCATION cxyz[8]
void TNodeArray::GetNodeSizes( long int ndx, LOCATION cxyz[2] )
{
    LOCATION maxl;
    long int i, j, k;

    i = indN( ndx );
    j = indM( ndx );
    k = indK( ndx );

    cxyz[0] = grid[ndx];
    cxyz[1] = getGrid( i+1, j+1, k+1 );
    /*
  cxyz[1] = getGrid( i, j, k+1 );
  cxyz[2] = getGrid( i, j+1, k );
  cxyz[3] = getGrid( i, j+1, k+1 );
  cxyz[4] = getGrid( i+1, j, k );
  cxyz[5] = getGrid( i+1, j, k+1 );
  cxyz[6] = getGrid( i+1, j+1, k );
  cxyz[7] = getGrid( i+1, j+1, k+1 );
  */
}

// get full mass of all particles of the type ptype in node ndx
// ndx    -  (absolute) index of the node
// ptype  -  particle type index ( 1 to 255 )
// tcode  -  particle transport mechanism code (see enum PTCODE)
// ips   - DataBr index of phase or species to which this particle is connected
double TNodeArray::GetNodeMass( long int ndx,
                                char /*ptype*/, char tcode, unsigned char iips )
{
    double mass = 0.;
    //   DATABR* dbr = NodT1[ndx];
    DATACH* dch = pCSD();
    long int xWatCH, ips = (long int)iips;

    switch( tcode )
    {
    case DISSOLVED: // mass of dissolved matter in aqueous solution
        xWatCH = dch->nDCinPH[dch->xph[0]]-1; // CH index of water
        //                       mass = node1_mPS(ndx,ips); // - node1_xPA(ndx,ips)*dch->DCmm[xWatCH];
        mass = node1_xPA(ndx,ips)*dch->DCmm[xWatCH]; // Mass of aq-solvent
        break;
    case ADVECTIVE: // mass of aq solution
        mass = node1_mPS( ndx, ips );
        break;
    case COLLOID:  // mass of phase - solid solution, sorption or pure solid
        if( ips < dch->nPSb )
            mass = node1_mPS( ndx, ips );
        else
            mass = node1_mPH( ndx, ips );
        break;
    case DIFFUSIVE: // mass of the diffusing species
        mass = DCmm( ips ) * node1_xDC( ndx, ips );
        break;
    default: break;
    }
    return mass;
}

// move mass m_v from node ndx_from to node ndx_to, one particle move
// ndx_from    -  (absolute) index of the old node
// ndx_to     -  (absolute) index of the new  node
// type  -  particle type index ( 1 to 255 )
// COmpMode: true: transport of DCs; false - transport of ICs 
// tcode  -  particle transport mechanism code (see enum PTCODE)
// iips   - DataBr index of phase or species to which this particle is connected
// m_v -  mass or volume of the particle (depending on ptype and mmode)
void TNodeArray::MoveParticleMass( long int ndx_from, long int ndx_to,
                                   char /*type*/, char CompMode, char tcode, unsigned char iips, double m_v )
{
    double mass = 0., coeff, mol, mWat=0., fmolal=1., aji;
    DATABR* dbr = NodT1[ndx_from];
    DATACH* dch = pCSD();
    long int xWatCH=0, ic, ips = (long int)iips;
    if( tcode == DISSOLVED || tcode == ADVECTIVE || tcode == DIFFUSIVE )
    {
        xWatCH = pCSD()->nDCinPH[pCSD()->xph[0]]-1; // CH index of water
        //	   mWat = node1_xDC( ndx_from, xWatCH )* CSD->DCmm[xWatCH];
        mWat = node1_xPA(ndx_from, ips) * pCSD()->DCmm[xWatCH];  // Mass of water-solvent
        fmolal = 1.0; // 1000./mWat;              // molality conversion factor
    }

    switch( tcode )
    {
    case DISSOLVED: // mass of dissolved matter in aqueous solution
        mass = mWat;  // trying normalization over mass of water-solvent
        //    				mass = 1000/fmolal; // grams of water in the node
        //    				fmolal = 1.0;
        //    	mass = node1_mPS(ndx_from,ips) - mWat;
        break;
    case ADVECTIVE: // mass of liquid phase for full advection
        mass = node1_mPS( ndx_from, ips );
        //  - node1_xPA(ndx_from,ips) * dch->DCmm[xWatCH];
        break;
    case COLLOID:  // colloid particles
        if( ips < dch->nPSb )
            mass = node1_mPS( ndx_from, ips );
        else
            mass = node1_mPH( ndx_from, ips );
        break;
    case DIFFUSIVE: // gets the mass of diffusing species in the phase
        mass = DCmm( ips ) * node1_xDC( ndx_from, ips );
        break;
    }
    coeff = m_v/mass; // mass of particle/mass of phase (solvent). Is this reasonable?

    if( CompMode == true )
    { // Moving dependent components
        for(long int jc=0; jc < pCSD()->nDC; jc++ )
        {
            mol = 0.; // moles of DC transported in the particle
            switch( tcode )
            {
            case DISSOLVED: // moving only dissolved DC (-H2O)
                //             if( jc == xWatCH )
                //            	 continue;  // H2O is ignored - not moved with the particle
            case ADVECTIVE: // moving DC of the whole aq phase
                if( jc > xWatCH )
                    continue;     // ignoring non-aqueous species
                mol = node1_xDC( ndx_from, jc ) * coeff * fmolal;
                break;
                //             ( node1_bPS( ndx_from, ips, ie )
                //                   - nodeCH_A( xWatCH, ie)
                //                   * node1_xPA(ndx_from,ips)) * coeff;
            case COLLOID:  // moving DC of solid particle - to be completed!
                // 	    	 if( ips < dch->nPSb )
                //                  mol = node1_bPS( ndx_from, ips, ie ) * coeff;
                //             else
                //                  mol = node1_bPH( ndx_from, ips, ie ) * coeff;
                break;
            case DIFFUSIVE: // moving DC - a diffusing species
                if( jc != ips )
                    continue;     // ignoring other diffusive species
                mol = node1_xDC( ndx_from, jc ) * coeff * fmolal;
                break;
            }
            if( tcode == DISSOLVED || tcode == ADVECTIVE || tcode == DIFFUSIVE )
                mol /= fmolal;       // back from molality to moles

            if( fabs(mol) > 1e-20 ) // mtp->cdv ) // Threshold for DC change carried over in the particle
            {
                if( NodT1[ndx_from]->NodeTypeHY != NBC3source )
                {
                    node1_xDC( ndx_from, jc ) -= mol;   // Correcting species amount in source node at T1
                    for( ic=0; ic<pCSD()->nICb; ic++)  // incrementing independent components
                    {
                        aji = DCaJI( jc, ic );
                        if( noZero( aji ) )
                            node1_bIC(ndx_from, ic) -= aji * mol;
                    }
                }
                if( ndx_to >= 0 && ndx_to < anNodes )
                {
                    if( NodT1[ndx_to]->NodeTypeHY != NBC3source )
                    {
                        node1_xDC( ndx_to, jc ) += mol;
                        for( ic=0; ic<pCSD()->nICb; ic++)  // incrementing independent components
                        {
                            aji = DCaJI( jc, ic );
                            if( noZero( aji ) )
                                node1_bIC(ndx_to, ic) += aji * mol;
                        }
                    }
                }
                else
                    if(dbr->NodeTypeHY != NBC3sink  && dbr->NodeTypeHY != NBC3source)
                        std::cout << "W002MTRW " << "Warning: Particle jumped outside the domain" << std::endl;
                // 	  			  Error( "W002MTRW", "Warning: Particle jumped outside the domain" );
                //	   }
            }
        } // loop jc
    }
    else {
        // Transport of independent components
        for(long int ie=0; ie < pCSD()->nICb; ie++ )
        {
            mol = 0.; // moles of IC in the particle
            switch( tcode )
            {
            case DISSOLVED: // moving only dissolved IC (-H2O)
                mol = ( node1_bPS( ndx_from, ips, ie )
                        //                              - nodeCH_A( xWatCH, ie)* node1_xPA(ndx_from,ips)
                        ) * coeff * fmolal;
                break;
            case ADVECTIVE: // moving IC of the whole aq phase
                mol = ( node1_bPS( ndx_from, ips, ie )
                        // - nodeCH_A( xWatCH, ie) * node1_xPA(ndx_from,ips)
                        ) * coeff * fmolal;
                break;
            case COLLOID:   // moving IC of solid particle
                if( ips < dch->nPSb )
                    mol = node1_bPS( ndx_from, ips, ie ) * coeff;
                else
                    mol = node1_bPH( ndx_from, ips, ie ) * coeff;
                break;
            case DIFFUSIVE: // moving IC of diffusing species
                mol = node1_xDC( ndx_from, ips ) * DCaJI( ips, ie )
                        * coeff * fmolal;
                break;
            }
            if( tcode == DISSOLVED || tcode == ADVECTIVE || tcode == DIFFUSIVE )
                mol /= fmolal;
            if( dbr->NodeTypeHY != NBC3source )
                dbr->bIC[ie] -= mol;
            if( ndx_to >= 0 && ndx_to < anNodes )
            {
                if( NodT1[ndx_to]->NodeTypeHY != NBC3source )
                    NodT1[ndx_to]->bIC[ie] += mol;
            }
            else
                if(dbr->NodeTypeHY != NBC3sink  && dbr->NodeTypeHY != NBC3source)
                    std::cout << "W002MTRW " << "Warning: Particle jumped outside the domain" << std::endl;
            //        	 Error( "W002MTRW", "Warning: Particle jumped outside the domain" );

        } // loop ie
    } // else
    // End of function
}

//==========================================================================
// Data collection for monitoring differences

// Prints difference increments in a all nodes (cells) for time point t / at
//
void TNodeArray::logDiffsIC( FILE* diffile, long int t, double at, long int nx, long int every_t )
{
    double dc;
    long int i, ie;

    if( t % every_t )
        return;

    fprintf( diffile, "\nStep= %-8ld  Time= %-12.4g\nNode#   ", t, at );
    for( ie=0; ie < (pCSD()->nICb); ie++ )
        fprintf( diffile, "%-12.4s ", pCSD()->ICNL[ie] );
    for (i=0; i<nx; i++)    // node iteration
    {
        fprintf( diffile, "\n%5ld   ", i );
        for( ie=0; ie < (pCSD()->nICb); ie++ )
        {
            dc = NodT1[i]->bIC[ie] - NodT0[i]->bIC[ie];
            fprintf( diffile, "%-12.4g ", dc );
        }
    }
    fprintf( diffile, "\n" );
}

// Data collection for monitoring 1D profiles in debugging FMT models
//
// Prints dissolved species molarities in all cells for time point t / at
//
void TNodeArray::logProfileAqDC( FILE* logfile, long int t, double at, long int nx, long int every_t )
{
    double pm;
    long int i, is;
    if( t % every_t )
        return;
    fprintf( logfile, "\nStep= %-8ld\tTime= %-12.4g, s\tDissolved species concentrations, M\n", t, at );
    fprintf(logfile, "%s","Node#   ");
    for( is=0; is < (pCSD()->nDCb); is++ )
        fprintf( logfile, "%-12.4s ", pCSD()->DCNL[is] );
    for (i=0; i<nx; i++)    // node iteration
    {
        fprintf( logfile, "\n%5ld   ", i );
        for( is=0; is < (pCSD()->nDCinPH[0]); is++ )
        {
            pm = NodT1[i]->xDC[is]/NodT1[i]->vPS[0]/1000.;  // Assumes there is aq phase!
            // dissolved species molarity
            fprintf( logfile, "%-12.4g ", pm );
        }
    }
    fprintf( logfile, "\n" );
}

// Prints dissolved elemental molarities in all cells for time point t / at
//
void TNodeArray::logProfileAqIC( FILE* logfile, long int t, double at, long int nx, long int every_t )
{
    double pm;
    long int i, ie;
    if( t % every_t )
        return;
    fprintf( logfile, "\nStep= %-8ld\tTime= %-12.4g,s\tDissolved IC total concentrations, M\n", t, at );
    fprintf(logfile, "%s","Node#   ");
    for( ie=0; ie < (pCSD()->nICb); ie++ )
        fprintf( logfile, "%-12.4s ", pCSD()->ICNL[ie] );
    for (i=0; i<nx; i++)    // node iteration
    {
        fprintf( logfile, "\n%5ld   ", i );
        for( ie=0; ie < (pCSD()->nICb); ie++ )
        {
            pm = NodT1[i]->bPS[ie]/NodT1[i]->vPS[0]/1000.;  // Assumes there is aq phase!
            // total dissolved element molarity
            fprintf( logfile, "%-12.4g ", pm );
        }
    }
    fprintf( logfile, "\n" );
}

// Data collection for monitoring 1D profiles
// Prints total elemental amounts in all cells for time point t / at
//
void TNodeArray::logProfileTotIC( FILE* logfile, long int t, double at, long int nx, long int every_t )
{
    double pm;
    long int i, ie;
    if( t % every_t )
        return;
    fprintf( logfile, "\nStep= %-8ld\tTime= %-12.4g,s\tBulk IC amounts, moles\n", t, at );
    fprintf(logfile, "%s","Node#   ");
    for( ie=0; ie < (pCSD()->nICb); ie++ )
        fprintf( logfile, "%-12.4s ", pCSD()->ICNL[ie] );
    for (i=0; i<nx; i++)    // node iteration
    {
        fprintf( logfile, "\n%5ld   ", i );
        for( ie=0; ie < (pCSD()->nICb); ie++ )
        {
            pm = NodT1[i]->bIC[ie];
            fprintf( logfile, "%-12.4g ", pm );
        }
    }
    fprintf( logfile, "\n" );
}

// Prints amounts of reactive phases in all cells for time point t / at
void TNodeArray::logProfilePhMol( FILE* logfile, TParticleArray* pa, long int t, double at, long int nx, long int every_t )
{
    double pm;
    long int i, ip;
    if( t % every_t )
        return;
    fprintf( logfile, "\nStep= %-8ld\tTime= %-12.4g,s\tAmounts of reactive phases, moles\n", t, at );
    fprintf(logfile, "%s","Node#   ");
    for( ip=0; ip < (pCSD()->nPHb); ip++ )
        fprintf( logfile, "%-12.12s ", pCSD()->PHNL[ip] );
    for (i=0; i<nx; i++)    // node iteration
    {
        fprintf( logfile, "\n%5ld   ", i );
        for( ip=0; ip < (pCSD()->nPHb); ip++ )
        {
            //       pm = NodT1[i]->xPH[ip];
            pm = node1_xPH( i, ip );
            fprintf( logfile, "%-12.4g ", pm );
        }
#ifndef NOPARTICLEARRAY
        if( pa )
            pa->logProfilePhMol( logfile, i );
#endif
    }
    fprintf( logfile, "\n" );
}

// Prints volumes of reactive phases in all cells for time point t / at
// in nodearray layer C1
//
void TNodeArray::logProfilePhVol( FILE* logfile, long int t, double at, long int nx, long int every_t )
{
    double pm;
    long int i, ip;
    if( t % every_t )
        return;
    fprintf( logfile, "\nStep= %-8ld\tTime= %-12.4g,s\tVolumes of reactive phases, moles\n", t, at );
    fprintf(logfile, "%s","Node#   ");
    for( ip=0; ip < (pCSD()->nPHb); ip++ )
        fprintf( logfile, "%-12.12s ", pCSD()->PHNL[ip] );
    for (i=0; i<nx; i++)    // node iteration
    {
        fprintf( logfile, "\n%5ld  ", i );
        for( ip=0; ip < (pCSD()->nPSb); ip++ )
        {   // Multi-component phases
            pm = node1_vPS( i, ip );
            fprintf( logfile, "%-12.4g ", pm );
        }
        for( ip=(pCSD()->nPSb); ip < (pCSD()->nPHb); ip++ )
        {  // Single-component phases
            pm = node1_vPH( i, ip );
            fprintf( logfile, "%-12.4g ", pm );
        }
    }
    fprintf( logfile, "\n" );
}

void TNodeArray::databr_to_vtk( std::fstream& ff, const char*name, double time, long int  cycle,
                                long int  nFilds, long int  (*Flds)[2])
{
    bool all = false;
    long int kk, ii, nf, nel, nel2;
    long int i, j,k;

    // write header of file
    kk = sizeM;
    if(sizeM==1 && sizeK==1) // 05.12.2012 workaround for 2D paraview
        kk=2;
    calcNode->databr_head_to_vtk( ff, name, time, cycle, sizeN, kk, sizeK );

    if( nFilds < 1 || !Flds )
    {  all = true;
        nFilds = 51;
    }

    for( kk=0; kk<nFilds; kk++)
    {
        if( all )
            nf = kk;
        else
            nf= Flds[kk][0];

        calcNode->databr_size_to_vtk(  nf, nel, nel2 );

        if( all )
        { ii=0;}
        else
        { ii = Flds[kk][1];
            nel = ii+1;
        }

        for( ; ii<nel; ii++ )
        {
            calcNode->databr_name_to_vtk( ff, nf, ii, nel2 );

            // cycle for TNode array
            for( i = 0; i < sizeN; i++ )
                for( j = 0; j < sizeM; j++ )
                    for( k = 0; k < sizeK; k++ )
                    {
                        int ndx = iNode( i, j, k );
                        CopyWorkNodeFromArray( calcNode, ndx, anNodes,  pNodT0() );
                        calcNode->databr_element_to_vtk( ff, calcNode->pCNode()/*pNodT0()[(ndx)]*/, nf, ii );
                    }
            if( sizeM==1 && sizeK==1)  // 05.12.2012 workaround for 2D paraview
            {
                for( i = 0; i < sizeN; i++ )
                    for( j = 0; j < sizeM; j++ )
                        for( k = 0; k < sizeK; k++ )
                        {
                            int ndx = iNode( i, j, k );
                            CopyWorkNodeFromArray( calcNode, ndx, anNodes,  pNodT0() );
                            calcNode->databr_element_to_vtk( ff, calcNode->pCNode()/*pNodT0()[(ndx)]*/, nf, ii );
                        }
            }
        }
    }
}

#endif
//-----------------------End of nodearray.cpp--------------------------

