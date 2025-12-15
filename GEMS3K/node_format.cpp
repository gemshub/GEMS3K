//--------------------------------------------------------------------
// $Id$
//
/// \file node_format.cpp
/// Interface for writing DBR to vtk files format
//
// Copyright (c) 2006-2023 S.Dmytriyeva, D.Kulik
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
#include "io_keyvalue.h"
#include "node.h"

namespace  dbr_dch_api {
extern std::vector<io_formats::outField> DataBR_fields;
}

void TNode::databr_element_to_vtk( std::fstream& ff, DATABR *CNode_, long int nfild, long int ndx )
{
    io_formats::KeyValueWrite prar( ff );

    switch( nfild )
    {
    case f_NodeHandle: prar.writeValue( CNode_->NodeHandle);
        break;
    case f_NodeTypeHY: prar.writeValue( CNode_->NodeTypeHY);
        break;
    case f_NodeTypeMT: prar.writeValue( CNode_->NodeTypeMT);
        break;
    case f_NodeStatusFMT: prar.writeValue( CNode_->NodeStatusFMT);
        break;
    case f_NodeStatusCH: prar.writeValue( CNode_->NodeStatusCH);
        break;
    case f_IterDone: prar.writeValue( CNode_->IterDone);
        break;
    case f_TK: prar.writeValue( CNode_->TK);
        break;
    case f_P: prar.writeValue( CNode_->P);
        break;
    case f_Vs: prar.writeValue( CNode_->Vs);
        break;
    case f_Vi: prar.writeValue( CNode_->Vi);
        break;
    case f_Ms: prar.writeValue( CNode_->Ms);
        break;
    case f_Mi: prar.writeValue( CNode_->Mi);
        break;
    case f_Hs: prar.writeValue( CNode_->Hs);
        break;
    case f_Hi: prar.writeValue( CNode_->Hi);
        break;
    case f_Gs: prar.writeValue( CNode_->Gs);
        break;
    case f_IS: prar.writeValue( CNode_->IC);
        break;
    case f_pH:prar.writeValue( CNode_->pH);
        break;
    case f_pe: prar.writeValue( CNode_->pe);
        break;
    case f_Eh: prar.writeValue( CNode_->Eh);
        break;
    case f_Tm: prar.writeValue( CNode_->Tm);
        break;
    case f_dt:prar.writeValue( CNode_->dt);
        break;
#ifndef NO_NODEARRAYLEVEL
    case f_Dif: prar.writeValue( CNode_->Dif);
        break;
    case f_Vt: prar.writeValue( CNode_->Vt);
        break;
    case f_vp: prar.writeValue( CNode_->vp);
        break;
    case f_eps: prar.writeValue( CNode_->eps);
        break;
    case f_Km:  prar.writeValue( CNode_->Km);
        break;
    case f_Kf:  prar.writeValue( CNode_->Kf);
        break;
    case f_S:  prar.writeValue( CNode_->S);
        break;
    case f_Tr:  prar.writeValue( CNode_->Tr);
        break;
    case f_h:  prar.writeValue( CNode_->h);
        break;
    case f_rho:  prar.writeValue( CNode_->rho);
        break;
    case f_al:  prar.writeValue( CNode_->al);
        break;
    case f_at:  prar.writeValue( CNode_->at);
        break;
    case f_av:  prar.writeValue( CNode_->av);
        break;
    case f_hDl:  prar.writeValue( CNode_->hDl);
        break;
    case f_hDt:  prar.writeValue( CNode_->hDt);
        break;
    case f_hDv: prar.writeValue( CNode_->hDv);
        break;
    case f_nto: prar.writeValue( CNode_->nto);
        break;
#endif
    case f_bIC: prar.writeValue(   CNode_->bIC[ndx] );
        break;
    case f_rMB: prar.writeValue(   CNode_->rMB[ndx] );
        break;
    case f_uIC: prar.writeValue(   CNode_->uIC[ndx] );
        break;
    case f_xDC: prar.writeValue(   CNode_->xDC[ndx] );
        break;
    case f_gam: prar.writeValue(   CNode_->gam[ndx] );
        break;
    case f_dll: prar.writeValue(   CNode_->dll[ndx] );
        break;
    case f_dul: prar.writeValue(   CNode_->dul[ndx] );
        break;
    case f_aPH: prar.writeValue(   CNode_->aPH[ndx] );
        break;
    case f_xPH: prar.writeValue(   CNode_->xPH[ndx] );
        break;
    case f_vPS: prar.writeValue(   CNode_->vPS[ndx] );
        break;
    case f_mPS: prar.writeValue(   CNode_->mPS[ndx] );
        break;
    case f_bPS: prar.writeValue(   CNode_->bPS[ndx] );
        break;
    case f_xPA: prar.writeValue(   CNode_->xPA[ndx] );
        break;
    case f_bSP: prar.writeValue(   CNode_->bSP[ndx] );
        break;
    case f_amru: prar.writeValue(   CNode_->amru[ndx] );
        break;
    case f_amrl: prar.writeValue(   CNode_->amrl[ndx] );
        break;
    case f_omph: prar.writeValue(   CNode_->omPH[ndx] );
        break;
        // CNode_ must be pointer to a work node data bridge structure CNode
    case f_mPH: prar.writeValue(   Ph_Mass(ndx) );
        break;
    case f_vPH: prar.writeValue(   Ph_Volume(ndx) );
        break;
    case f_m_t: prar.writeValue(   Get_mIC( ndx ) );
        break;
    case f_con: prar.writeValue(   Get_cDC(ndx) );
        break;
    case f_mju: prar.writeValue(   Get_muDC(ndx, false ) );
        break;
    case f_lga: prar.writeValue(   Get_aDC(ndx, false ) );
        break;
    default: break;
    }
    ff << " " << std::endl;
}


void TNode::databr_name_to_vtk( std::fstream& ff, long int nfild, long int ndx, long int nx2 )
{
  std::string str="", str2="";
  // full name of data fiels (add index name)
  ff << "SCALARS " << dbr_dch_api:: DataBR_fields[nfild].name.c_str();

  switch( dbr_dch_api::DataBR_fields[nfild].indexation )
  {
    case 1: break;
    case nICbi: str = char_array_to_string( CSD->ICNL[ IC_xDB_to_xCH( ndx ) ],MaxICN );
                break;
    case nDCbi: str = char_array_to_string( CSD->DCNL[ DC_xDB_to_xCH( ndx ) ],MaxDCN );
              break;
    case nPHbi:
    case nPSbi: str = char_array_to_string(CSD->PHNL[ Ph_xDB_to_xCH( ndx ) ],MaxPHN );
            break;
    case nPSbnICbi:
                str = char_array_to_string(  CSD->PHNL[ Ph_xDB_to_xCH( ndx/nx2 ) ],MaxPHN );
                str2 = char_array_to_string( CSD->ICNL[ IC_xDB_to_xCH( ndx%nx2 ) ], MaxICN );
          break;
    default: str = std::string( "UNDEFINED");
  }

  strip(str);
  strip(str2);

  if( !str.empty() )
  { ff << "_" << str.c_str();
    if( !str2.empty() )
      ff << "_" << str2.c_str();
    ff << "_";
  }

  if( nfild < 6)
     ff << " long";
  else
     ff << " double";
  ff << " 1" << std::endl;

  ff << "LOOKUP_TABLE default" << std::endl;
}

void TNode::databr_size_to_vtk(  long int nfild, long int& nel, long int& nel2 )
{

    nel = 1;
    nel2 = 1;
    switch( dbr_dch_api::DataBR_fields[nfild].indexation )
    {
      case 1: break;
      case nICbi: nel = CSD->nICb;
                  break;
      case nDCbi: nel = CSD->nDCb;
                  break;
      case nPHbi: nel = CSD->nPHb;
                  break;
      case nPSbi: nel = CSD->nPSb;
                  break;
      case nPSbnICbi: nel = CSD->nPSb; nel2 = CSD->nICb;
        break;
   }

}

void TNode::databr_head_to_vtk( std::fstream& ff, const char*name, double time, long cycle,
                               long nx, long ny, long nz )
{
 ff << "# vtk DataFile Version 3.0" <<  std::endl;
 ff << "GEM2MT " << name <<  std::endl;
 ff << "ASCII" <<  std::endl;
 ff << "DATASET STRUCTURED_POINTS" <<  std::endl;
 ff << "DIMENSIONS " <<  nx << " " << ny << " " << nz << std::endl;
 ff << "ORIGIN 0.0 0.0 0.0" <<  std::endl;
 ff << "SPACING 1.0 1.0 1.0" <<  std::endl;
 //????
 ff << "FIELD TimesAndCycles 2" << std::endl;
 ff << "TIME 1 1 double" <<  std::endl << time << std::endl;
 ff << "CYCLE 1 1 long" <<  std::endl << cycle << std::endl;
 ff << "POINT_DATA " << nx*ny*nz << std::endl;
}

void TNode::databr_to_vtk( std::fstream& ff, const char*name, double time, long int  cycle,
                          long int  nFilds, long int  (*Flds)[2])
{
   bool all = false;
   long int kk, ii, nf, nel, nel2;

   // write header of file
   databr_head_to_vtk( ff, name, time, cycle, 1, 1, 1 );

   if( nFilds < 1 || !Flds )
   {
      all = true;
      nFilds = 51;
   }

   for( kk=0; kk<nFilds; kk++)
   {
       if( all )
         nf = kk;
       else nf= Flds[kk][0];

       databr_size_to_vtk(  nf, nel, nel2 );

       if( all )
         { ii=0; }
       else
         { ii = Flds[kk][1];
           nel = ii+1;
         }

       for( ; ii<nel; ii++ )
       {
        databr_name_to_vtk( ff, nf, ii, nel2 );
        // cycle for TNode array
        databr_element_to_vtk( ff, CNode, nf, ii );
       }
   }
}


//-----------------------End of node_format.cpp--------------------------
