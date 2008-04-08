#include <math.h>
#include <stdio.h>
#include <iomanip>
#include "io_arrays.h"

#ifndef IPMGEMPLUGIN

#include "m_unspace.h"

#else

#include "ms_unspace.h"

#endif
// !!! internal using  syp->GEX, mup->Laq, pm->tpp->G


//=========================================================================
// Output part (must be changed all )

// out information: criteria and statistics
void TUnSpace::out_QT( int Ngr  )
{
  double Ca,Ca1,sCa,sCa1, quanLapl,quanHur,quanWald,quanHom;
  float st=2.,st1=2.; // coeff. Studenta
  int i, k,ii,l,/*ca=0,*/kf, kf1;

 //  paragen( ); // put data to nPhA, PhAndx, PhAfreq ....


//==============================================
// UiDC

  ii = usp->N;
#ifndef IPMGEMPLUGIN
  if( TProfil::pm->mup->Laq )
#else
  if(mup_Laq )
#endif
      ii--;

  for( l=0; l<ii ; l++ )
  {
   Ca = sCa = Ca1 = sCa1= 0.;
   kf=0; kf1=0;
   for( i=0; i<usp->Q; i++ )
   { if(usp->sv[i]>=0)
     { Ca += usp->vMol[i*usp->N+l];
       kf++;
     }
     if( abs(usp->sv[i]) == Ngr )
     {  Ca1 += usp->vMol[i*usp->N+l];
        kf1++;
     }
   }
   if(kf)
    Ca /= kf;
   if(kf1)
    Ca1 /= kf1;

   for( i=0; i<usp->Q; i++ )
   {
     if(usp->sv[i]>=0)
      sCa += (usp->vMol[i*usp->N+l]-Ca) * (usp->vMol[i*usp->N+l]-Ca);
     if( abs(usp->sv[i]) == Ngr )
      sCa1 += (usp->vMol[i*usp->N+l]-Ca1) * (usp->vMol[i*usp->N+l]-Ca1);
   }
   if(kf>1)
     sCa /= (kf-1);
   if(kf1>1)
     sCa1 /= (kf1-1);
   if(sCa>0)
       sCa = sqrt(sCa);
   if(sCa1>0)
       sCa1 = sqrt(sCa1);

   quanLapl=quanHur=quanWald=quanHom=0.;
   for(k=0;k<usp->qQ;k++)
   { quanLapl += usp->vMol[usp->quanCx[k][0]*usp->N+l];
     quanHur +=  usp->vMol[usp->quanCx[k][1]*usp->N+l];
     quanWald += usp->vMol[usp->quanCx[k][2]*usp->N+l];
     quanHom  += usp->vMol[usp->quanCx[k][3]*usp->N+l];
   }
   usp->UnIC[l][0] = usp->vMol[l];
   usp->UnIC[l][1] = Ca;
   if(kf>0)
    usp->UnIC[l][2] = st*sCa / sqrt(double(kf));
   usp->UnIC[l][3] = Ca1;
   if(kf1>0)
    usp->UnIC[l][4] = st1*sCa1 / sqrt(double(kf1));
   usp->UnIC[l][5] = quanLapl/usp->qQ;
   usp->UnIC[l][6] = quanHur/usp->qQ;
   usp->UnIC[l][7] = quanWald/usp->qQ;
   usp->UnIC[l][8] = quanHom/usp->qQ;
   usp->UnIC[l][9] = 0.;

 }

//==============================================
// UaDC

  for( l=0; l<usp->Ls ; l++ )
  {
   Ca = sCa = Ca1 = sCa1= 0.;
   kf=0; kf1=0;

   for( i=0; i<usp->Q; i++ )
   { if(usp->sv[i]>=0)
     { Ca += pow( 10., usp->vFug[i*usp->Ls+l]);
       kf++;
     }
     if( abs(usp->sv[i]) == Ngr )
     {  Ca1 += pow( 10., usp->vFug[i*usp->Ls+l]);
        kf1++;
     }
   }
   if(kf)
    Ca /= kf;
   if(kf1)
    Ca1 /= kf1;

   for( i=0; i<usp->Q; i++ )
   {
     if(usp->sv[i]>=0)
      sCa += (pow( 10., usp->vFug[i*usp->Ls+l]-Ca)) *
                            pow( 10., (usp->vFug[i*usp->Ls+l]-Ca));
     if( abs(usp->sv[i]) == Ngr )
      sCa1 += pow( 10., (usp->vFug[i*usp->Ls+l]-Ca1)) *
                           pow( 10., (usp->vFug[i*usp->Ls+l]-Ca1));
   }
   if(kf>1)
     sCa /= (kf-1);
   if(kf1>1)
     sCa1 /= (kf1-1);
   if(sCa>0)
       sCa = sqrt(sCa);
   if(sCa1>0)
       sCa1 = sqrt(sCa1);

   quanLapl=quanHur=quanWald=quanHom=0.;
   for(k=0;k<usp->qQ;k++)
   { quanLapl += pow( 10., usp->vFug[usp->quanCx[k][0]*usp->Ls+l]);
     quanHur +=  pow( 10., usp->vFug[usp->quanCx[k][1]*usp->Ls+l]);
     quanWald += pow( 10., usp->vFug[usp->quanCx[k][2]*usp->Ls+l]);
     quanHom  += pow( 10., usp->vFug[usp->quanCx[k][3]*usp->Ls+l]);
   }
   usp->UaDC[l][0] = pow( 10., usp->vFug[l]);
   usp->UaDC[l][1] = Ca;
   if(kf>0)
    usp->UaDC[l][2] = st*sCa / sqrt(double(kf));
   usp->UaDC[l][3] = Ca1;
   if(kf1>0)
    usp->UaDC[l][4] = st1*sCa1 / sqrt(double(kf1));
   usp->UaDC[l][5] = quanLapl/usp->qQ;
   usp->UaDC[l][6] = quanHur/usp->qQ;
   usp->UaDC[l][7] = quanWald/usp->qQ;
   usp->UaDC[l][8] = quanHom/usp->qQ;
   usp->UaDC[l][9] = 0.;
 }


//==============================================
// UgDC
for( int ii=0; ii<usp->nPG; ii++ )
{
  double var1=0, var3=0, var_L=0, var_H=0, var_W=0, var_Hom=0;
  int sz1 = 0, sz3 = 0;

  usp->UgDC[ii][0] = value_nPG( ii, 0 );

  for(int q=0; q<usp->Q; q++)
  {
    if( usp->sv[q] >0 )
    { var1 += value_nPG( ii, q );
      sz1++;
    }
    if( abs(usp->sv[q]) == Ngr )
    { var3 += value_nPG( ii, q );
      sz3++;
    }
  }

   if( sz1)
     usp->UgDC[ii][1] = var1/sz1;
   else
     usp->UgDC[ii][1] = 0;

   if( sz3>0 )
     usp->UgDC[ii][3] = var3/sz3;
   else
     usp->UgDC[ii][3] = 0;

  var1 = 0;
  var3 = 0;
  for(int q=0; q<usp->Q; q++)
  {
    if( usp->sv[q] >0 )
     var1 += pow( usp->UgDC[ii][1]- value_nPG( ii, q ), 2);
    if( abs(usp->sv[q]) == Ngr )
     var3 += pow( usp->UgDC[ii][3]- value_nPG( ii, q ), 2);
  }

   if( sz1>1 )
     usp->UgDC[ii][2] = st*sqrt(var1/(sz1-1)) / sqrt(double(sz1));
   else
     usp->UgDC[ii][2] = 0;

   if( sz3>1 )
     usp->UgDC[ii][4] = st*sqrt(var3/(sz3-1)) / sqrt(double(sz3));
   else
     usp->UgDC[ii][4] = 0;

  for(int t=0; t<usp->qQ; t++)
  {
    var_L += value_nPG( ii, usp->quanCx[t][0] );
    var_H += value_nPG( ii, usp->quanCx[t][1] );
    var_W += value_nPG( ii, usp->quanCx[t][2] );
    var_Hom += value_nPG( ii, usp->quanCx[t][3] );
    // calc statistics
  }

  usp->UgDC[ii][5] = var_L/usp->qQ;
  usp->UgDC[ii][6] = var_H/usp->qQ;
  usp->UgDC[ii][7] = var_W/usp->qQ;
  usp->UgDC[ii][8] = var_Hom/usp->qQ;

  var1 = 0;
  for(int t=0; t<usp->qQ; t++)
  {
    switch( usp->Pa_Crit )
    {
     case UNSP_CRIT_LAPL_QAN:
        var1 += pow( usp->UgDC[ii][5]- value_nPG( ii, usp->quanCx[t][0] ), 2);
        break;
     case UNSP_CRIT_HUR_QAN:
        var1 += pow( usp->UgDC[ii][6]- value_nPG( ii, usp->quanCx[t][1] ), 2);
        break;
     case UNSP_CRIT_WALD_QAN:
        var1 += pow( usp->UgDC[ii][7]- value_nPG( ii, usp->quanCx[t][2] ), 2);
        break;
     case UNSP_CRIT_HOME_QAN:
        var1 += pow( usp->UgDC[ii][8]- value_nPG( ii, usp->quanCx[t][3] ), 2);
        break;
   }
  }
   if( usp->qQ >1 )
     usp->UgDC[ii][9] = st*sqrt(var1/(usp->qQ-1)) / sqrt(double(usp->qQ));
   else
     usp->UgDC[ii][9] = 0;
}
//=================================================

}

//calculation number of unspace components
double TUnSpace::value_nPG( int line, int q )
{
  int j, count=0;

  if( usp->PsGen[3]== S_ON &&  usp->NgT > 0 )
  {   if( line == count )
       return usp->vT[q];
      count++;
  }
  if( usp->PsGen[4]== S_ON &&  usp->NgP > 0 )
  {   if( line == count )
       return usp->vP[q];
      count++;
  }
  if( /*usp->PsGen[]== S_ON && */ usp->NgV > 0 )
  {   if( line == count )
       return usp->vV[q];
      count++;
  }

  if( usp->PsGen[2]== S_ON )
     for( j=0; j<usp->N; j++)
       if( usp->NgNb[j] > 0  )
       {   if( line == count )
          return usp->vB[q*usp->N+j];
         count++;
        }

  if( usp->PsGen[0]== S_ON )
     for( j=0; j<usp->L; j++)
       if( usp->NgLg[j] > 0  )
       {   if( line == count )
          return usp->vG[q*usp->L+j];
         count++;
        }

  if( usp->PsGen[1]== S_ON )
     for( j=0; j<usp->L; j++)
       if( usp->NgLs[j] > 0  )
       {   if( line == count )
          return usp->vS[q*usp->L+j];
         count++;
        }

  if( usp->PsGen[5]== S_ON )
     for( j=0; j<usp->L; j++)
       if( usp->NgLv[j] > 0  )
       {   if( line == count )
          return usp->vmV[q*usp->L+j];
         count++;
        }

  if( usp->PsGen[6]== S_ON )
     for( j=0; j<usp->Ls; j++)
       if( usp->NgGam[j] > 0  )
       {   if( line == count )
          return usp->vNidP[q*usp->Ls+j];
          count++;
        }

  return 0;
}


//  correction values G298 and intervals
void TUnSpace::AdapG()
// Iter- iteration of adaptation
{
  int Ngr;
 // Ngr - number set phases from criteria
 Ngr = sel_parg( usp->sv );


//===============================================
// UnDCA
for( int ii=0; ii<usp->nPG; ii++ )
{
  double val, sg_val, sr=0, min=0, max=0, sg_sr=0;
  int i1 = 0, sr_sz = 0;


  for(int q=0; q<usp->Q; q++)
  {
    if( usp->sv[q] == Ngr )      // filters and phase selections
    {
      val = value_nPG( ii, q );
      sr += val;
      sr_sz++;
      if(!i1)  // first
      {  min = max = val; i1++; }
      else
      {
        if( val < min ) min = val;
        if( val > max ) max = val;
      }
    }
  }

  if( sr_sz )
     sr = sr/sr_sz;
  else
     sr = 0;

  for(int q=0; q<usp->Q; q++)
  {
    if( usp->sv[q] == Ngr )
     sg_sr += pow( sr - value_nPG( ii, q ), 2);
  }

  if( sr_sz >1 )
     sg_sr = 3*sqrt(sg_sr/(sr_sz-1)) / sqrt(double(sr_sz));
  else
     sg_sr = 0;

  switch( usp->Pa_Crit )
   {
    case UNSP_CRIT_PA:         val = sr;   sg_val = sg_sr;  break;
    case UNSP_CRIT_LAPL_QAN:   val = usp->UgDC[ii][5];
                               sg_val = usp->UgDC[ii][9];
                               break;
    case UNSP_CRIT_HUR_QAN:    val = usp->UgDC[ii][6];
                               sg_val = usp->UgDC[ii][9];
                               break;
    case UNSP_CRIT_WALD_QAN:   val = usp->UgDC[ii][7];
                               sg_val = usp->UgDC[ii][9];
                               break;
    case UNSP_CRIT_HOME_QAN:   val = usp->UgDC[ii][8];
                               sg_val = usp->UgDC[ii][9];
    default:                   break;
   }
   adapt_nPG( ii, val, sg_val );

   usp->UnDCA[ii][0] = val;
   usp->UnDCA[ii][1] = sg_val;
   usp->UnDCA[ii][2] = min;
   usp->UnDCA[ii][3] = max;
   usp->UnDCA[ii][4] = sr;
   usp->UnDCA[ii][5] = sg_sr;
   usp->UnDCA[ii][6] = 0.;
   usp->UnDCA[ii][7] = 0.;
 }
}

//set up new data and inervals for  unspace components
void TUnSpace::adapt_nPG( int line, double new_val, double new_int )
{
  int j, count=0;

  if( usp->PsGen[3]== S_ON &&  usp->NgT > 0 )
  {   if( line == count )
      {
       usp->T[0] = new_val;
       usp->IntT[0] = new_int;
       return;
      }
      count++;
  }
  if( usp->PsGen[4]== S_ON &&  usp->NgP > 0 )
  {   if( line == count )
      {
       usp->P[0] = new_val;
       usp->IntP[0] = new_int;
       return;
      }
      count++;
  }
  if( /*usp->PsGen[]== S_ON && */ usp->NgV > 0 )
  {   if( line == count )
      {
       usp->V[0] = new_val;
       usp->IntV[0] = new_int;
       return;
      }
      count++;
  }

  if( usp->PsGen[2]== S_ON )
     for( j=0; j<usp->N; j++)
       if( usp->NgNb[j] > 0  )
       {  if( line == count )
          {
            usp->Bs[j][0] = new_val;
            usp->IntNb[j][0] = new_int;
            return;
          }
         count++;
     }

  if( usp->PsGen[0]== S_ON )
    for( j=0; j<usp->L; j++)
       if( usp->NgLg[j] > 0  )
       {   if( line == count )
           {
            usp->Gs[j][0] = (float)new_val;
            usp->IntLg[j][0] = new_int;
            //syu->GEX[j] = // 16/02/2007
            // usp->Gs[j][0]-float(TProfil::pm->tpp->G[j]);
            for( int jj=0; jj<pmu->L; jj++)
             if( pmu->muj[jj] == j )
             {
#ifndef IPMGEMPLUGIN
                 pmu->GEX[jj] = (usp->Gs[j][0]-
                  float(TProfil::pm->tpp->G[j]))/pmu->RT;
#else
                  pmu->GEX[jj] = (usp->Gs[j][0]-
                   float(pmu->tpp_G[j]))/pmu->RT;
#endif
               break;
             }
            return;
           }
         count++;
        }

  if( usp->PsGen[1]== S_ON )
     for( j=0; j<usp->L; j++)
       if( usp->NgLs[j] > 0  )
       {   if( line == count )
           {
             usp->Ss[j][0] = new_val;
             usp->IntLs[j][0] = new_int;
             return;
           }
         count++;
        }

  if( usp->PsGen[5]== S_ON )
     for( j=0; j<usp->L; j++)
       if( usp->NgLv[j] > 0  )
       {   if( line == count )
           {
             usp->Vs[j][0] = new_val;
             usp->IntLv[j][0] = new_int;
             return;
           }
         count++;
        }

  if( usp->PsGen[6]== S_ON )
     for( j=0; j<usp->Ls; j++)
       if( usp->NgGam[j] > 0  )
       {   if( line == count )
           {
             usp->GAMs[j][0] = new_val;
             usp->IntGam[j][0] = new_int;
             return;
           }
          count++;
        }
}

//====================================================================
extern bool _comment;

void TUnSpace::to_text_file( fstream& ff, bool with_comments )
{
  _comment = with_comments;
  
  TPrintArrays  prar(ff);

   if( _comment )
   {  ff << "# GEMIPM2K v. 2.2.3" << endl;
      ff << "# Prototype 31.03.2008" << endl;
      ff << "# Comments can be marked with # $ ;" << endl << endl;
      ff << "# Template for the UnSpace data" << endl;
      ff << "# (should be read only after the DATACH, the IPM-DAT and DATABR files)" << endl << endl;
      ff << "#Section (scalar): Controls and dimensionalities of the UnSpace operation" << endl;
   }
   if( _comment )
      ff << "# ( + on; - off ) probe type  [0]-G; [1]-S; [2]-b; [3]-T; [4]-P; [5]-V [6] pGam";
    prar.writeArray(  "PsGen", usp->PsGen, 7, 1 );
    if( _comment )
       ff << "\n# code for adaptive loops for G(298): 0-no; >=1-yes (max number of loops)" << endl;
    ff << left << setw(17) << "<Pa_Adapt> " << "\'"<<  usp->Pa_Adapt << "\'"<< endl;
    if( _comment )
       ff << "# code of function to construct the payoff matrix (A or 0; B or 1; C or 2; D or 3; E or 4; F or 5)" << endl;
    ff << left << setw(17) << "<Pa_OF> " << "\'"<<  usp->Pa_OF << "\'"<< endl;
    if( _comment )
    {
    	ff << "# (0 1 2 3 4) OSP criterion: 0- PA frequency statistics, 1- Laplace(quantile)" << endl;
        ff << "#  2-Homeniuk(quantile) 3-Laplace (single point), 4-Homeniuk (single point)" << endl;
    }
    ff << left << setw(17) << "<Pa_Crit> " << "\'"<<  usp->Pa_Crit << "\'" << endl;
    if( _comment )
       ff << "# mode of Laplace function calculation: + mean of abs.values; - abs. mean of values" << endl;
    ff << left << setw(17) << "<Pa_Zcp> " << "\'" <<  usp->Pa_Zcp << "\'" << endl;
    if( _comment )
       ff << "# code of function to construct the payoff matrix (A or 0; B or 1; C or 2; D or 3; E or 4; F or 5)" << endl;
    ff << left << setw(17) << "<PvSi> " << "\'" <<  usp->PvSi << "\'" << endl;

    if( _comment )
       ff << "\n#  Flags that control input of criteria and analysis of sampled results" << endl;
    if( _comment )
       ff << "# flag of input for the filter on phase association" << endl;
    ff << left << setw(17) << "<Pa_f_pha> " << "\'" <<  usp->Pa_f_pha << "\'" << endl;
    if( _comment )
       ff << "#  flag of input for filters on total IC molality" << endl;
    ff << left << setw(17) << "<Pa_f_mol> " << "\'" <<  usp->Pa_f_mol << "\'" << endl;
    if( _comment )
       ff << "#  flag of input for the filters on fugacity/activity " << endl;
    ff << left << setw(17) << "<Pa_f_fug> " << "\'" <<  usp->Pa_f_fug << "\'" << endl;
    if( _comment )
       ff << "# flag of input for the filters on mole fraction" << endl;
    ff << left << setw(17) << "<Pa_f_mfr> " << "\'" <<  usp->Pa_f_mfr << "\'" << endl;
    if( _comment )
       ff << "# flag of input for the filter on pH" << endl;
    ff << left << setw(17) << "<Pa_f_pH> " << "\'" <<  usp->Pa_f_pH << "\'" << endl;
    if( _comment )
       ff << "# flag of input for the filter on Eh" << endl;
    ff << left << setw(17) << "<Pa_f_Eh> " << "\'" <<  usp->Pa_f_Eh << "\'" << endl;
    if( _comment )
       ff << "# flag of input for the filter on Independent Components" << endl;
    ff << left << setw(17) << "<Pa_f_IC> " << "\'" <<  usp->Pa_f_IC << "\'" << endl;

    if( _comment )
       ff << "\n#  Allocation flags" << endl;
    if( _comment )
       ff << "# allocation flag for payoff matrix (if Q <= 1001)" << endl;
    ff << left << setw(17) << "<PvPOM> " << "\'" <<  usp->PvPOM << "\'" << endl;
    if( _comment )
       ff << "#  allocation flag for payoff matrix row (if Q > 1001)" << endl;
    ff << left << setw(17) << "<PvPOR> " << "\'" <<  usp->PvPOR << "\'" << endl;
   
    if( _comment )
       ff << "\n#  Dimensionalities related to the UnSpace problem (Project/System indexation)" << endl;
    if( _comment )
       ff << "# input number of sample GEM calculations to be generated" << endl;
    ff << left << setw(7) << "<Q> " <<  usp->Q << endl;
    if( _comment )
       ff << "#  quantile level (default 0.05)" << endl;
    ff << left << setw(7) << "<quan_lev> " <<  usp->quan_lev << endl;
    if( _comment )
       ff << "# Total number of uncertain input parameters (for all unspace groups)" << endl;
    ff << left << setw(7) << "<nPG> " <<  usp->nPG << endl;
    if( _comment )
       ff << "#  number of UnSpace parameter groups for sampling by Belov's grid" << endl;
    ff << left << setw(7) << "<nGB> " <<  usp->nGB << endl;
    if( _comment )
       ff << "# number of UnSpace groups used in Monte Carlo with normal distribution" << endl;
    ff << left << setw(7) << "<nGN> " <<  usp->nGN << endl;
    if( _comment )
       ff << "#  number of UnSpace groups used in Monte Carlo with uniform distribution" << endl;
    ff << left << setw(7) << "<nGR> " <<  usp->nGR << endl;

    if( _comment )
       ff << "# Number of independent components   (from pmu->N)" << endl;
    ff << left << setw(7) << "<N> " <<  /*usp*/pmu->N << endl;
    if( _comment )
       ff << "#  Number of dependent components     (from pmu->L)" << endl;
    ff << left << setw(7) << "<L> " <<  /*usp*/pmu->L << endl;
    if( _comment )
       ff << "# Total number of DC in multi-comp.phases (from pmu->Ls)" << endl;
    ff << left << setw(7) << "<Ls> " <<  /*usp*/pmu->Ls << endl;
    if( _comment )
       ff << "# Number of phases                   (from pmu->Fi)" << endl;
    ff << left << setw(7) << "<Fi> " <<  /*usp*/pmu->FI << endl;
 
    ff<< "\n<END_DIM>\n";

  // dynamic arrays - must follow static data
    if( _comment )
    {   ff << "\n## Task configuration section ";
        ff << "\n# List of UnSpace group name";
    }
   prar.writeArray(  "SGp", usp->SGp[0], usp->nG, MAXPHNAME);
    if( _comment )
      ff << "\n# generation codes for unspace groups: 0- Belov; 1-uniform; 2-normale";
   prar.writeArray(  "PbD", usp->PbD, usp->nG);
   if( _comment )
     ff << "\n# pseudo-random numbers for Belov algorithm";
   prar.writeArray(  "OVB", usp->OVB, usp->nGB);
   if( _comment )
     ff << "\n# pseudo-random numbers by Monte Carlo with uniform distribution";
   prar.writeArray(  "OVR", usp->OVR, usp->nGR);
   if( _comment )
      ff << "\n#  pseudo-random numbers by Monte Carlo with normal distribution";
   prar.writeArray(  "OVN", usp->OVN, usp->nGN);

   if( _comment )
      ff << "\n\n## Input filters for sampled solution variants ";
   if( _comment )
      ff << "\n# pH lower filter limit (default 0)" << endl;
   ff << left << setw(7) << "<pH_lo> " <<  usp->pH_lo << endl;
   if( _comment )
      ff << "#  pH upper filter limit (default 14)" << endl;
   ff << left << setw(7) << "<pH_up> " <<  usp->pH_up << endl;
   if( _comment )
      ff << "# Eh lower filter limit (default -1)" << endl;
   ff << left << setw(7) << "<Eh_lo> " <<  usp->Eh_lo << endl;
   if( _comment )
      ff << "# Eh upper filter limit (default 1)" << endl;
   ff << left << setw(7) << "<Eh_up> " <<  usp->Eh_up << endl;
   if( _comment )
      ff << "# Ionic strength lower limit (default 0)" << endl;
   ff << left << setw(7) << "<IC_lo> " <<  usp->IC_lo << endl;
   if( _comment )
      ff << "# Ionic strength lower limit (default 3))" << endl;
   ff << left << setw(7) << "<IC_up> " <<  usp->IC_up << endl;
   if( _comment )
     ff << "# total IC molality filter - lower limits (default 0/-20) (units defined in PsUnFltI flag)";
   prar.writeArray(  "m_t_lo", usp->m_t_lo, pmu->N/*usp->N*/, pmu->mui);
   if( _comment )
    ff << "\n# total IC molality filter - upper limits (default 10/+1) (units defined in PsUnFltI flag)";
   prar.writeArray(  "m_t_up", usp->m_t_up, pmu->N/*usp->N*/, pmu->mui);
   if( _comment )
    ff << "\n# fugacity/activity filter - lower limits (default 0/-20) (units defined in PsUnFltD flag)";
   prar.writeArray(  "fug_lo", usp->fug_lo, pmu->Ls/*usp->Ls*/, pmu->muj);
   if( _comment )
     ff << "\n#  fugacity/activity filter - upper limits (default 10/+1) (units defined in PsUnFltD flag)";
   prar.writeArray(  "fug_up", usp->fug_up, pmu->Ls/*usp->Ls*/, pmu->muj);
   if( _comment )
     ff << "\n#  input phase indexes of known equilibrium phase association for a primary filter (Pa_f_pha)";
   prar.writeArray(  "f_PhA", usp->f_PhA, pmu->N/*usp->N*/, pmu->mui);
   
   if( _comment )
      ff << "\n\n## Input intervals for filtering sampled GEM solution variants ";
   if( _comment )
      ff << "\n# Unspace group index for uncertain temperature T; 0- T is not uncertain" << endl;
   ff << left << setw(7) << "<NgT> " <<  usp->NgT << endl;
   if( _comment )
     ff << "# Adapted/Initial Deterministic value of T (Celsius)";
   prar.writeArray(  "T", usp->T, 2);
   if( _comment )
    ff << "\n# Adapted/Initial 0 or uncertainty half-interval for T (C)";
   prar.writeArray(  "IntT", usp->IntT, 2);
 
   if( _comment )
      ff << "\n# Unspace group index for uncertain pressure P; 0- P not uncertain" << endl;
   ff << left << setw(7) << "<NgP> " <<  usp->NgP << endl;
   if( _comment )
     ff << "# Adapted/Initial Deterministic value of P (bar)";
   prar.writeArray(  "P", usp->P, 2);
   if( _comment )
    ff << "\n# Adapted/Initial 0 or uncertainty half-interval for P (bar)";
   prar.writeArray(  "IntP", usp->IntP, 2);

   if( _comment )
      ff << "\n# Unspace group index for uncertain volume constraint V, 0 - V not uncertain (reserved)" << endl;
   ff << left << setw(7) << "<NgV> " <<  usp->NgV << endl;
   if( _comment )
     ff << "# Adapted/Initial Deterministic value of V (in cm3? L?)";
   prar.writeArray(  "V", usp->V, 2);
   if( _comment )
    ff << "\n# Adapted/Initial 0 or uncertainty half-interval for V (in cm3? L?)";
   prar.writeArray(  "IntV", usp->IntV, 2);

   if( _comment )
     ff << "\n# 0-deterministic DC; value >0 is index of a group of varying G to which this DC is assigned";
   prar.writeArray(  "NgLg", usp->NgLg, pmu->L/*usp->L*/, pmu->muj);
   if( _comment )
     ff << "\n# adapted/initial  copy of deterministic values of G298 for DCs";
   prar.writeArray(  "Gs", usp->Gs[0], pmu->L/*usp->L*/, pmu->muj, 2);
   if( _comment )
    ff << "\n# adapted/initial  uncertainty half-intervals for G variation (flag PsUnInt defines units)";
   prar.writeArray(  "IntLg", usp->IntLg[0], pmu->L/*usp->L*/, pmu->muj, 2);

   if( _comment )
     ff << "\n# 0-deterministic DC; value >0 is index of a group varying gamma params (reserved)";
   prar.writeArray(  "NgGam", usp->NgGam, pmu->Ls/*usp->Ls*/, pmu->muj);
   if( _comment )
     ff << "\n# adapted/initial copy of deterministic Par Gamma (reserved)";
   prar.writeArray(  "GAMs", usp->GAMs[0], pmu->Ls/*usp->Ls*/, pmu->muj, 2);
   if( _comment )
    ff << "\n# adapted/initial uncertainty half-intervals for gamma params (reserved)";
   prar.writeArray(  "IntGam", usp->IntGam[0], pmu->Ls/*usp->Ls*/, pmu->muj, 2);

   if( _comment )
      ff << "\n# 0-deterministic DC; value >0 is index of a group varying S";
    prar.writeArray(  "NgLs", usp->NgLs, pmu->L/*usp->L*/, pmu->muj);
    if( _comment )
      ff << "\n# adapted/initial  copy of deterministic values of S298 for DCs (reserved)";
    prar.writeArray(  "Ss", usp->Ss[0], pmu->L/*usp->L*/, pmu->muj, 2);
    if( _comment )
     ff << "\n#  adapted/initial uncertainty half-intervals for S (reserved)";
    prar.writeArray(  "IntLs", usp->IntLs[0], pmu->L/*usp->L*/, pmu->muj, 2);
   
    if( _comment )
       ff << "\n#  0-deterministic DC; value >0 is index of a group varying V";
    prar.writeArray(  "NgLv", usp->NgLv, pmu->L/*usp->L*/, pmu->muj);
    if( _comment )
       ff << "\n# adapted/initial  copy of deterministic values of V298 for DCs";
    prar.writeArray(  "Vs", usp->Vs[0], pmu->L/*usp->L*/, pmu->muj, 2);
    if( _comment )
      ff << "\n#  adapted/initial uncertainty half-intervals for V variation (flag PsUnInt defines units)";
    prar.writeArray(  "IntLv", usp->IntLv[0], pmu->L/*usp->L*/, pmu->muj, 2);
    
    if( _comment )
       ff << "\n#  0-deterministic DC; value >0 is index of a group varying b to which this IC is assigned";
    prar.writeArray(  "NgNb", usp->NgNb, pmu->N/*usp->N*/, pmu->mui);
    if( _comment )
       ff << "\n# adapted/initial  copy of deterministic values of bulk composition vector b";
    prar.writeArray(  "Bs", usp->Bs[0], pmu->N/*usp->N*/, pmu->mui, 2);
    if( _comment )
      ff << "\n#  adapted/initial uncertainty half-intervals for b variation (flag PsUnInB defines units)";
    prar.writeArray(  "IntNb", usp->IntNb[0], pmu->N/*usp->N*/, pmu->mui, 2);
    
    if( _comment )
      ff << "\n## Names of columns for output array UnIC ";
    prar.writeArray(  "UnICn", usp->UnICn[0], UNSP_SIZE1, NAME_SIZE);
    if( _comment )
      ff << "\n## Names of columns for output array UgDC  ";
    prar.writeArray(  "UgDCn", usp->UgDCn[0], UNSP_SIZE1, NAME_SIZE);
    if( _comment )
      ff << "\n## Names of columns for output array  UaDc ";
    prar.writeArray(  "UaDCn", usp->UaDCn[0], UNSP_SIZE1, NAME_SIZE);
    if( _comment )
      ff << "\n## Names of columns for output array UnDCA ";
    prar.writeArray(  "UnDCAn", usp->UnICn[0], UNSP_SIZE2, NAME_SIZE);
    if( _comment )
      ff << "\n## List of insertain input parameters names ";
    prar.writeArray(  "ParNames", usp->ParNames[0], usp->nPG, PARNAME_SIZE);

//    if( _comment )
//     ff << "\n# table of normalised coordinates of points in uncertainty space";
//    prar.writeArray(  "ncp", usp->ncp, usp->Q*usp->nG );
    
   if( _comment )
         ff << "\n# End of file"<< endl;
}

void TUnSpace::result_to_text_file( fstream& ff, bool with_comments )
{
  _comment = with_comments;
  
  TPrintArrays  prar(ff);

  if( _comment )
   ff << "\n# table of normalised coordinates of points in uncertainty space";
  prar.writeArray(  "ncp", usp->ncp, usp->Q*usp->nG );
  if( _comment )
    ff << "\n# index of optimal sample point";
   prar.writeArray(  "vY", usp->vY, usp->Q*usp->L, usp->L );
   if( _comment )
     ff << "\n# G0 values used in sample input data variants (indexes q j)";
   prar.writeArray(  "vG", usp->vG, usp->Q*usp->L,  usp->L );
   if( _comment )
     ff << "\n# m_t values of total IC molality from sample GEM solution variants (indexes t i)";
  prar.writeArray(  "vMol", usp->vMol, usp->Q*usp->N, usp->N );
  if( _comment )
     ff << "\n#   lnGam values for DC from sample GEM solution variants (indexes t j)";
   prar.writeArray(  "vGam", usp->vGam, usp->Q*usp->L, usp->L );
   if( _comment )
      ff << "\n# XF (phase mole amounts) from sample GEM solution variants (indexes t k)";
    prar.writeArray(  "vYF", usp->vYF, usp->Q*usp->Fi, usp->Fi );
  if( _comment )
       ff << "\n#  u values of dual chemical potentials from sample GEM solution variants";
   prar.writeArray(  "vU", usp->vU, usp->Q*usp->N,usp->N );
   if( _comment )
      ff << "\n# lga values of log10 fugacity/activity from sample GEM solution variants";
    prar.writeArray(  "vFug", usp->vFug, usp->Q*usp->Ls );
  if( _comment )
       ff << "\n#  vector of pH, Eh, IC  values from sample GEM solution variants";
   prar.writeArray(  "vpH", usp->vpH[0], usp->Q*3, 3 );
  
   if( _comment )
        ff << "\n#   list of phase assemblage titles (made of Phase names";
   prar.writeArray(  "PhAlst", usp->PhAlst[0], usp->nPhA, 80 );
   prar.writeArray(  "PhNum", usp->PhNum, usp->nPhA );
   prar.writeArray(  "PhAfreq", usp->PhAfreq, usp->nPhA );
   prar.writeArray(  "sv", usp->sv, usp->Q );

   if( _comment )
    ff << "\n\n# index of optimal sample point";
   prar.writeArray(  "Lapl", &usp->Lapl, 4 );
   if( _comment )
    ff << "\n# number of sample points with the same phase assemblage as that selected by";
   prar.writeArray(  "nl", &usp->nl, 4 );
   if( _comment )
    ff << "\n# single criterion value (from analysing payoff matrix)";
   prar.writeArray(  "CrL", &usp->CrL, 4 );
  
   if( _comment )
    ff << "\n# indices of sample GEM variants taken into quantile Laplace, Hurtvitz, Wald, Homenyuk (columns";
   prar.writeArray(  "quanCv", usp->quanCv[0], usp->qQ*4, 4 );
   if( _comment )
    ff << "\n# Values taken into quantile Laplace. Hurtvitz, Wald, Homenyuk (columns )";
   prar.writeArray(  "quanCx", usp->quanCx[0],  usp->qQ*4, 4);

   if( usp->PvPOM == S_ON )
  {  if( _comment )
     ff << "\n# payoff matrix";
    prar.writeArray(  "POM", usp->POM, usp->Q*usp->Q, usp->Q );
  }
  if( usp->PvPOR == S_ON )
    {  if( _comment )
       ff << "\n# payoff matrix row";
      prar.writeArray(  "POR", usp->POR, usp->Q );
    }

  if( _comment )
     ff << "\n# statistics over independent components";
  prar.writeArray(  "UnIC", usp->UnIC[0],  usp->N*UNSP_SIZE1, UNSP_SIZE1 );
  if( _comment )
     ff << "\n# V statistics set 1 over dependent components";
  prar.writeArray(  "UgDC", usp->UgDC[0],  usp->nPG*UNSP_SIZE1, UNSP_SIZE1);
  if( _comment )
     ff << "\n# statistics set 2 over dependent components";
  prar.writeArray(  "UaDC", usp->UaDC[0],  usp->Ls*UNSP_SIZE1, UNSP_SIZE1 );
  if( _comment )
     ff << "\n# statistics set 1 over dependent components";
  prar.writeArray(  "UnDCA", usp->UnDCA[0],  usp->nPG*UNSP_SIZE2, UNSP_SIZE2);


//  if( _comment )
//     ff << "\n# matrix A - test";
//  prar.writeArray(  "A", usp->A,  usp->L*usp->N, usp->L);

  if( _comment )
        ff << "\n# End of file"<< endl;
}

// end of m_unsp_out.cpp
