//-------------------------------------------------------------------
// $Id$
//
/// \file tsolmod_multi_file.cpp
/// Implementation of writing IPM trace file and allocation model arrays
//
// Copyright (c) 2023 S.Dmytriyeva, D.Kulik
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

#include "tsolmod_multi.h"
#include "v_service.h"
#include "io_template.h"
#include "io_keyvalue.h"
#include "gdatastream.h"
#include "jsonconfig.h"

void TSolModMulti::getLsModsum( long int& LsModSum, long int& LsIPxSum )
{  LsModSum = 0;
    LsIPxSum = 0;
    for(long int i=0; i<pm.FIs; i++)
    {
        LsModSum += (pm.LsMod[i*3]*pm.LsMod[i*3+2]);
        LsIPxSum += (pm.LsMod[i*3]*pm.LsMod[i*3+1]);
    }
}

void TSolModMulti::getLsMdcsum( long int& LsMdcSum,long int& LsMsnSum,long int& LsSitSum )
{  LsMdcSum = 0;
    LsMsnSum = 0;
    LsSitSum = 0;

    for(long int i=0; i<pm.FIs; i++)
    {
        LsMdcSum += (pm.LsMdc[i*3]*pm.L1[i]);
        LsMsnSum += (pm.LsMdc[i*3+1]*pm.LsMdc[i*3+2]*pm.L1[i]);
        LsSitSum += (pm.LsMdc[i*3+1]*pm.LsMdc[i*3+2]);
    }
}

// dimensions from LsPhl array
void TSolModMulti::getLsPhlsum( long int& PhLinSum,long int& lPhcSum )
{  PhLinSum = 0;
    lPhcSum = 0;

    for(long int i=0; i<pm.FI; i++)
    {
        PhLinSum += (pm.LsPhl[i*2]);
        lPhcSum += (/*pm.LsPhl[i*2]**/pm.LsPhl[i*2+1]);

    }
}

// dimensions from LsMdc2 array
void TSolModMulti::getLsMdc2sum( long int& DQFcSum,long int& rcpcSum )
{  DQFcSum = 0;
    rcpcSum = 0;

    for(long int i=0; i<pm.FIs; i++)
    {
        DQFcSum += (pm.LsMdc2[i*3]*pm.L1[i]);
        //       rcpcSum += (pm.LsMdc2[i*3+1]*pm.L1[i]);
    }
}

//---------------------------------------------------------//
/// Writing structure MULTI ( free format file  )
void TSolModMulti::to_text_file( const char *path, bool append )
{
    std::ios::openmode mod = std::ios::out;
    if( append )
        mod = std::ios::out|std::ios::app;
    std::fstream ff(GemsSettings::with_directory(path), mod );
    ErrorIf( !ff.good() , path, "Fileopen error");

    io_formats::KeyValueWrite out_format( ff );
    out_format.put_head( "", "ipm");
    io_formats::TPrintArrays<io_formats::KeyValueWrite>  prar( 0, {}, out_format );

    if( append )
        prar.writeComment( true,"\nNext record" );
    prar.writeComment( true, char_array_to_string(pm.stkey, EQ_RKLEN)+"\n" );
    //  TProfil::pm->pa.p.write(ff);

    //prar.writeArray( "Short_PARAM",  &base_param()->PC, 10L );
    //prar.writeArray( "Double_PARAM",  &base_param()->DG, 28L );
    prar.writeArray( "Short_Const",  &pm.N, 39L );
    prar.writeArray(  "Double_Const",  &pm.TC, 53, 20 );
    // prar.writeArray(  "Add_Double_Const",  &pm.XwMinM, 12, 20 );
    prar.writeArray(  "EpsW", pm.epsW, 5);
    prar.writeArray(  "EpsWg", pm.epsWg, 5);
    prar.writeArray(  "DenW", pm.denW, 5);
    prar.writeArray(  "DenWg", pm.denWg, 5);
    prar.writeComment( true, std::string("Error Code ")+ pm.errorCode);
    prar.writeComment( true, std::string("Error Message") + pm.errorBuf);

    //dynamic values

    // Part 1
    prar.writeArray(  "L1", pm.L1,  pm.FI);
    prar.writeArray(  "Vol", pm.Vol,  pm.L);
    prar.writeArray(  "Pparc", pm.Pparc,  pm.L);
    prar.writeArray(  "MM", pm.MM,  pm.L);
    prar.writeArray(  "Awt", pm.Awt, pm.N);
    prar.writeArray(  "A", pm.A,  pm.N*pm.L);
    prar.writeArray(  "G", pm.G,  pm.L);
    prar.writeArray(  "G0", pm.G0,  pm.L);
    prar.writeArray(  "lnGam", pm.lnGam,  pm.L);
    prar.writeArray(  "lnGmo", pm.lnGmo,  pm.L);
    prar.writeArray(  "B", pm.B,  pm.N);
    prar.writeArray(  "X", pm.X,  pm.L);
    prar.writeArray(  "YOF", pm.YOF,  pm.FI);
    prar.writeArray(  "lnGmM", pm.lnGmM,  pm.L);
    prar.writeArray(  "fDQF", pm.fDQF,  pm.L);
    prar.writeArray(  "FVOL", pm.FVOL,  pm.FI);
    prar.writeArray(  "FWGT", pm.FWGT,  pm.FI);

    if( pm.L > 0 )
    {
        prar.writeArray(  "Wx", pm.Wx,  pm.L);
        prar.writeArray(  "Gamma", pm.Gamma,  pm.L);
        prar.writeArray(  "lnGmf", pm.lnGmf,  pm.L);
    }

    // Part 2  not always required arrays
    if( pm.LO > 1 )
    {
        prar.writeArray(  "Y_m", pm.Y_m,  pm.L);
    }
    if( pm.E )
    {
        prar.writeArray(  "EZ", pm.EZ,  pm.L);
    }

    // Part 3  new Phase definition
    if( pm.FIs > 0 && pm.Ls > 0 )
    {
        prar.writeArray(  "sMod", &pm.sMod[0][0], pm.FIs,8L);
        prar.writeArray(  "LsMod", pm.LsMod, pm.FIs*3);
        long int LsModSum;
        long int LsIPxSum;
        getLsModsum( LsModSum, LsIPxSum );
        prar.writeArray(  "IPxPH", pm.IPx,  LsIPxSum);
        prar.writeArray(  "PMc", pm.PMc,  LsModSum);
        long int LsMdcSum;
        long int LsMsnSum;
        long int LsSitSum;
        prar.writeArray(  "LsMdc", pm.LsMdc, pm.FIs*3);
        getLsMdcsum( LsMdcSum,LsMsnSum, LsSitSum );
        prar.writeArray(  "DMc", pm.DMc,  LsMdcSum);
        prar.writeArray(  "MoiSN", pm.MoiSN,  LsMsnSum);
        prar.writeArray(  "SitFr", pm.SitFr,  LsSitSum);
        long int DQFcSum, rcpcSum;
        getLsMdc2sum( DQFcSum, rcpcSum );
        prar.writeArray(  "LsMdc2", pm.LsMdc2, pm.FIs*3);
        prar.writeArray(  "DQFc", pm.DQFc,  DQFcSum);
        //      prar.writeArray(  "rcpc", pm.rcpc,  rcpcSum);
        long int PhLinSum, lPhcSum;
        getLsPhlsum( PhLinSum,lPhcSum );
        prar.writeArray(  "LsPhl", pm.LsPhl, pm.FI*2);
        prar.writeArray(  "PhLin", &pm.PhLin[0][0], PhLinSum*2);
        prar.writeArray(  "lPhc", pm.lPhc,  lPhcSum);

        prar.writeArray(  "lnDQFt", pm.lnDQFt, pm.Ls);
        prar.writeArray(  "lnRcpt", pm.lnRcpt, pm.Ls);
        prar.writeArray(  "lnExet", pm.lnExet, pm.Ls);
        prar.writeArray(  "lnCnft", pm.lnCnft, pm.Ls);
    }

    if(pm.H0)
        prar.writeArray("H0",pm.H0, pm.L);
    if(pm.A0)
        prar.writeArray("A0",pm.A0, pm.L);
    if(pm.U0)
        prar.writeArray("U0",pm.U0, pm.L);
    if(pm.S0)
        prar.writeArray("S0",pm.S0, pm.L);
    if(pm.Cp0)
        prar.writeArray("Cp0",pm.Cp0, pm.L);

    prar.writeArray(  "VPh", &pm.VPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "GPh", &pm.GPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "HPh", &pm.HPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "SPh", &pm.SPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "CPh", &pm.CPh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "APh", &pm.APh[0][0], pm.FIs*MIXPHPROPS);
    prar.writeArray(  "UPh", &pm.UPh[0][0], pm.FIs*MIXPHPROPS);
}

void TSolModMulti::alloc_IPx( long int LsIPxSum )
{
    if( pm.IPx ) delete[] pm.IPx;
    pm.IPx = new long int[ LsIPxSum];
}

void TSolModMulti::alloc_PMc( long int LsModSum )
{
    if( pm.PMc ) delete[] pm.PMc;
    pm.PMc = new double[LsModSum];
}

void TSolModMulti::alloc_DMc( long int LsMdcSum )
{
    if( pm.DMc ) delete[] pm.DMc;
    pm.DMc = new double[LsMdcSum];
}

void TSolModMulti::alloc_MoiSN( long int LsMsnSum )
{
    if(pm.MoiSN) delete[] pm.MoiSN;
    pm.MoiSN = new double[LsMsnSum];
}

void TSolModMulti::alloc_SitFr( long int LsSitSum )
{
    if(pm.SitFr) delete[] pm.SitFr;
    pm.SitFr = new double[LsSitSum];
}

void TSolModMulti::alloc_DQFc( long int DQFcSum )
{
    if(pm.DQFc) delete[] pm.DQFc;
    pm.DQFc = new double[DQFcSum];
}

void TSolModMulti::alloc_PhLin( long int PhLinSum )
{
    if(pm.PhLin) delete[] pm.PhLin;
    pm.PhLin = new long int[PhLinSum][2];
}

void TSolModMulti::alloc_lPhc( long int lPhcSum )
{
    if(pm.lPhc) delete[] pm.lPhc;
    pm.lPhc = new double[lPhcSum];
}

//--------------------- end of tsolmod_multi_file.cpp ---------------------------


