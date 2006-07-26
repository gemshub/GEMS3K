#include <iomanip>
#include  <iostream>

#include "m_param.h"
#include "node.h"

#ifdef IPMGEMPLUGIN
  istream& f_getline(istream& is, gstring& str, char delim);
#endif

bool _comment = true;

// skip  ' ',  '\n', '\t' and comments (from '#' to end of line)
void  skipSpace( fstream& ff )
{
  char input;
  ff.get( input );
  while( input == '#' || input == ' ' ||
        input == '\n' || input == '\t')
 {
   if( input == '#' )
    do{
         ff.get( input );
      }while( input != '\n' && input != '\0');
   if( input == '\0' )
     return;
   ff.get( input );
  }
 ff.putback(input);
}

void inArray( fstream& ff, char *name, short* arr, int size )
{
 char buf[200];
 skipSpace( ff );
 ff >> buf;
 if( memcmp( name, buf+1, strlen(name)))
    Error( buf, "DataBR text read 01: Invalid name of (int*2) array");
 for( int ii=0; ii<size; ii++  )
 {
   skipSpace( ff );
   ff >> arr[ii];
 }
}

void inArray( fstream& ff, char *name, float* arr, int size )
{
 char buf[200];
 skipSpace( ff );
 ff >> buf;
 if( memcmp( name, buf+1, strlen(name)))
    Error( buf, "DataBR text read 02: Invalid name of float (real*4) array");
 for( int ii=0; ii<size; ii++  )
 {
     skipSpace( ff );
     ff >> arr[ii];
 }
}

void inArray( fstream& ff, char *name, double* arr, int size )
{
 char buf[200];
 skipSpace( ff );
 ff >> buf;
 if( memcmp( name, buf+1, strlen(name)))
    Error( buf, "DataBR text read 03: Invalid name of double (real*8) array");
 for( int ii=0; ii<size; ii++  )
 {
     skipSpace( ff );
     ff >> arr[ii];
 }
}

void inArray( fstream& ff, char *name, char* arr,
                              int size, int el_size )
{
 char ch;
 char buf[200];
 skipSpace( ff );
 ff >> buf;
 if( memcmp( name, buf+1, strlen(name)))
    Error( buf, "DataBR text read 04: Invalid name of char[el_size] array");

 for( int ii=0; ii<size; ii++  )
 {
   skipSpace( ff );
   ff.get(ch);
   while( ff.good() && ch != '\'' )
       ff.get(ch);
   ff.getline( buf, el_size+1, '\'');
   strncpy( arr +(ii*el_size), buf, el_size );
 }

}

//---------------------------------------------------------//
// for test out data
void outArray( fstream& ff, char *name, char* arr,
                              int size, int arr_siz )
{
 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == 10)
    { jj=0;  ff << endl;}
    gstring str = gstring( arr +(ii*arr_siz), 0, arr_siz );
    str.strip();
    ff  << "\'" << str.c_str() << "\'" << " ";
 }
}

void outArray( fstream& ff, char *name, short* arr,
                 int size, int l_size  )
{
  int sz = 10;
  if( l_size > 0 )
        sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    ff << arr[ii] << " ";
 }
}

void outArray( fstream& ff, char *name,  float* arr,
            int size, int l_size )
{
 int sz = 10;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    ff << setprecision(10) << scientific << arr[ii] << " ";
 }
}

void outArray( fstream& ff, char *name,  double* arr,
            int size, int l_size )
{
 int sz = 10;
 if( l_size > 0 )
       sz = l_size;

 ff << endl << "<" << name << ">" << endl;
 for( int ii=0, jj=0; ii<size; ii++, jj++  )
 {
    if(jj == sz)
    { jj=0;  ff << endl;}
    ff << setprecision(18) << scientific << arr[ii] << " ";
 }
}

//===================================================================

void TMulti::to_text_file_gemipm( const char *path )
{
  SPP_SETTING *pa = &TProfil::pm->pa;

   //static values
   char PAalp;
   char PSigm;
   float EpsW;
   float RoW;

#ifndef IPMGEMPLUGIN
   PAalp = syp->PAalp;
   PSigm = syp->PSigm;
   EpsW = TProfil::pm->tpp->EpsW;
   RoW = TProfil::pm->tpp->RoW;
#else
   PAalp = PAalp_;
   PSigm = PSigm_;
   EpsW = EpsW_;
   RoW = RoW_;
#endif
  fstream ff( path, ios::out );
  ErrorIf( !ff.good() , path, "Fileopen error");

if( _comment )
{   ff << "# GEMIPM2K v. 0.725" << endl;
   ff << "# Prototype 12.07.2006" << endl;
   ff << "# Comments marked with #" << endl << endl;
   ff << "# Template for the ipm-dat text input file for the internal MULTI data" << endl;
   ff << "# (should be read after the DATACH file and before DATABR files)" << endl << endl;
   ff << "# ID key of the initial chemical system definition" << endl;
}
  ff << "\"" << pmp->stkey << "\"" << endl << endl;
// static data
  if( _comment )
  { ff << "## (1) Controls of the numerical behavior of the GEM IPM algorithm" << endl;
    ff << "# Minimum amount of independent component in bulk composition (except charge Zz), moles" << endl;
  }
   ff << left << setw(12) << "<pa_DB> " <<  right << setw(8) << pa->p.DB << endl;
   if( _comment )
     ff << "\n# Maximum allowed mass balance residual (moles) for major independent components" << endl;
   ff << left << setw(12) << "<pa_DHB> " << right << setw(8) <<  pa->p.DHB << endl;
   if( _comment )
     ff << "# Precision criterion of the simplex() procedure to obtain automatic initial approximation" << endl;
   ff << left << setw(12) << "<pa_EPS> " <<  right << setw(8) << pa->p.EPS << endl;
   if( _comment )
     ff << "\n# IPM convergence threshold for the Dikin criterion (1e-6 < DK < 1e-4)" << endl;
   ff << left << setw(12) << "<pa_DK> " <<  right << setw(8) << pa->p.DK << endl;
   if( _comment )
     ff << "\n# Threshold for the application of the Karpov phase stability criterion f_alpha" << endl;
   ff << left << setw(12) << "<pa_DF> " <<  right << setw(8) << pa->p.DF << endl;
   if( _comment )
     ff << "\n# Maximum allowed number of iterations in the EnterFeasibleDomain() procedure" << endl;
   ff << left << setw(12) << "<pa_DP> " << right << setw(8) << pa->p.DP << endl;
   if( _comment )
     ff << "\n# Maximum allowed number of iterations in the MainIPM_Descent() procedure" << endl;
   ff << left << setw(12) << "<pa_IIM> " << right << setw(8) <<  pa->p.IIM << endl;
   if( _comment )
     ff << "\n# Control on calling built-in Debye-Hueckel() and other models for aqueous activity coefficients" << endl;
   ff << left << setw(12) << "<pa_PD> " <<  right << setw(8) << pa->p.PD << endl;
   if( _comment )
   { ff << "\n# Positive number: control of calling IPM_gamma() on iterations of GEM EFD and IPM (default 3)" << endl;
     ff << "# Negative number: the number of additional EFD-IPM loops to improve the GEM final solution" << endl;
   }
   ff << left << setw(12) << "<pa_PRD> " <<  right << setw(8) << pa->p.PRD << endl;
   if( _comment )
     ff << "\n# Smoothing parameters controlling convergence in highly non-ideal systems " << endl;
   ff << left << setw(12) << "<pa_AG> " <<  right << setw(8) << pa->p.AG << endl;
   ff << left << setw(12) << "<pa_DGC> " <<  right << setw(8) << pa->p.DGC << endl;
   if( _comment )
     ff << "\n# Flag for using initial activity coefficients for simplex() initial approximation (1-enable, 0-disable)" << endl;
   ff << left << setw(12) << "<pa_PSM> " <<  right << setw(8) << pa->p.PSM << endl;
   if( _comment )
     ff << "# Activity coefficient values for simplex(): GAR for major and GAH for minor components" << endl;
   ff << left << setw(12) << "<pa_GAR> " <<  right << setw(8) << pa->p.GAR << endl;
   ff << left << setw(12) << "<pa_GAH> " <<  right << setw(8) << pa->p.GAH << endl;
   if( _comment )
     ff << "\n# Cutoff threshold for the amount of phase for phase elimination " << endl;
   ff << left << setw(12) << "<pa_DS> " << right << setw(8) <<  pa->p.DS << endl;
   if( _comment )
   {  ff << "\n# Cutoffs for elimination of: Xw water solvent; Sc - solid sorbent; Dc - solution or surface species; " << endl;
      ff << "# Ph - non-electrolyte solution phases" << endl;
   }
   ff << left << setw(12) << "<pa_XwMin> " <<  right << setw(8) << pa->p.XwMin << endl;
   ff << left << setw(12) << "<pa_ScMin> " <<  right << setw(8) << pa->p.ScMin << endl;
   ff << left << setw(12) << "<pa_DcMin> " <<  right << setw(8) << pa->p.DcMin << endl;
   ff << left << setw(12) << "<pa_PhMin> " <<  right << setw(8) << pa->p.PhMin << endl;
   if( _comment )
     ff << "\n# Minimal ionic strength below which the activity coefficients for aqueous species are set to 1" << endl;
   ff << left << setw(12) << "<pa_ICmin> " <<  right << setw(8) << pa->p.ICmin << endl;
   if( _comment )
     ff << "\n# Mode of Selekt2() procedure operation" << endl;
   ff << left << setw(12) << "<pa_PC> " <<  right << setw(8) << pa->p.PC << endl;
   if( _comment )
     ff << "# Threshold of Karpov stability criterion for insertion of a phase in Selekt2() procedure" << endl;
   ff << left << setw(12) << "<pa_DFM> " <<  right << setw(8) << pa->p.DFM << endl;
   if( _comment )
   {  ff << "# Insertion amounts used after the simplex() initial approximation and in Selekt2() algorithm" << endl;
      ff << "# DFYw - water solvent;" << endl;
   }
   ff << left << setw(12) << "<pa_DFYw> " <<  right << setw(8) << pa->p.DFYw << endl;
   if( _comment )
     ff << "# DFYaq - aqueous species;" << endl;
   ff << left << setw(12) << "<pa_DFYaq> " <<  right << setw(8) << pa->p.DFYaq << endl;
   if( _comment )
     ff << "\n# DFYid - ideal solution components;" << endl;
//   ff << left << setw(12) << "<pa_PC> " <<  right << setw(8) << pa->p.PC << endl;
   ff << left << setw(12) << "<pa_DFYid> " <<  right << setw(8) << pa->p.DFYid << endl;
   if( _comment )
     ff << "# DFYr - major solution components;" << endl;
   ff << left << setw(12) << "<pa_DFYr> " <<  right << setw(8) << pa->p.DFYr << endl;
   if( _comment )
     ff << "# DFYh - minor solution components;" << endl;
   ff << left << setw(12) << "<pa_DFYh> " <<  right << setw(8) << pa->p.DFYh << endl;
   if( _comment )
     ff << "# DFYc - single-component phase; " << endl;
   ff << left << setw(12) << "<pa_DFYc> " <<  right << setw(8) << pa->p.DFYc << endl;
   if( _comment )
     ff << "# Insertion amount of single-component phase (Selekt2() algorithm only)" << endl;
   ff << left << setw(12) << "<pa_DFYs> " << right << setw(8) <<  pa->p.DFYs << endl;
   if( _comment )
   {  ff << "\n# Parameters for high-accuracy IPM-2 algorithm " << endl;
      ff << "# Number of the IPM-2 enhancement loops for high-accuracy mass balance (from 0 to 14)" << endl;
   }
   ff << left << setw(12) << "<pa_DW> " << right << setw(8) << pa->p.DW  << endl;
   if( _comment )
     ff << "# Exponent for dual-thermodynamic restoring of low amounts of solutes (+2 to -5), relative to DHBM" << endl;
   ff << left << setw(12) << "<pa_DT> " << right << setw(8) << pa->p.DT  << endl;
   if( _comment )
     ff << "\n# IPM2 balance accuracy control factor DHBM[i]/b[i]" << endl;
   ff << left << setw(12) << "<pa_GAS> " << right << setw(8) <<  pa->p.GAS << endl;
//   ff << "# 'pa_DG'        1e-30  now internal in LU decomposition procedure from JAMA-TNT" << endl << endl;
//  ff << left << setw(12) << "<pa_DG> " <<  right << setw(8) << pa->p.DG << endl;
   if( _comment )
     ff << "# Standard surface density (nm-2) for calculating activity of surface species" << endl;
   ff << left << setw(12) << "<pa_DNS> " <<  right << setw(8) << pa->p.DNS << endl;
   if( _comment )
     ff << "# Control parameter of SACT calculation in sorption/surface complexation models" << endl;
   ff << left << setw(12) << "<pa_IEPS> " <<  right << setw(8) << pa->p.IEPS << endl;
   if( _comment )
     ff << "\n# Flag for using metastability/kinetic constraints on calculated DC amounts (see dll, dul arrays)" << endl;
   ff << left << setw(12) << "<pKin> " <<  right << setw(8) << pmp->PLIM << endl;
   if( _comment )
     ff << "# Tolerance on amount of DC with two-side metastability constraints set in dll, dul (moles) " << endl;
   ff << left << setw(12) << "<pa_DKIN> " <<  right << setw(8) << pa->p.DKIN << endl;
//   ff << "# 'pa_PLLG'          0 used only in GEMS-PSI shell" << endl;
//   ff << left << setw(12) << "<pa_PLLG> " <<  right << setw(8) << pa->p.PLLG << endl;
   if( _comment )
     ff << "\n# Flag for using electroneutrality condition in GEM IPM calculations " << endl;
   ff << left << setw(12) << "<pa_PE> " <<  right << setw(8) << pa->p.PE << endl;
//   ff << "# 'E'                1" << endl;
//   ff << left << setw(12) << "<E> " <<  right << setw(8) << pmp->E << endl;
   if( _comment )
     ff << "\n# Flag for the volume balance constraint (on Vol IC)" << endl;
   ff << left << setw(12) << "<PV> " <<  right << setw(8) << pmp->PV << endl;
//   ff << "# These dimensions can be calculated from the DATACH information" << endl;
//   ff << "# 'Ls'              23" << endl;
//   ff << "# 'LO'              18" << endl;
//   ff << "# 'PG'               4" << endl;
//   ff << left << setw(12) << "<Ls> " <<  right << setw(8) << pmp->Ls << endl;
//   ff << left << setw(12) << "<LO> " <<  right << setw(8) << pmp->LO << endl;
//   ff << left << setw(12) << "<PG> " <<  right << setw(8) << pmp->PG << endl;
   if( _comment )
       ff << "\n# Total number of DCs in liquid hydrocarbon phases" << endl;
   ff << left << setw(12) << "<PSOL> " <<  right << setw(8) << pmp->PSOL << endl;
//   ff << "# Do not know if this stuff is really necessary" << endl;
//   ff << "# 'GWAT'         55.51" << endl;
//   ff << "# 'EpsW'       78.2451" << endl;
//   ff << "# 'RoW'       0.997061" << endl << endl;
//   ff << left << setw(12) << "<GWAT> " <<  right << setw(8) << pmp->GWAT << endl;
//   ff << left << setw(12) << "<EpsW> " <<  right << setw(8) << EpsW << endl;
//   ff << left << setw(12) << "<RoW>  " <<  right << setw(8) << RoW << endl;
   if( _comment )
     ff << "\n# Flag for using (+) or ignoring (-) specific surface areas of phases " << endl;
   ff << left << setw(12) << "<PAalp> " <<  right << setw(6) <<
      "\'" << PAalp << "\'" << endl;
   if( _comment )
    ff << "\n# Flag for using (+) or ignoring (-) specific surface free energies of phase interfaces " << endl;
   ff << left << setw(12) << "<PSigm> " <<  right << setw(6) <<
      "\'" << PSigm << "\'" << endl;
   if( _comment )
   {  ff << "\n## (2) Important dimensionalities that affect memory allocation" << endl;
      ff << "# Total number of dependent components in sorption phases included in this system" << endl;
   }
   ff << left << setw(12) << "<Lads> " <<  right << setw(8) << pmp->Lads << endl;
   if( _comment )
     ff << "# Number of sorption phases included in this system" << endl;
   ff << left << setw(12) << "<FIa> " <<  right << setw(8) << pmp->FIa << endl;
   if( _comment )
     ff << "# Allowed number of surface types per adsorption phase (default: 6 if FIa > 0)" << endl;
   ff << left << setw(12) << "<FIat> " <<  right << setw(8) << pmp->FIat << endl << endl;
//   ff << left << setw(12) << "<FIat> " <<  right << setw(8) << pmp->FIat << endl;
//   ff << left << setw(12) << "<sitNc> " <<  right << setw(8) << pmp->sitNcat << endl;
//   ff << left << setw(12) << "<sitNa> " <<  right << setw(8) << pmp->sitNan << endl;

//dynamic arrays
if( pm.FIs > 0 && pm.Ls > 0 )
{
  if( _comment )
  {   ff << "## (3) Initial data for multicomponent phases (see DATACH file for dimension nPHs)" << endl;
      ff << "# Codes for mixing models of multicomponent phases";
  }
  outArray( ff, "sMod", pmp->sMod[0], pmp->FIs, 6 );
  if( _comment )
  {  ff << "\n\n# Dimensions for parameters of non-ideal mixing models for each multicomponent phase" << endl;
     ff << "# Number of parameters per phase";
  }
   outArray( ff, "LsMod", pmp->LsMod, pmp->FIs);
   if( _comment )
     ff << "\n\n# Number of parameters per component of the phase";
   outArray( ff, "LsMdc", pmp->LsMdc, pmp->FIs);
   int LsModSum = 0;
   int LsMdcSum = 0;
   for(int i=0; i<pmp->FIs; i++)
   {
     LsModSum += pmp->LsMod[i];
     LsMdcSum += (pmp->LsMdc[i]*pmp->L1[i]);
   }
   if(LsModSum )
   {
     if( _comment )
        ff << "\n\n# Parameters of non-ideal mixing models for multicomponent phases ";
     outArray( ff, "PMc", pmp->PMc,  LsModSum);
   }
   if(LsMdcSum )
   {   if( _comment )
          ff << "\n\n# Parameters of non-ideal mixing models for components in the phase ";
     outArray( ff, "DMc", pmp->DMc,  LsMdcSum);
   }
}
  if( _comment )
  {  ff << "\n\n## (4) Some data arrays which are not provided in DATACH and DATABR files" << endl;
     ff << "# Full total bulk composition of the initial system (vector b) - see DATACH file for dimension nIC";
  }
   outArray( ff, "B", pmp->B,  pmp->N);
   if( _comment )
   {  ff << "\n\n# Initial data for DCs - see DATACH file for dimensions nDC, nDCs" << endl;
      ff << "# generic DC classes (asymmetric, solvent, ideal, single)";
   }
   outArray( ff, "DCCW", pmp->DCCW,  pmp->L, 1);
   if( _comment )
      ff << "\n\n# Partial pressures (fugacities) of dependent components (for setting constant chemical potentials)";
   outArray( ff, "Pparc", pmp->Pparc,  pmp->L);
 //  ff << "\n\n# This is not necessary - can be calculated from G0 ???????????";
 //  outArray( ff, "G0", pmp->G0,  pmp->L);
   if( _comment )
      ff << "\n\n# DC G0 increments for adjustments";
   outArray( ff, "GEX", pmp->GEX,  pmp->L);
   if( _comment )
      ff << "\n\n# DC Fixed (start) activity coefficients";
   outArray( ff, "lnGmf", pmp->lnGmf,  pmp->L);
   if( pmp->E )
   {
     if( _comment )
        ff << "\n\n# DC Unit formula charges - can be extracted from the stoich. matrix ????";
     outArray( ff, "EZ", pmp->EZ,  pmp->L);
   }
   if( _comment )
   {  ff << "\n\n# (5) Section for metastability/ kinetic constraints" << endl;
      ff << "# Code of metastability/kinetic constraints for DCs";
   }
   outArray( ff, "RLC", pmp->RLC, pmp->L, 1 );
   if( _comment )
     ff << "\n\n# Units of metastability/kinetic constraints for DCs (see vectors dul, dll)";
   outArray( ff, "RSC", pmp->RSC, pmp->L, 1 );
   if( _comment )
     ff << "\n\n# Full vector of lower metastability constraints on DC amounts in the system";
   outArray( ff, "DLL", pmp->DLL,  pmp->L);
   if( _comment )
     ff << "\n\n# Full vector of upper metastability constraints on DC amounts in the system";
   outArray( ff, "DUL", pmp->DUL,  pmp->L);
   if( _comment )
   {  ff << "\n\n# (6) Initial data for phases" << endl;
      ff << "\n# Specific surface areas of phases (whole list)";
   }
   outArray( ff, "Aalp", pmp->Aalp,  pmp->FI);
   if( PSigm != S_OFF )
   {
      if( _comment )
         ff << "\n\n# Specific surface free energy for phase-water interface";
      outArray( ff, "Sigw", pmp->Sigw,  pmp->FI);
      if( _comment )
         ff << "\n\n# Specific surface free energy for phase-gas interface";
      outArray( ff, "Sigg", pmp->Sigg,  pmp->FI);
   }
   if( _comment )
      ff << "\n\n# Surface energy or metastability parameters for phases";
   outArray( ff, "YOF", pmp->YOF,  pmp->FI);
   if( pm.FIat > 0 && /*pm.Lads > 0 &&Sveta 12/09/99*/ pm.FIs > 0 )
    { /* ADSORPTION AND ION EXCHANGE */
      if( _comment )
      {  ff << "\n\n# (7) Initial data for sorption" << endl;
         ff << "\n# Function of sorbent surface allocated to surface types";
      }
      outArray( ff, "Nfsp", &pmp->Nfsp[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      if( _comment )
        ff << "\n# Total maximum site  density per surface type, mkmol/g";
      outArray( ff, "MASDT", &pmp->MASDT[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      if( _comment )
        ff << "\n# Inner capacitance density parameter C1 (TLM, BSM, CCM), F/m2";
      outArray( ff, "C1", &pmp->XcapA[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      if( _comment )
        ff << "\n# Outer capacitance density parameter C2 (TLM, 3LM), F/m2";
      outArray( ff, "C2", &pmp->XcapB[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      if( _comment )
        ff << "\n# Third capacitance density parameter C3 (reserved)";
      outArray( ff, "C3", &pmp->XcapF[0][0], pmp->FIs*pmp->FIat, pmp->FIat);
      if( _comment )
        ff << "\n# Density of permanent surface type charge (mkeq/m2)";
      outArray( ff, "pCh", &pmp->Xetaf[0][0], pmp->FIs*pmp->FIat, pmp->FIat);

     if( _comment )
     {   ff << "\n# Setup of surface sites and specres: link to";
         ff << "\n# [0] surface type; [1] sorbent emd member;";
         ff << "\n# [2] surface site in surf. type; [3] surface EDL plane";
     }
      outArray( ff, "SATX", &pmp->SATX[0][0], pmp->Lads*4, 4);
      if( _comment )
      {  ff << "\n# Parameters of surface binding model:";
         ff << "\n# [0] max site density mmol/g; [1] charge allocated to 0 plane;";
         ff << "\n# [2] charge allocated to beta -or third plane; [3] Frumuin interaction parameter;";
         ff << "\n# [4] dentateness or CN; [5] reserved isoterm parameter.";
      }
      outArray( ff, "MASDJ", &pmp->MASDJ[0][0], pmp->Lads*DFCN, DFCN);
      if( _comment )
         ff << "\n# Classifier of EDL models applied to surface types.";
      outArray( ff, "SCM", pmp->SCM[0], pmp->FIs, pmp->FIat );
      if( _comment )
         ff << "\n# Classifier of applied SACT terms.";
      outArray( ff, "SACT", pmp->SATT, pmp->Lads, 1 );
      if( _comment )
         ff << "\n# Classifier of cpecies in sorption phase.";
      outArray( ff, "DCads", pmp->DCC3, pmp->Lads, 1 );
    }
//outArray( ff, "Vol", pmp->Vol,  pmp->L);
//outArray( ff, "G0", pmp->G0,  pmp->L);
//outArray( ff, "PUL", pmp->PUL,  pmp->L);
//outArray( ff, "PLL", pmp->PLL,  pmp->L);
//outArray( ff, "lnGam", pmp->lnGam,  pmp->L);
//outArray( ff, "F0", pmp->F0,  pmp->L);

/*
   if( pm.sitNcat*pm.sitNcat )
     outArray( ff, "sitE", pmp->sitE, pmp->sitNcat*pmp->sitNan );
   if( pm.sitNcat )
     outArray( ff, "sitXc", pmp->sitXcat, pmp->sitNcat );
   if( pm.sitNan )
      outArray( ff, "sitXa", pmp->sitXan, pmp->sitNan );
*/
 if( _comment )
   ff << "\n\n# End of file" << endl;

}

void TMulti::from_text_file_gemipm( const char *path )
{
  SPP_SETTING *pa = &TProfil::pm->pa;
  DATACH  *dCH = TNode::na->pCSD();
  int ii;

   //static values
   char PAalp;
   char PSigm;
   float EpsW;
   float RoW;

  memset( &pm.N, 0, 38*sizeof(short));
  memset( &pm.TC, 0, 55*sizeof(double));
  // get sizes from DATACH
  pmp->TC = pmp->P = 0;
  pmp->N = pmp->NR = dCH->nIC;
  pmp->L = dCH->nDC;
  pmp->FI = dCH->nPH;
  pmp->FIs = dCH->nPS;
  pmp->Ls = 0;
  for( ii=0; ii<dCH->nPS; ii++)
  {
    pmp->Ls += dCH->nDCinPH[ii];
    if( dCH->ccPH[ii] == 'a' )
     pmp->LO = pmp->Ls-1;
    if( dCH->ccPH[ii] == 'g' || dCH->ccPH[ii] == 'p' || dCH->ccPH[ii] == 'f')
      pmp->PG = dCH->nDCinPH[ii];
  }
  // read sizes and constants from txt file
  fstream ff( path, ios::in );
  ErrorIf( !ff.good() , path, "Fileopen error");

   gstring str;
   skipSpace( ff );
   f_getline( ff, str, '\n');
   memcpy( pmp->stkey, str.c_str(), EQ_RKLEN );
// static data
//   ff << left << setw(12) << "<pa_PC> " <<  right << setw(8) << pa->p.PC << endl;
   inArray( ff,"pa_DB" , &pa->p.DB, 1);
   inArray( ff,"pa_DHB", &pa->p.DHB, 1);
   inArray( ff,"pa_EPS" , &pa->p.EPS, 1);
   inArray( ff,"pa_DK" , &pa->p.DK, 1);
   inArray( ff,"pa_DF" , &pa->p.DF, 1);
   inArray( ff,"pa_DP", &pa->p.DP, 1);
   inArray( ff,"pa_IIM", &pa->p.IIM, 1);
   inArray( ff,"pa_PD" , &pa->p.PD, 1);
   inArray( ff,"pa_PRD" , &pa->p.PRD, 1);
   inArray( ff,"pa_AG" , &pa->p.AG, 1);
   inArray( ff,"pa_DGC" , &pa->p.DGC, 1);
   inArray( ff,"pa_PSM" , &pa->p.PSM, 1);
   inArray( ff,"pa_GAR" , &pa->p.GAR, 1);
   inArray( ff,"pa_GAH" , &pa->p.GAH, 1);
   inArray( ff,"pa_DS", &pa->p.DS, 1);
   inArray( ff,"pa_XwMin" , &pa->p.XwMin, 1);
   inArray( ff,"pa_ScMin" , &pa->p.ScMin, 1);
   inArray( ff,"pa_DcMin" , &pa->p.DcMin, 1);
   inArray( ff,"pa_PhMin" , &pa->p.PhMin, 1);
   inArray( ff,"pa_ICmin" , &pa->p.ICmin, 1);
   inArray( ff,"pa_PC" , &pa->p.PC, 1);
   inArray( ff,"pa_DFM" , &pa->p.DFM, 1);
   inArray( ff,"pa_DFYw" , &pa->p.DFYw, 1);
   inArray( ff,"pa_DFYaq" , &pa->p.DFYaq, 1);
   inArray( ff,"pa_DFYid" , &pa->p.DFYid, 1);
   inArray( ff,"pa_DFYr" , &pa->p.DFYr, 1);
   inArray( ff,"pa_DFYh" , &pa->p.DFYh, 1);
   inArray( ff,"pa_DFYc" , &pa->p.DFYc, 1);
   inArray( ff,"pa_DFYs", &pa->p.DFYs, 1);
   inArray( ff,"pa_DW", &pa->p.DW , 1);
   inArray( ff,"pa_DT", &pa->p.DT , 1);
   inArray( ff,"pa_GAS", &pa->p.GAS, 1);
//  inArray( ff,"pa_DG" , &pa->p.DG, 1);
   inArray( ff,"pa_DNS" , &pa->p.DNS, 1);
   inArray( ff,"pa_IEPS" , &pa->p.IEPS, 1);
   inArray( ff,"pKin" , &pmp->PLIM, 1);
   inArray( ff,"pa_DKIN" , &pa->p.DKIN, 1);
//   inArray( ff,"pa_PLLG" , &pa->p.PLLG, 1);
   inArray( ff,"pa_PE" , &pa->p.PE, 1);
   pmp->E = pa->p.PE;
//    inArray( ff,"E" , &pmp->E, 1);
   inArray( ff,"PV" , &pmp->PV, 1);
//   inArray( ff,"Ls" , &pmp->Ls, 1);
//   inArray( ff,"LO" , &pmp->LO, 1);
//   inArray( ff,"PG" , &pmp->PG, 1);
   inArray( ff,"PSOL" , &pmp->PSOL, 1);
//   inArray( ff,"GWAT" , &pmp->GWAT, 1);
//   inArray( ff,"EpsW" , &EpsW, 1);
//   inArray( ff,"RoW" , &RoW, 1);
   inArray( ff,"PAalp" , &PAalp, 1, 1);
   inArray( ff,"PSigm" , &PSigm, 1, 1);
   inArray( ff,"Lads" , &pmp->Lads, 1);
   inArray( ff,"FIa" , &pmp->FIa, 1);
   inArray( ff,"FIat" , &pmp->FIat, 1);
//   inArray( ff,"sitNc" , &pmp->sitNcat, 1);
//   inArray( ff,"sitNa" , &pmp->sitNan, 1);

//   if( dCH->ccPH[0] == PH_AQUEL )
//   {
//     RoW = dCH->roW[0];
//     EpsW = dCH->epsW[0];
//   }
//   else
//  {
    RoW = 0.99706137180;
    EpsW = 78.245147705;
//  }

#ifndef IPMGEMPLUGIN
//   syp->PAalp = PAalp;
//   syp->PSigm = PSigm;
#else
   PAalp_ = PAalp;
   PSigm_ = PSigm;
   EpsW_ = EpsW;
   RoW_ =  RoW;
#endif

   //realloc memory
#ifdef IPMGEMPLUGIN
   multi_realloc( PAalp, PSigm );
#endif

// get dynamic data from DATACH file
  for( ii=0; ii<dCH->nPH; ii++)
    pmp->L1[ii] = dCH->nDCinPH[ii];
  memcpy( pmp->A, dCH->A, dCH->nIC*dCH->nDC*sizeof(float));
  for( ii=0; ii< dCH->nIC; ii++ )
   pmp->Awt[ii]  = dCH->ICmm[ii];
  memcpy( pmp->MM, dCH->DCmm, dCH->nDC*sizeof(double));

  memset( pmp->SB, ' ', MaxICN*dCH->nIC*sizeof(char));
  for( ii=0; ii< dCH->nIC; ii++ )
  {   memcpy( pmp->SB[ii], dCH->ICNL[ii], MaxICN*sizeof(char) );
      pmp->SB[ii][MaxICN] = dCH->ccIC[ii];
  }

  memcpy( pmp->SM, dCH->DCNL, MaxDCN*dCH->nDC*sizeof(char));

  memset( pmp->SF, ' ', MaxPHN*dCH->nPH*sizeof(char) );
  for( ii=0; ii< dCH->nPH; ii++ )
  {  memcpy( pmp->SF[ii]+4, dCH->PHNL[ii], MaxPHN*sizeof(char));
     pmp->SF[ii][0] = dCH->ccPH[ii];
  }

  memcpy( pmp->ICC, dCH->ccIC, dCH->nIC*sizeof(char));
  memcpy( pmp->DCC, dCH->ccDC, dCH->nDC*sizeof(char));
// !!!!  memcpy( pmp->DCCW, dCH->ccDCW, dCH->nDC*sizeof(char));
  memcpy( pmp->PHC, dCH->ccPH, dCH->nPH*sizeof(char));

//read dynamic values from txt file
if( pm.FIs > 0 && pm.Ls > 0 )
{
   inArray( ff, "sMod", pmp->sMod[0], pmp->FIs, 6 );
   inArray( ff, "LsMod", pmp->LsMod, pmp->FIs);
   inArray( ff, "LsMdc", pmp->LsMdc, pmp->FIs);
   int LsModSum = 0;
   int LsMdcSum = 0;
   for(int i=0; i<pmp->FIs; i++)
   {
     LsModSum += pmp->LsMod[i];
     LsMdcSum += (pmp->LsMdc[i]*pmp->L1[i]);
   }
   if(LsModSum )
      inArray( ff, "PMc", pmp->PMc,  LsModSum);
   if(LsMdcSum )
     inArray( ff, "DMc", pmp->DMc,  LsMdcSum);
   }
   inArray( ff, "B", pmp->B,  pmp->N);
   inArray( ff, "DCCW", pmp->DCCW,  pmp->L, 1);
   inArray( ff, "Pparc", pmp->Pparc,  pmp->L);
//   inArray( ff, "G0", pmp->G0,  pmp->L);
   inArray( ff, "GEX", pmp->GEX,  pmp->L);
   inArray( ff, "lnGmf", pmp->lnGmf,  pmp->L);
   if( pmp->E )
     inArray( ff, "EZ", pmp->EZ,  pmp->L);
   inArray( ff, "RLC", pmp->RLC, pmp->L, 1 );
   inArray( ff, "RSC", pmp->RSC, pmp->L, 1 );
   inArray( ff, "DLL", pmp->DLL,  pmp->L);
   inArray( ff, "DUL", pmp->DUL,  pmp->L);
   inArray( ff, "Aalp", pmp->Aalp,  pmp->FI);
   if( PSigm != S_OFF )
   {
      inArray( ff, "Sigw", pmp->Sigw,  pmp->FI);
      inArray( ff, "Sigg", pmp->Sigg,  pmp->FI);
   }
   inArray( ff, "YOF", pmp->YOF,  pmp->FI);
   if( pm.FIat > 0 && /*pm.Lads > 0 &&Sveta 12/09/99*/ pm.FIs > 0 )
    { /* ADSORPTION AND ION EXCHANGE */
      inArray( ff, "Nfsp", &pmp->Nfsp[0][0], pmp->FIs*pmp->FIat);
      inArray( ff, "MASDT", &pmp->MASDT[0][0], pmp->FIs*pmp->FIat);
      inArray( ff, "C1", &pmp->XcapA[0][0], pmp->FIs*pmp->FIat);
      inArray( ff, "C2", &pmp->XcapB[0][0], pmp->FIs*pmp->FIat);
      inArray( ff, "C3", &pmp->XcapF[0][0], pmp->FIs*pmp->FIat);
      inArray( ff, "pCh", &pmp->Xetaf[0][0], pmp->FIs*pmp->FIat);

      inArray( ff, "SATX", &pmp->SATX[0][0], pmp->Lads*4);
      inArray( ff, "MASDJ", &pmp->MASDJ[0][0], pmp->Lads*DFCN);
      inArray( ff, "SCM", pmp->SCM[0], pmp->FIs, pmp->FIat );
      inArray( ff, "SACT", pmp->SATT, pmp->Lads, 1 );
      inArray( ff, "DCads", pmp->DCC3, pmp->Lads, 1 );
    }
//    inArray( ff, "Vol", pmp->Vol,  pmp->L);
//    inArray( ff, "G0", pmp->G0,  pmp->L);
//    inArray( ff, "PUL", pmp->PUL,  pmp->L);
//    inArray( ff, "PLL", pmp->PLL,  pmp->L);
//    inArray( ff, "lnGam", pmp->lnGam,  pmp->L);
//    inArray( ff, "F0", pmp->F0,  pmp->L);

/*
   if( pm.sitNcat*pm.sitNcat )
     inArray( ff, "sitE", pmp->sitE, pmp->sitNcat*pmp->sitNan );
   if( pm.sitNcat )
     inArray( ff, "sitXc", pmp->sitXcat, pmp->sitNcat );
   if( pm.sitNan )
     inArray( ff, "sitXa", pmp->sitXan, pmp->sitNan );
*/
}

//=============================================================================
