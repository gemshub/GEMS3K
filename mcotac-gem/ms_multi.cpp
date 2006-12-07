//-------------------------------------------------------------------
// $Id: ms_multi.cpp 710 2006-06-07 14:58:15Z gems $
//
// Implementation of TMulti class, config functions
//
// Rewritten from C to C++ by S.Dmytriyeva
// Copyright (C) 1995-2001 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
// Uses: GEM-Vizor GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://les.web.psi.ch/Software/GEMS-PSI for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//
#include "m_param.h"
#include "v_object.h"
#include "visor.h"


TMulti::TMulti( int nrt, SYSTEM* sy_, MTPARM *tp_, RMULTS *mu_ ):
        TSubModule( nrt ), syp(sy_), tpp(tp_), mup(mu_)
{
    set_def();
    pmp = &pm;
}


// link values to objects
void TMulti::ods_link( int /*q*/)
{
    //aObj[ o_wd_stkey].SetPtr( pm.stkey );
    //aObj[ o_wd_nv_].SetPtr( pm.NV_ );
    aObj[ o_wi_pun ].SetPtr( &pm.PunE ); /*a 4*/
    aObj[ o_wi_nnr ].SetPtr( &pm.N );    /*i 2*/
    aObj[ o_wi_ldim ].SetPtr( &pm.L );   /*i 6*/
    aObj[ o_wi_fdim ].SetPtr( &pm.FI );  /*i 3*/
    aObj[ o_wo_fdim ].SetPtr( &pm.FI1 ); /*i 3*/
    aObj[ o_wo_it ].SetPtr( &pm.IT );
    aObj[ o_wi_mode ].SetPtr( &pm.E );   /*i 4*/
    aObj[ o_wo_ipar ].SetPtr( &pm.Ec );  /*i 17*/
    aObj[ o_wio_tc ].SetPtr( &pm.TC );   /*d 2*/
    aObj[ o_wio_t ].SetPtr( &pm.T );     /*d 2*/
    aObj[ o_wio_p ].SetPtr( &pm.P );     /*d 2*/
    aObj[ o_wo_vx ].SetPtr( &pm.VX_ );   /*d 2*/
    aObj[ o_wo_gx ].SetPtr( &pm.GX_ );   /*d 2*/
    aObj[ o_wo_ax ].SetPtr( &pm.AX_ );   /*d 2*/
    aObj[ o_wo_ux ].SetPtr( &pm.UX_ );   /*d 2*/
    aObj[ o_wo_hx ].SetPtr( &pm.HX_ );   /*d 2*/
    aObj[ o_wo_sx ].SetPtr( &pm.SX_ );   /*d 2*/
    aObj[ o_wo_cpx ].SetPtr( &pm.CpX_ ); /*d 2*/
    aObj[ o_wo_cvx ].SetPtr( &pm.CvX_ ); /*d 2*/
    aObj[ o_wi_t0 ].SetPtr( &pm.T0 );
    aObj[ o_wi_ve ].SetPtr( &pm.VE );
    aObj[ o_wi_mbx ].SetPtr( &pm.MBX );
    aObj[ o_wo_fx ].SetPtr( &pm.FX );
    aObj[ o_wo_ic ].SetPtr( &pm.IC );
    aObj[ o_wd_ph ].SetPtr( &pm.pH );
    aObj[ o_wd_pe ].SetPtr( &pm.pe );
    aObj[ o_wd_eh ].SetPtr( &pm.Eh );
    aObj[ o_wi_mol ].SetPtr( &pm.DHBM );
    aObj[ o_wi_gwat ].SetPtr( &pm.GWAT );
    aObj[ o_wi_ymet ].SetPtr( &pm.YMET );
    aObj[ o_wio_denw ].SetPtr(&pm.denW ); /*d 2*/
    aObj[ o_wio_epsw ].SetPtr(&pm.epsW ); /*d 2*/
    aObj[ o_wo_psi ].SetPtr( &pm.PCI );   /*d 2*/
    aObj[ o_wo_lnp ].SetPtr( &pm.lnP );
    aObj[ o_wo_rt ].SetPtr( &pm.RT );
    aObj[ o_wo_frt ].SetPtr( &pm.FRT );
    aObj[ o_wo_yw ].SetPtr( &pm.Yw );
    aObj[ o_wio_cons ].SetPtr(&pm.ln5551 );/*d 6*/
    aObj[ o_wio_fitv ].SetPtr(pm.FitVar ); /*d 5*/
    aObj[ o_wd_sitn ].SetPtr(&pm.sitNcat );

    // dynamic part 1
    aObj[ o_wi_l1 ].SetPtr( pm.L1 );
    aObj[ o_wi_l1 ].SetDim( pm.FI, 1 );
    aObj[ o_wi_lsmod ].SetPtr(pm.LsMod );
    aObj[ o_wi_lsmod ].SetDim( pm.FIs, 1 );
    aObj[ o_wi_lsmdc ].SetPtr(pm.LsMdc );
    aObj[ o_wi_lsmdc ].SetDim( pm.FIs, 1 );
    aObj[ o_wi_mui ].SetPtr( pm.mui );
    aObj[ o_wi_mui ].SetDim( pm.N, 1 );
    aObj[ o_wi_muk ].SetPtr( pm.muk );
    aObj[ o_wi_muk ].SetDim( pm.FI, 1 );
    aObj[ o_wi_muj ].SetPtr( pm.muj );
    aObj[ o_wi_muj ].SetDim( pm.L, 1 );
    aObj[ o_wi_dul ].SetPtr( pm.DUL );
    aObj[ o_wi_dul ].SetDim( pm.L, 1 );
    aObj[ o_wi_dll ].SetPtr( pm.DLL );
    aObj[ o_wi_dll ].SetDim( pm.L, 1 );
    aObj[ o_wi_gex ].SetPtr( pm.GEX );
    aObj[ o_wi_gex ].SetDim( pm.L, 1 );
    aObj[ o_wi_pul ].SetPtr( pm.PUL );
    aObj[ o_wi_pul ].SetDim( pm.FIs, 1 );
    aObj[ o_wi_pll ].SetPtr( pm.PLL );
    aObj[ o_wi_pll ].SetDim( pm.FIs, 1 );
    aObj[ o_wi_yof ].SetPtr( pm.YOF );
    aObj[ o_wi_yof ].SetDim( pm.FI, 1 );
    /* SIZE OF PMc- CALC M1=1 */
    aObj[ o_wi_pmc ].SetPtr( pm.PMc );  // OBsetDim( o_wi_pmc,  M1, 1 );
    /* SIZE OF DMc- CALC M2=1 */
    aObj[ o_wi_dmc ].SetPtr( pm.DMc );  // OBsetDim( o_wi_Dmc,  M2, 1 );
    aObj[ o_wio_vol ].SetPtr( pm.Vol );
    aObj[ o_wio_vol ].SetDim( pm.L, 1 );
    aObj[ o_w_hym ].SetPtr( pm.HYM );
    aObj[ o_w_hym ].SetDim( pm.L, 1 );
    aObj[ o_w_vl ].SetPtr( pm.VL );
    aObj[ o_w_vl ].SetDim( pm.L, 1 );
    aObj[ o_wi_mm ].SetPtr( pm.MM );
    aObj[ o_wi_mm ].SetDim( pm.L, 1 );
    aObj[ o_wio_h0 ].SetPtr( pm.H0 );
    aObj[ o_wio_h0 ].SetDim( pm.L, 1 );
    aObj[ o_wio_a0 ].SetPtr( pm.A0 );
    aObj[ o_wio_a0 ].SetDim( pm.L, 1 );
    aObj[ o_wio_u0 ].SetPtr( pm.U0 );
    aObj[ o_wio_u0 ].SetDim( pm.L, 1 );
    aObj[ o_wio_s0 ].SetPtr( pm.S0 );
    aObj[ o_wio_s0 ].SetDim( pm.L, 1 );
    aObj[ o_wio_cp0 ].SetPtr( pm.Cp0 );
    aObj[ o_wio_cp0 ].SetDim( pm.L, 1 );
    aObj[ o_wio_cv0 ].SetPtr( pm.Cv0 );
    aObj[ o_wio_cv0 ].SetDim( pm.L, 1 );
    aObj[ o_wi_pparc ].SetPtr(pm.Pparc );
    aObj[ o_wi_pparc ].SetDim( pm.L, 1 );

    // part 2
    aObj[o_wi_satndx].SetPtr( pm.SATX );  // changed 27.10.2004 for CD MUSIC
    if( pm.FIat > 0 ) aObj[ o_wi_satndx ].SetDim( pm.Lads, 4 );
//
    aObj[ o_wd_ym].SetPtr( pm.Y_m );
    aObj[ o_wd_ym ].SetDim( pm.L, 1 );
    aObj[ o_wd_yla].SetPtr( pm.Y_la );
    aObj[ o_wd_yla ].SetDim( pm.L, 1 );
    aObj[ o_wd_yw].SetPtr( pm.Y_w );
    aObj[ o_wd_yw ].SetDim( pm.L, 1 );
    aObj[ o_wd_gamma].SetPtr( pm.Gamma );
    aObj[ o_wd_gamma ].SetDim( pm.L, 1);
    aObj[ o_wi_lngmf].SetPtr( pm.lnGmf );
    aObj[ o_wi_lngmf ].SetDim( pm.L, 1);
    aObj[ o_wio_lngmm].SetPtr(pm.lnGmM );
    aObj[ o_wio_lngmm ].SetDim(pm.L, 1);
    aObj[ o_wi_ez].SetPtr( pm.EZ );
    aObj[ o_wi_ez ].SetDim( pm.L, 1 );
    aObj[ o_wi_aalp].SetPtr( pm.Aalp );
    aObj[ o_wi_aalp ].SetDim( pm.FI, 1 );
    aObj[ o_wi_sigw].SetPtr( pm.Sigw );
    aObj[ o_wi_sigw ].SetDim( pm.FI, 1 );
    aObj[ o_wi_sigg].SetPtr( pm.Sigg );
    aObj[ o_wi_sigg ].SetDim( pm.FI, 1 );
    aObj[ o_wi_nfsp].SetPtr( pm.Nfsp );
    aObj[o_wi_nfsp].SetDim(pm.FIs,pm.FIat);
    aObj[ o_wi_masdt].SetPtr( pm.MASDT );
    aObj[o_wi_masdt].SetDim(pm.FIs,pm.FIat);
    //    aObj[ o_wi_misdt].SetPtr( pm.MISDT );  aObj[o_wi_misdt].SetDim(pm.FIs,pm.FIat);
    aObj[ o_wd_fvol].SetPtr( pm.FVOL );
    aObj[ o_wd_fvol ].SetDim( pm.FI, 1 );
    aObj[ o_wd_fwgt].SetPtr( pm.FWGT );
    aObj[ o_wd_fwgt ].SetDim( pm.FI, 1 );
    aObj[ o_w_xcapa].SetPtr( pm.XcapA );
    aObj[o_w_xcapa].SetDim(pm.FIs,pm.FIat);
    aObj[ o_w_xcapb].SetPtr( pm.XcapB );
    aObj[o_w_xcapb].SetDim(pm.FIs,pm.FIat);
    aObj[ o_w_xcapd].SetPtr( pm.XcapD );
    aObj[o_w_xcapd].SetDim(pm.FIs,pm.FIat);
    aObj[ o_w_xcapf].SetPtr( pm.XcapF );
    aObj[o_w_xcapf].SetDim(pm.FIs,pm.FIat);
    aObj[ o_wo_xdla].SetPtr( pm.XdlA );
    aObj[ o_wo_xdla].SetDim(pm.FIs,pm.FIat);
    aObj[ o_wo_xdlb].SetPtr( pm.XdlB );
    aObj[ o_wo_xdlb].SetDim(pm.FIs,pm.FIat);
    aObj[ o_wo_xdld].SetPtr( pm.XdlD );
    aObj[ o_wo_xdld].SetDim(pm.FIs,pm.FIat);
    aObj[ o_wo_psia].SetPtr( pm.XpsiA );
    aObj[o_wo_psia].SetDim(pm.FIs,pm.FIat);
    aObj[ o_wo_psib].SetPtr( pm.XpsiB );
    aObj[o_wo_psib].SetDim(pm.FIs,pm.FIat);
    aObj[ o_wo_psid].SetPtr( pm.XpsiD );
    aObj[o_wo_psid].SetDim(pm.FIs,pm.FIat);
    aObj[ o_wi_xr0h0].SetPtr( pm.Xr0h0 );
    aObj[o_wi_xr0h0 ].SetDim(pm.FI, 2);
    aObj[ o_wi_xlama].SetPtr( pm.XlamA );
    aObj[o_wi_xlama].SetDim(pm.FIs,pm.FIat);
    aObj[ o_wi_xetaf].SetPtr( pm.Xetaf );
    aObj[o_wi_xetaf].SetDim(pm.FIs,pm.FIat);
    aObj[ o_wi_xcond].SetPtr( pm.Xcond );
    aObj[ o_wi_xcond ].SetDim(pm.FI, 1 );
    aObj[ o_wi_xeps].SetPtr( pm.Xeps );
    aObj[ o_wi_xeps ].SetDim( pm.FI, 1 );
    aObj[ o_wi_awt].SetPtr( pm.Awt );
    aObj[ o_wi_awt ].SetDim( pm.N, 1 );
    aObj[ o_wi_a].SetPtr( pm.A );
    aObj[ o_wi_a ].SetDim( pm.L, pm.N );
    aObj[ o_wd_xfs].SetPtr( pm.XFs );
    aObj[ o_wd_xfs ].SetDim( pm.FI, 1 );
    aObj[ o_wd_falps].SetPtr( pm.Falps );
    aObj[ o_wd_falps ].SetDim( pm.FI, 1 );
    aObj[ o_wo_g].SetPtr( pm.G );
    aObj[ o_wo_g ].SetDim( pm.L, 1 );
    aObj[ o_wi_g0].SetPtr( pm.G0 );
    aObj[ o_wi_g0 ].SetDim( pm.L, 1 );
    aObj[ o_wo_lngam].SetPtr( pm.lnGam );
    aObj[ o_wo_lngam ].SetDim( pm.L, 1);
    aObj[ o_w_lngmo].SetPtr( pm.lnGmo );
    aObj[ o_w_lngmo ].SetDim( pm.L, 1 );
    aObj[ o_wi_b].SetPtr( pm.B );
    aObj[ o_wi_b ].SetDim( pm.N, 1 );
    aObj[ o_wo_u].SetPtr( pm.U );
    aObj[ o_wo_u ].SetDim( pm.N, 1 );
    aObj[ o_wd_ur].SetPtr( pm.U_r );
    aObj[ o_wd_ur ].SetDim( pm.N, 1 );
    aObj[ o_wo_c].SetPtr( pm.C );
    aObj[ o_wo_c ].SetDim( pm.N, 1 );
    aObj[ o_wd_icm].SetPtr( pm.IC_m );
    aObj[ o_wd_icm ].SetDim( pm.N, 1 );
    aObj[ o_wd_iclcm].SetPtr( pm.IC_lm );
    aObj[ o_wd_iclcm ].SetDim( pm.N, 1 );
    aObj[ o_wd_icwm].SetPtr( pm.IC_wm );
    aObj[ o_wd_icwm ].SetDim( pm.N, 1 );

    // part 3
    aObj[ o_w_xfts].SetPtr( pm.XFTS );
    aObj[ o_w_xfts ].SetDim(pm.FIs,pm.FIat);
//  Changed 27.10.2004 for CD MUSIC
    aObj[ o_wo_masdj].SetPtr( pm.MASDJ );
    aObj[ o_wo_masdj ].SetDim( pm.Lads, DFCN );
//    aObj[ o_wo_masdj ].SetDim( pm.Ls, 1 );
    aObj[ o_wo_bf].SetPtr( pm.BF );
    aObj[ o_wo_bf ].SetDim( pm.FIs, pm.N );
    aObj[ o_w_xf].SetPtr( pm.XF );
    aObj[ o_w_xf ].SetDim( pm.FI, 1 );
    aObj[ o_wo_yf].SetPtr( pm.YF );
    aObj[ o_wo_yf ].SetDim( pm.FI, 1 );
    aObj[ o_ww_xfa].SetPtr( pm.XFA );
    aObj[ o_ww_xfa ].SetDim( pm.FIs, 1 );
    aObj[ o_ww_yfa].SetPtr( pm.YFA );
    aObj[ o_ww_yfa ].SetDim( pm.FIs, 1 );
    aObj[ o_wo_falp].SetPtr( pm.Falp );
    aObj[ o_wo_falp ].SetDim( pm.FI, 1 );
    aObj[ o_w_xetaa].SetPtr( pm.XetaA );
    aObj[ o_w_xetaa].SetDim(pm.FIs,pm.FIat);
    aObj[ o_w_xetab].SetPtr( pm.XetaB );
    aObj[ o_w_xetab].SetDim(pm.FIs,pm.FIat);
aObj[ o_w_xetad].SetPtr( pm.XetaD );    // Added 12.09.05 KD
aObj[ o_w_xetad].SetDim(pm.FIs,pm.FIat);
    aObj[ o_w_x].SetPtr( pm.X );
    aObj[ o_w_x ].SetDim( pm.L, 1 );
    aObj[ o_wo_y].SetPtr( pm.Y );
    aObj[ o_wo_y ].SetDim( pm.L, 1 );
    aObj[ o_w_xy].SetPtr( pm.XY );
    aObj[ o_w_xy ].SetDim( pm.L, 1 );
    aObj[ o_wi_qp].SetPtr( pm.Qp );        // aObj[ o_wi_qp ].SetDim( pm., 1 );
    aObj[ o_wi_qd].SetPtr( pm.Qd );        // aObj[ o_wi_qd ].SetDim( pm., 1 );
    aObj[ o_wo_mu].SetPtr( pm.MU );
    aObj[ o_wo_mu ].SetDim( pm.L, 1 );
    aObj[ o_w_emu].SetPtr( pm.EMU );
    aObj[ o_w_emu ].SetDim( pm.L, 1 );
    aObj[ o_w_nmu].SetPtr( pm.NMU );
    aObj[ o_w_nmu ].SetDim( pm.L, 1 );
    aObj[ o_w_w].SetPtr( pm.W );
    aObj[ o_w_w ].SetDim( pm.L, 1 );
    aObj[ o_wo_fxa].SetPtr( pm.Fx );
    aObj[ o_wo_fxa ].SetDim( pm.L, 1 );
    aObj[ o_wo_wx].SetPtr( pm.Wx );
    aObj[ o_wo_wx ].SetDim( pm.L, 1 );
    aObj[ o_wo_f].SetPtr( pm.F );
    aObj[ o_wo_f ].SetDim( pm.L, 1 );
    aObj[ o_w_f0].SetPtr( pm.F0 );
    aObj[ o_w_f0 ].SetDim( pm.L, 1 );
// Changed for CD MUSIC 27.10.2004 KD
    aObj[ o_w_d ].SetPtr( pm.D );
    aObj[ o_w_d ].SetDim( MST, MST );
//    aObj[ o_w_d].SetPtr( pm.D );
//    aObj[ o_w_d ].SetDim( pm.L, 1 );
//    aObj[ o_w_r].SetPtr( pm.R );
//    aObj[ o_w_r ].SetDim( pm.NR, pm.NR+1 );
//    aObj[ o_w_r1].SetPtr( pm.R1 );
//    aObj[ o_w_r1 ].SetDim( pm.NR, pm.NR+1 );
    aObj[ o_wi_smod ].SetPtr( pm.sMod );
    aObj[ o_wi_smod ].SetDim( pm.FIs, 1 );
    aObj[ o_wd_sb ].SetPtr( pm.SB );
    aObj[ o_wd_sb ].SetDim( pm.N, 1 );
    aObj[ o_w_sbh ].SetPtr( pm.SB1 );
    aObj[ o_w_sbh ].SetDim( 1, pm.N );
    aObj[ o_wd_sm ].SetPtr( pm.SM );
    aObj[ o_wd_sm ].SetDim( pm.L, 1 );
    aObj[ o_wd_sf ].SetPtr( pm.SF );
    aObj[ o_wd_sf ].SetDim( pm.FI, 1 );
    aObj[ o_wd_sm2 ].SetPtr( pm.SM2 );
    aObj[ o_wd_sm2 ].SetDim( pm.Ls, 1 );
// Added 27.10.2004 for CD MUSIC models
    aObj[ o_wd_sm3 ].SetPtr( pm.SM3 );
    aObj[ o_wd_sm3 ].SetDim( pm.Lads, 1 );
//
    aObj[ o_wd_sf2].SetPtr( pm.SF2 );
    aObj[ o_wd_sf2 ].SetDim( pm.FIs, 1 );
    aObj[ o_wd_sfs].SetPtr( pm.SFs );
    aObj[ o_wd_sfs ].SetDim( pm.FI, 1 );
    aObj[ o_wd_txtb].SetPtr( pm.pbuf );    //aObj[o_wd_txtb].SetDim(1,size(pm.pbuf));
    aObj[ o_wi_rlc].SetPtr( pm.RLC );
    aObj[ o_wi_rlc ].SetDim( pm.L, 1 );
    aObj[ o_wi_rsc].SetPtr( pm.RSC );
    aObj[ o_wi_rsc ].SetDim( pm.L, 1 );
    aObj[ o_wi_rflc].SetPtr( pm.RFLC );
    aObj[ o_wi_rflc ].SetDim( pm.FIs, 1 );
    aObj[ o_wi_rfsc].SetPtr( pm.RFSC );
    aObj[ o_wi_rfsc ].SetDim( pm.FIs, 1 );
    aObj[ o_wi_icc].SetPtr( pm.ICC );
    aObj[ o_wi_icc ].SetDim( pm.N, 1 );
    aObj[ o_wi_dcc].SetPtr( pm.DCC );
    aObj[ o_wi_dcc ].SetDim( pm.L, 1 );
//  Added 27.10.2004 for CD MUSIC models
    aObj[ o_wi_dcc3 ].SetPtr( pm.DCC3 );
    aObj[ o_wi_dcc3 ].SetDim( pm.Lads, 1 );
//
    aObj[ o_wi_scm].SetPtr( pm.SCM );
    aObj[ o_wi_scm ].SetDim(pm.FIs,pm.FIat);
// Changed 27.10.2004
    aObj[ o_wi_satt].SetPtr( pm.SATT );
    aObj[ o_wi_satt ].SetDim( pm.Lads, 1 );
//
    aObj[ o_wi_dccw].SetPtr( pm.DCCW );
    aObj[ o_wi_dccw ].SetDim( pm.L, 1 );
    aObj[ o_wi_phc].SetPtr( pm.PHC );
    aObj[ o_wi_phc ].SetDim( pm.FI, 1 );
    aObj[ o_wd_fug].SetPtr( pm.Fug );
    aObj[ o_wd_fug ].SetDim( pm.PG, 1 );
    aObj[ o_wd_fugl].SetPtr( pm.Fug_l );
    aObj[ o_wd_fugl ].SetDim( pm.PG, 1 );
    aObj[ o_wd_ppgl].SetPtr( pm.Ppg_l );
    aObj[ o_wd_ppgl ].SetDim( pm.PG, 1 );
    aObj[ o_wio_wb].SetPtr( pm.Wb );
    aObj[ o_wio_wb ].SetDim( pm.Ls, 1 );
    aObj[ o_wio_wabs].SetPtr( pm.Wabs );
    aObj[ o_wio_wabs ].SetDim( pm.Ls, 1 );
    aObj[ o_wi_rion].SetPtr( pm.Rion );
    aObj[ o_wi_rion ].SetDim( pm.Ls, 1 );
//    aObj[ o_wo_lnsat].SetPtr( pm.lnSAT );
//    aObj[ o_wo_lnsat ].SetDim( pm.Ls, 1 );
aObj[ o_wo_lnsat].SetPtr( pm.lnSAC );
aObj[ o_wo_lnsat ].SetDim( pm.Lads, 4 );

//  Added 16.11.2004 by Sveta
    aObj[ o_wd_sitxcat ].SetPtr( pm.sitXcat );
    aObj[ o_wd_sitxcat ].SetDim( pm.sitNcat, 1 );
    aObj[ o_wd_sitxan ].SetPtr( pm.sitXan );
    aObj[ o_wd_sitxan ].SetDim( 1, pm.sitNan );
    aObj[ o_wd_site ].SetPtr( pm.sitE );
    aObj[ o_wd_site ].SetDim( pm.sitNcat, pm.sitNan );
}

// set dynamic Objects ptr to values
void TMulti::dyn_set(int /*q*/)
{

   pm.L1    = (short *)aObj[ o_wi_l1 ].GetPtr();
   pm.LsMod = (short *)aObj[ o_wi_lsmod ].GetPtr();
    pm.LsMdc = (short *)aObj[ o_wi_lsmdc ].GetPtr();
    pm.mui   = (short *)aObj[ o_wi_mui ].GetPtr();
    pm.muk   = (short *)aObj[ o_wi_muk ].GetPtr();
    pm.muj   = (short *)aObj[ o_wi_muj ].GetPtr();
pm.SATX =(short (*)[4])aObj[ o_wi_satndx ].GetPtr();
    pm.DUL   = (double *)aObj[ o_wi_dul ].GetPtr();
    pm.DLL   = (double *)aObj[ o_wi_dll ].GetPtr();
    pm.GEX   = (double *)aObj[ o_wi_gex ].GetPtr();
    pm.PUL   = (double *)aObj[ o_wi_pul ].GetPtr();
    pm.PLL   = (double *)aObj[ o_wi_pll ].GetPtr();
    pm.YOF   = (double *)aObj[ o_wi_yof ].GetPtr();
    pm.PMc   = (float *)aObj[ o_wi_pmc ].GetPtr();
    pm.DMc   = (float *)aObj[ o_wi_dmc ].GetPtr();
    pm.Vol   = (double *)aObj[ o_wio_vol ].GetPtr();
    pm.HYM   = (float *)aObj[ o_w_hym ].GetPtr();
    pm.VL    = (float *)aObj[ o_w_vl ].GetPtr();
    pm.MM    = (double *)aObj[ o_wi_mm ].GetPtr();
    pm.H0    = (float *)aObj[ o_wio_h0 ].GetPtr();
    pm.A0    = (float *)aObj[ o_wio_a0 ].GetPtr();
    pm.U0    = (float *)aObj[ o_wio_u0 ].GetPtr();
    pm.S0    = (float *)aObj[ o_wio_s0 ].GetPtr();
    pm.Cp0   = (float *)aObj[ o_wio_cp0 ].GetPtr();
    pm.Cv0   = (float *)aObj[ o_wio_cv0 ].GetPtr();
    pm.Pparc = (double *)aObj[ o_wi_pparc ].GetPtr();
    pm.Y_m   = (double *)aObj[ o_wd_ym ].GetPtr();
    pm.Y_la  = (double *)aObj[ o_wd_yla ].GetPtr();
    pm.Y_w   = (double *)aObj[ o_wd_yw ].GetPtr();
    pm.Gamma = (double *)aObj[ o_wd_gamma ].GetPtr();
    pm.lnGmf = (double *)aObj[ o_wi_lngmf ].GetPtr();
    pm.lnGmM = (double *)aObj[ o_wio_lngmm ].GetPtr();
    pm.EZ    = (double *)aObj[ o_wi_ez ].GetPtr();
    pm.Wb    = (float *)aObj[ o_wio_wb ].GetPtr();
    pm.Wabs  = (float *)aObj[ o_wio_wabs ].GetPtr();
    pm.Rion  = (float *)aObj[ o_wi_rion ].GetPtr();
    pm.Aalp  = (float *)aObj[ o_wi_aalp ].GetPtr();
    pm.Sigw  = (float *)aObj[ o_wi_sigw ].GetPtr();
    pm.Sigg  = (float *)aObj[ o_wi_sigg ].GetPtr();
    pm.Nfsp  = (float (*)[MST])aObj[ o_wi_nfsp ].GetPtr();
    pm.MASDT = (float (*)[MST])aObj[ o_wi_masdt ].GetPtr();
    pm.FVOL  = (double *)aObj[ o_wd_fvol ].GetPtr();
    pm.FWGT  = (double *)aObj[ o_wd_fwgt ].GetPtr();
    pm.XcapA = (float (*)[MST])aObj[ o_w_xcapa ].GetPtr();
    pm.XcapB = (float (*)[MST])aObj[ o_w_xcapb ].GetPtr();
    pm.XcapD = (float (*)[MST])aObj[ o_w_xcapd ].GetPtr();
    pm.XcapF = (float (*)[MST])aObj[ o_w_xcapf ].GetPtr();
    pm.XdlA  = (float (*)[MST])aObj[ o_wo_xdla ].GetPtr();
    pm.XdlB  = (float (*)[MST])aObj[ o_wo_xdlb ].GetPtr();
    pm.XdlD  = (float (*)[MST])aObj[ o_wo_xdld ].GetPtr();
    pm.XpsiA = (double (*)[MST])aObj[ o_wo_psia ].GetPtr();
    pm.XpsiB = (double (*)[MST])aObj[ o_wo_psib ].GetPtr();
    pm.XpsiD = (double (*)[MST])aObj[ o_wo_psid ].GetPtr();
    pm.Xr0h0 = (float (*)[2])aObj[ o_wi_xr0h0 ].GetPtr();
    pm.XlamA = (float (*)[MST])aObj[ o_wi_xlama ].GetPtr();
    pm.Xetaf = (float (*)[MST])aObj[ o_wi_xetaf ].GetPtr();
    pm.Xcond = (float *)aObj[ o_wi_xcond ].GetPtr();
    pm.Xeps  = (float *)aObj[ o_wi_xeps ].GetPtr();
    pm.Awt   = (float *)aObj[ o_wi_awt ].GetPtr();
    pm.A     = (float *)aObj[ o_wi_a ].GetPtr();
    pm.XFs   = (float *)aObj[ o_wd_xfs ].GetPtr();
    pm.Falps = (float *)aObj[ o_wd_falps ].GetPtr();
    pm.Fug   = (float *)aObj[ o_wd_fug ].GetPtr();
    pm.Fug_l = (float *)aObj[ o_wd_fugl ].GetPtr();
    pm.Ppg_l = (float *)aObj[ o_wd_ppgl ].GetPtr();
    pm.XFTS  = (double (*)[MST])aObj[ o_w_xfts ].GetPtr();
    pm.MASDJ = (float (*)[DFCN])aObj[ o_wo_masdj ].GetPtr();
    pm.G     = (double *)aObj[ o_wo_g ].GetPtr();
    pm.G0    = (double *)aObj[ o_wi_g0 ].GetPtr();
    pm.lnGam = (double *)aObj[ o_wo_lngam ].GetPtr();
    pm.lnGmo = (double *)aObj[ o_w_lngmo ].GetPtr();
pm.lnSAC = (double (*)[4])aObj[ o_wo_lnsat ].GetPtr();
//    pm.lnSAT = (double *)aObj[ o_wo_lnsat ].GetPtr();
    pm.B     = (double *)aObj[ o_wi_b ].GetPtr();
    pm.U     = (double *)aObj[ o_wo_u ].GetPtr();
    pm.U_r   = (double *)aObj[ o_wd_ur ].GetPtr();
    pm.C     = (double *)aObj[ o_wo_c ].GetPtr();
    pm.IC_m  = (double *)aObj[ o_wd_icm ].GetPtr();
    pm.IC_lm = (double *)aObj[ o_wd_iclcm ].GetPtr();
    pm.IC_wm = (double *)aObj[ o_wd_icwm ].GetPtr();
    pm.BF    = (double *)aObj[ o_wo_bf ].GetPtr();
    pm.XF    = (double *)aObj[ o_w_xf ].GetPtr();
    pm.YF    = (double *)aObj[ o_wo_yf ].GetPtr();
    pm.XFA   = (double *)aObj[ o_ww_xfa ].GetPtr();
    pm.YFA   = (double *)aObj[ o_ww_yfa ].GetPtr();
    pm.Falp  = (double *)aObj[ o_wo_falp ].GetPtr();
    pm.XetaA = (double (*)[MST])aObj[ o_w_xetaa ].GetPtr();
    pm.XetaB = (double (*)[MST])aObj[ o_w_xetab ].GetPtr();
pm.XetaD = (double (*)[MST])aObj[ o_w_xetad ].GetPtr(); // added 12.09.05 KD
    pm.X     = (double *)aObj[ o_w_x ].GetPtr();
    pm.Y     = (double *)aObj[ o_wo_y ].GetPtr();
    pm.XY    = (double *)aObj[ o_w_xy ].GetPtr();
    pm.Qp    = (double *)aObj[ o_wi_qp ].GetPtr();
    pm.Qd    = (double *)aObj[ o_wi_qd ].GetPtr();
    pm.MU    = (double *)aObj[ o_wo_mu ].GetPtr();
    pm.EMU   = (double *)aObj[ o_w_emu ].GetPtr();
    pm.NMU   = (double *)aObj[ o_w_nmu ].GetPtr();
    pm.W     = (double *)aObj[ o_w_w ].GetPtr();
    pm.Fx    = (double *)aObj[ o_wo_fxa ].GetPtr();
    pm.Wx    = (double *)aObj[ o_wo_wx ].GetPtr();
    pm.F     = (double *)aObj[ o_wo_f ].GetPtr();
    pm.F0    = (double *)aObj[ o_w_f0 ].GetPtr();
pm.D     = (double (*)[MST])aObj[ o_w_d ].GetPtr();
//    pm.R     = (double *)aObj[ o_w_r ].GetPtr();
//    pm.R1    = (double *)aObj[ o_w_r1 ].GetPtr();
    pm.sMod  = (char (*)[6])aObj[ o_wi_smod ].GetPtr();
    pm.SB    = (char (*)[MAXICNAME+MAXSYMB])aObj[ o_wd_sb ].GetPtr();
    pm.SB1    = (char (*)[MAXICNAME])aObj[ o_w_sbh ].GetPtr();
    pm.SM    = (char (*)[MAXDCNAME])aObj[ o_wd_sm ].GetPtr();
    pm.SF    = (char (*)[MAXPHNAME+MAXSYMB])aObj[ o_wd_sf ].GetPtr();
    pm.SM2    = (char (*)[MAXDCNAME])aObj[ o_wd_sm2 ].GetPtr();
pm.SM3    = (char (*)[MAXDCNAME])aObj[ o_wd_sm3 ].GetPtr();
    pm.SF2    = (char (*)[MAXPHNAME+MAXSYMB])aObj[ o_wd_sf2 ].GetPtr();
    pm.SFs   = (char (*)[MAXPHNAME+MAXSYMB])aObj[ o_wd_sfs ].GetPtr();
    pm.pbuf  = (char *)aObj[ o_wd_txtb ].GetPtr();
    pm.RLC   = (char *)aObj[ o_wi_rlc ].GetPtr();
    pm.RSC   = (char *)aObj[ o_wi_rsc ].GetPtr();
    pm.RFLC  = (char *)aObj[ o_wi_rflc ].GetPtr();
    pm.RFSC  = (char *)aObj[ o_wi_rfsc ].GetPtr();
    pm.ICC   = (char *)aObj[ o_wi_icc ].GetPtr();
    pm.DCC   = (char *)aObj[ o_wi_dcc ].GetPtr();
pm.DCC3   = (char *)aObj[ o_wi_dcc3 ].GetPtr();
    pm.PHC   = (char *)aObj[ o_wi_phc ].GetPtr();
    pm.SCM   = (char (*)[MST])aObj[ o_wi_scm ].GetPtr();
    pm.SATT  = (char *)aObj[ o_wi_satt ].GetPtr();
    pm.DCCW  = (char *)aObj[ o_wi_dccw ].GetPtr();
//  Added 16.11.2004 by Sveta
    pm.sitXcat = (short *)aObj[ o_wd_sitxcat ].GetPtr();
    pm.sitXan = (short *)aObj[ o_wd_sitxan ].GetPtr();
    pm.sitE = (float *)aObj[ o_wd_site ].GetPtr();

}


// free dynamic memory in objects and values
void TMulti::dyn_kill(int /*q*/)
{
    pm.L1    = (short *)aObj[ o_wi_l1 ].Free();
    pm.LsMod = (short *)aObj[ o_wi_lsmod ].Free();
    pm.LsMdc = (short *)aObj[ o_wi_lsmdc ].Free();
    pm.mui   = (short *)aObj[ o_wi_mui ].Free();
    pm.muk   = (short *)aObj[ o_wi_muk ].Free();
    pm.muj   = (short *)aObj[ o_wi_muj ].Free();
pm.SATX =(short (*)[4])aObj[ o_wi_satndx ].Free();
    pm.DUL   = (double *)aObj[ o_wi_dul ].Free();
    pm.DLL   = (double *)aObj[ o_wi_dll ].Free();
    pm.GEX   = (double *)aObj[ o_wi_gex ].Free();
    pm.PUL   = (double *)aObj[ o_wi_pul ].Free();
    pm.PLL   = (double *)aObj[ o_wi_pll ].Free();
    pm.YOF   = (double *)aObj[ o_wi_yof ].Free();
    pm.PMc   = (float *)aObj[ o_wi_pmc ].Free();
    pm.DMc   = (float *)aObj[ o_wi_dmc ].Free();
    pm.Vol   = (double *)aObj[ o_wio_vol ].Free();
    pm.HYM   = (float *)aObj[ o_w_hym ].Free();
    pm.VL    = (float *)aObj[ o_w_vl ].Free();
    pm.MM    = (double *)aObj[ o_wi_mm ].Free();
    pm.H0    = (float *)aObj[ o_wio_h0 ].Free();
    pm.A0    = (float *)aObj[ o_wio_a0 ].Free();
    pm.U0    = (float *)aObj[ o_wio_u0 ].Free();
    pm.S0    = (float *)aObj[ o_wio_s0 ].Free();
    pm.Cp0   = (float *)aObj[ o_wio_cp0 ].Free();
    pm.Cv0   = (float *)aObj[ o_wio_cv0 ].Free();
    pm.Pparc = (double *)aObj[ o_wi_pparc ].Free();
    pm.Y_m   = (double *)aObj[ o_wd_ym ].Free();
    pm.Y_la  = (double *)aObj[ o_wd_yla ].Free();
    pm.Y_w   = (double *)aObj[ o_wd_yw ].Free();
    pm.Gamma = (double *)aObj[ o_wd_gamma ].Free();
    pm.lnGmf = (double *)aObj[ o_wi_lngmf ].Free();
    pm.lnGmM = (double *)aObj[ o_wio_lngmm ].Free();
    pm.EZ    = (double *)aObj[ o_wi_ez ].Free();
    pm.Wb    = (float *)aObj[ o_wio_wb ].Free();
    pm.Wabs  = (float *)aObj[ o_wio_wabs ].Free();
    pm.Rion  = (float *)aObj[ o_wi_rion ].Free();
    pm.Aalp  = (float *)aObj[ o_wi_aalp ].Free();
    pm.Sigw  = (float *)aObj[ o_wi_sigw ].Free();
    pm.Sigg  = (float *)aObj[ o_wi_sigg ].Free();
    pm.Nfsp  = (float (*)[MST])aObj[ o_wi_nfsp ].Free();
    pm.MASDT = (float (*)[MST])aObj[ o_wi_masdt ].Free();
    pm.FVOL  = (double *)aObj[ o_wd_fvol ].Free();
    pm.FWGT  = (double *)aObj[ o_wd_fwgt ].Free();
    pm.XcapA = (float (*)[MST])aObj[ o_w_xcapa ].Free();
    pm.XcapB = (float (*)[MST])aObj[ o_w_xcapb ].Free();
    pm.XcapD = (float (*)[MST])aObj[ o_w_xcapd ].Free();
    pm.XcapF = (float (*)[MST])aObj[ o_w_xcapf ].Free();
    pm.XdlA  = (float (*)[MST])aObj[ o_wo_xdla ].Free();
    pm.XdlB  = (float (*)[MST])aObj[ o_wo_xdlb ].Free();
    pm.XdlD  = (float (*)[MST])aObj[ o_wo_xdld ].Free();
    pm.XpsiA = (double (*)[MST])aObj[ o_wo_psia ].Free();
    pm.XpsiB = (double (*)[MST])aObj[ o_wo_psib ].Free();
    pm.XpsiD = (double (*)[MST])aObj[ o_wo_psid ].Free();
    pm.Xr0h0 = (float (*)[2])aObj[ o_wi_xr0h0 ].Free();
    pm.XlamA = (float (*)[MST])aObj[ o_wi_xlama ].Free();
    pm.Xetaf = (float (*)[MST])aObj[ o_wi_xetaf ].Free();
    pm.Xcond = (float *)aObj[ o_wi_xcond ].Free();
    pm.Xeps  = (float *)aObj[ o_wi_xeps ].Free();
    pm.Awt   = (float *)aObj[ o_wi_awt ].Free();
    pm.A     = (float *)aObj[ o_wi_a ].Free();
    pm.XFs   = (float *)aObj[ o_wd_xfs ].Free();
    pm.Falps = (float *)aObj[ o_wd_falps ].Free();
    pm.Fug   = (float *)aObj[ o_wd_fug ].Free();
    pm.Fug_l = (float *)aObj[ o_wd_fugl ].Free();
    pm.Ppg_l = (float *)aObj[ o_wd_ppgl ].Free();
    pm.XFTS  = (double (*)[MST])aObj[ o_w_xfts ].Free();
    pm.MASDJ = (float (*)[DFCN])aObj[ o_wo_masdj ].Free();
    pm.G     = (double *)aObj[ o_wo_g ].Free();
    pm.G0    = (double *)aObj[ o_wi_g0 ].Free();
    pm.lnGam = (double *)aObj[ o_wo_lngam ].Free();
    pm.lnGmo = (double *)aObj[ o_w_lngmo ].Free();
//    pm.lnSAT = (double *)aObj[ o_wo_lnsat ].Free();
pm.lnSAC = (double (*)[4])aObj[ o_wo_lnsat ].Free();
    pm.B     = (double *)aObj[ o_wi_b ].Free();
    pm.U     = (double *)aObj[ o_wo_u ].Free();
    pm.U_r   = (double *)aObj[ o_wd_ur ].Free();
    pm.C     = (double *)aObj[ o_wo_c ].Free();
    pm.IC_m  = (double *)aObj[ o_wd_icm ].Free();
    pm.IC_lm = (double *)aObj[ o_wd_iclcm ].Free();
    pm.IC_wm = (double *)aObj[ o_wd_icwm ].Free();
    pm.BF    = (double *)aObj[ o_wo_bf ].Free();
    pm.XF    = (double *)aObj[ o_w_xf ].Free();
    pm.YF    = (double *)aObj[ o_wo_yf ].Free();
    pm.XFA   = (double *)aObj[ o_ww_xfa ].Free();
    pm.YFA   = (double *)aObj[ o_ww_yfa ].Free();
    pm.Falp  = (double *)aObj[ o_wo_falp ].Free();
    pm.XetaA = (double (*)[MST])aObj[ o_w_xetaa ].Free();
    pm.XetaB = (double (*)[MST])aObj[ o_w_xetab ].Free();
pm.XetaD = (double (*)[MST])aObj[ o_w_xetad ].Free();
    pm.X     = (double *)aObj[ o_w_x ].Free();
    pm.Y     = (double *)aObj[ o_wo_y ].Free();
    pm.XY    = (double *)aObj[ o_w_xy ].Free();
    pm.Qp    = (double *)aObj[ o_wi_qp ].Free();
    pm.Qd    = (double *)aObj[ o_wi_qd ].Free();
    pm.MU    = (double *)aObj[ o_wo_mu ].Free();
    pm.EMU   = (double *)aObj[ o_w_emu ].Free();
    pm.NMU   = (double *)aObj[ o_w_nmu ].Free();
    pm.W     = (double *)aObj[ o_w_w ].Free();
    pm.Fx    = (double *)aObj[ o_wo_fxa ].Free();
    pm.Wx    = (double *)aObj[ o_wo_wx ].Free();
    pm.F     = (double *)aObj[ o_wo_f ].Free();
    pm.F0    = (double *)aObj[ o_w_f0 ].Free();
pm.D     = (double (*)[MST])aObj[ o_w_d ].Free();
//    pm.R     = (double *)aObj[ o_w_r ].Free();
//    pm.R1    = (double *)aObj[ o_w_r1 ].Free();
    pm.sMod  = (char (*)[6])aObj[ o_wi_smod ].Free();
    pm.SB    = (char (*)[MAXICNAME+MAXSYMB])aObj[ o_wd_sb ].Free();
    pm.SB1   = (char (*)[MAXICNAME])aObj[ o_w_sbh ].Free();
    pm.SM    = (char (*)[MAXDCNAME])aObj[ o_wd_sm ].Free();
    pm.SF    = (char (*)[MAXPHNAME+MAXSYMB])aObj[ o_wd_sf ].Free();
    pm.SM2    = (char (*)[MAXDCNAME])aObj[ o_wd_sm2 ].Free();
    pm.SM3    = (char (*)[MAXDCNAME])aObj[ o_wd_sm3 ].Free();
    pm.SF2    = (char (*)[MAXPHNAME+MAXSYMB])aObj[ o_wd_sf2 ].Free();
    pm.SFs   = (char (*)[MAXPHNAME+MAXSYMB])aObj[ o_wd_sfs ].Free();
    pm.pbuf  = (char *)aObj[ o_wd_txtb ].Free(); //
    pm.RLC   = (char *)aObj[ o_wi_rlc ].Free();
    pm.RSC   = (char *)aObj[ o_wi_rsc ].Free();
    pm.RFLC  = (char *)aObj[ o_wi_rflc ].Free();
    pm.RFSC  = (char *)aObj[ o_wi_rfsc ].Free();
    pm.ICC   = (char *)aObj[ o_wi_icc ].Free();
    pm.DCC   = (char *)aObj[ o_wi_dcc ].Free();
pm.DCC3   = (char *)aObj[ o_wi_dcc3 ].Free();
    pm.PHC   = (char *)aObj[ o_wi_phc ].Free();
    pm.SCM   = (char (*)[MST])aObj[ o_wi_scm ].Free();
    pm.SATT  = (char *)aObj[ o_wi_satt ].Free();
    pm.DCCW  = (char *)aObj[ o_wi_dccw ].Free();
//  Added 16.11.2004 by Sveta
    pm.sitXcat = (short *)aObj[ o_wd_sitxcat ].Free();
    pm.sitXan = (short *)aObj[ o_wd_sitxan ].Free();
    pm.sitE = (float *)aObj[ o_wd_site ].Free();
}

void TMulti::sit_dyn_new()
{
//  Added 16.11.2004 by Sveta
   if( pm.sitNcat*pm.sitNan )
     pm.sitE = (float *)aObj[ o_wd_site].Alloc(pm.sitNcat, pm.sitNan, F_);
   else
     pm.sitE = (float *)aObj[ o_wd_site].Free();
   if( pm.sitNcat )
     pm.sitXcat = (short *)aObj[ o_wd_sitxcat].Alloc( pm.sitNcat, 1, I_ );
   else
     pm.sitXcat = (short *)aObj[ o_wd_sitxcat].Free();
   if( pm.sitNan )
     pm.sitXan = (short *)aObj[ o_wd_sitxan].Alloc( 1, pm.sitNan, I_ );
   else
     pm.sitXan = (short *)aObj[ o_wd_sitxan].Free();
}

// realloc dynamic memory
void TMulti::dyn_new(int /*q*/)
{
   if( pm.N < 2 || pm.L < 2 || pm.FI < 1 )
        Error( GetName(), "pm.N < 2 || pm.L < 2 || pm.FI < 1" );

    // Part 1

    /* need  always to alloc vectors */
    pm.L1 = (short *)aObj[ o_wi_l1].Alloc( pm.FI, 1, I_ );
    pm.muk = (short *)aObj[ o_wi_muk].Alloc( pm.FI, 1, I_ );
    pm.mui = (short *)aObj[ o_wi_mui].Alloc( pm.N, 1, I_);
    pm.muj = (short *)aObj[ o_wi_muj].Alloc( pm.L, 1, I_ );
    pm.DUL = (double *)aObj[ o_wi_dul].Alloc( pm.L, 1, D_ );
    pm.DLL = (double *)aObj[ o_wi_dll].Alloc( pm.L, 1, D_ );
    pm.Vol = (double *)aObj[ o_wio_vol].Alloc( pm.L, 1, D_ );
    pm.Pparc = (double *)aObj[ o_wi_pparc].Alloc( pm.L, 1, D_ );
    pm.MM = (double *)aObj[ o_wi_mm].Alloc( pm.L, 1, D_ );
    pm.Awt = (float *)aObj[ o_wi_awt].Alloc( pm.N, 1, F_ );
    pm.A = (float *)aObj[ o_wi_a].Alloc( pm.L, pm.N, F_ );
    pm.XFs = (float *)aObj[ o_wd_xfs].Alloc( pm.FI, 1, F_ );
    pm.Falps = (float *)aObj[ o_wd_falps].Alloc( pm.FI, 1, F_ );
    pm.G = (double *)aObj[ o_wo_g].Alloc( pm.L, 1, D_);
    pm.G0 = (double *)aObj[ o_wi_g0].Alloc( pm.L, 1, D_);
    pm.lnGam = (double *)aObj[ o_wo_lngam].Alloc( pm.L, 1, D_);
    pm.lnGmo = (double *)aObj[ o_w_lngmo].Alloc( pm.L, 1, D_);
    pm.B = (double *)aObj[ o_wi_b].Alloc( pm.N, 1, D_ );
    pm.U = (double *)aObj[ o_wo_u].Alloc( pm.N, 1, D_ );
    pm.U_r = (double *)aObj[ o_wd_ur].Alloc( pm.N, 1, D_ );
    pm.C = (double *)aObj[ o_wo_c].Alloc( pm.N, 1, D_ );
    pm.XF = (double *)aObj[ o_w_xf].Alloc( pm.FI, 1, D_ );
    pm.YF = (double *)aObj[ o_wo_yf].Alloc( pm.FI, 1, D_ );
    pm.Falp = (double *)aObj[ o_wo_falp].Alloc( pm.FI, 1, D_ );
    pm.X = (double *)aObj[ o_w_x].Alloc( pm.L, 1, D_);
    pm.Y = (double *)aObj[ o_wo_y].Alloc( pm.L, 1, D_);
    pm.XY = (double *)aObj[ o_w_xy].Alloc( pm.L, 1, D_);
    pm.MU = (double *)aObj[ o_wo_mu].Alloc( pm.L, 1, D_);
    pm.EMU = (double *)aObj[ o_w_emu].Alloc( pm.L, 1, D_);
    pm.NMU = (double *)aObj[ o_w_nmu].Alloc( pm.L, 1, D_);
    pm.W = (double *)aObj[ o_w_w].Alloc( pm.L, 1, D_);
    pm.F = (double *)aObj[ o_wo_f].Alloc( pm.L, 1, D_);
    pm.F0 = (double *)aObj[ o_w_f0].Alloc( pm.L, 1, D_);
    pm.YOF = (double *)aObj[ o_wi_yof].Alloc( pm.FI, 1, D_ );
    pm.SB = (char (*)[MAXICNAME+MAXSYMB])aObj[ o_wd_sb].Alloc(
                pm.N, 1, MAXICNAME+MAXSYMB );
    pm.SB1 = (char (*)[MAXICNAME])aObj[ o_w_sbh].Alloc( 1, pm.N, MAXICNAME );

    // Error in free memory
    //aObj[ o_w_sbh].SetPtr( pm.SB ); aObj[ o_w_sbh ].SetDim( pm.N, 1 );
    pm.SFs = (char (*)[MAXPHNAME+MAXSYMB])aObj[ o_wd_sfs].Alloc(
                 pm.FI, 1, MAXPHNAME+MAXSYMB );
    pm.SM = (char (*)[MAXDCNAME])aObj[ o_wd_sm].Alloc( pm.L, 1, MAXDCNAME );
    pm.SF = (char (*)[MAXPHNAME+MAXSYMB])aObj[ o_wd_sf].Alloc(
                pm.FI, 1, MAXPHNAME+MAXSYMB );
    pm.SM2 = (char (*)[MAXDCNAME])aObj[ o_wd_sm2].Alloc( pm.Ls, 1, MAXDCNAME );
    pm.SF2 = (char (*)[MAXPHNAME+MAXSYMB])aObj[ o_wd_sf2].Alloc(
                 pm.FIs, 1, MAXPHNAME+MAXSYMB );

    pm.RLC = (char *)aObj[ o_wi_rlc].Alloc( pm.L, 1, A_ );
    pm.RSC = (char *)aObj[ o_wi_rsc].Alloc( pm.L, 1, A_ );
    pm.ICC = (char *)aObj[ o_wi_icc].Alloc( pm.N, 1, A_ );
    pm.DCC = (char *)aObj[ o_wi_dcc].Alloc( pm.L, 1, A_ );
    pm.PHC = (char *)aObj[ o_wi_phc].Alloc( pm.FI, 1, A_ );
    pm.DCCW = (char *)aObj[ o_wi_dccw].Alloc( pm.L, 1, A_ );

    pm.lnGmM = (double *)aObj[ o_wio_lngmm].Alloc( pm.L, 1, D_ );
    pm.GEX = (double *)aObj[ o_wi_gex].Alloc( pm.L, 1, D_ );
    pm.FVOL = (double *)aObj[ o_wd_fvol].Alloc( pm.FI, 1, D_ );
    pm.FWGT = (double *)aObj[ o_wd_fwgt].Alloc( pm.FI, 1, D_ );

    if( pm.L > 0 )
    {
        pm.Y_la = (double *)aObj[ o_wd_yla].Alloc( pm.L, 1, D_ );
        pm.Y_w = (double *)aObj[ o_wd_yw].Alloc( pm.L, 1, D_ );
        pm.Fx = (double *)aObj[ o_wo_fxa].Alloc( pm.L, 1, D_);
        pm.Wx = (double *)aObj[ o_wo_wx].Alloc( pm.L, 1, D_);
        pm.VL = (float *)aObj[ o_w_vl].Alloc( pm.L, 1, F_ );
        pm.Gamma = (double *)aObj[ o_wd_gamma].Alloc( pm.L, 1, D_ );  /*!*/
        pm.lnGmf = (double *)aObj[ o_wi_lngmf].Alloc( pm.L, 1, D_ );  /*!*/
//        pm.D = (double *)aObj[ o_w_d].Alloc( pm.L, 1, D_);
    }
    else
    {
        pm.Y_la  = (double *)aObj[ o_wd_yla ].Free();
        pm.Y_w   = (double *)aObj[ o_wd_yw ].Free();
        pm.Fx    = (double *)aObj[ o_wo_fxa ].Free();
        pm.Wx    = (double *)aObj[ o_wo_wx ].Free();
        pm.VL    = (float *)aObj[ o_w_vl ].Free();
        pm.Gamma = (double *)aObj[ o_wd_gamma ].Free();
        pm.lnGmf = (double *)aObj[ o_wi_lngmf ].Free();
//        pm.D     = (double *)aObj[ o_w_d ].Free();
    }

   // Part 2  not always required arrays
    if( pm.FIs > 0 && pm.Ls > 0 )
    {
        pm.BF = (double *)aObj[ o_wo_bf].Alloc( pm.FIs, pm.N, D_ );
        pm.XFA = (double *)aObj[ o_ww_xfa].Alloc( pm.FIs, 1, D_ );
        pm.YFA = (double *)aObj[ o_ww_yfa].Alloc( pm.FIs, 1, D_ );
        pm.LsMod = (short *)aObj[ o_wi_lsmod].Alloc( pm.FIs, 1, I_ );
        pm.LsMdc = (short *)aObj[ o_wi_lsmdc].Alloc( pm.FIs, 1, I_ );
        /* PMs and DMS may be realloced after load arrays LsMod and LsMdc */
        pm.PMc = (float *)aObj[ o_wi_pmc].Alloc( pm.FIs, 1, F_);
        pm.DMc = (float *)aObj[ o_wi_dmc].Alloc( pm.Ls, 1, F_ );
        pm.PUL = (double *)aObj[ o_wi_pul].Alloc( pm.FIs, 1, D_);
        pm.PLL = (double *)aObj[ o_wi_pll].Alloc( pm.FIs, 1, D_);
        //     pm.YOF = (double *)aObj[ o_wi_yof].Alloc( pm.FI, 1, D_ );
        pm.sMod = (char (*)[6])aObj[ o_wi_smod].Alloc( pm.FIs, 1, 6 );
        pm.RFLC = (char *)aObj[ o_wi_rflc].Alloc( pm.FIs, 1, A_ );
        pm.RFSC = (char *)aObj[ o_wi_rfsc].Alloc( pm.FIs, 1, A_ );
    }
    else
    {
        pm.BF    = (double *)aObj[ o_wo_bf ].Free();
        pm.XFA   = (double *)aObj[ o_ww_xfa ].Free();
        pm.YFA   = (double *)aObj[ o_ww_yfa ].Free();
        pm.LsMod = (short *)aObj[ o_wi_lsmod ].Free();
        pm.LsMdc = (short *)aObj[ o_wi_lsmdc ].Free();
        pm.PUL   = (double *)aObj[ o_wi_pul ].Free();
        pm.PLL   = (double *)aObj[ o_wi_pll ].Free();
        //    pm.YOF   = (double *)aObj[ o_wi_yof ].Free();
        pm.PMc   = (float *)aObj[ o_wi_pmc ].Free();
        pm.DMc   = (float *)aObj[ o_wi_dmc ].Free();
        pm.sMod  = (char (*)[6])aObj[ o_wi_smod ].Free();
        pm.RFLC  = (char *)aObj[ o_wi_rflc ].Free();
        pm.RFSC  = (char *)aObj[ o_wi_rfsc ].Free();
    }

    if( pm.LO > 1 )
    {
        pm.Y_m = (double *)aObj[ o_wd_ym].Alloc( pm.L, 1, D_ );
        pm.IC_m = (double *)aObj[ o_wd_icm].Alloc( pm.N, 1, D_ );
        pm.IC_lm = (double *)aObj[ o_wd_iclcm].Alloc( pm.N, 1, D_ );
        pm.IC_wm = (double *)aObj[ o_wd_icwm].Alloc( pm.N, 1, D_ );
    }
    else
    {
        pm.IC_m  = (double *)aObj[ o_wd_icm ].Free();
        pm.IC_lm = (double *)aObj[ o_wd_iclcm ].Free();
        pm.IC_wm = (double *)aObj[ o_wd_icwm ].Free();
        pm.Y_m   = (double *)aObj[ o_wd_ym ].Free();
    }

    /* sorption phases */
    if( syp->PAalp != S_OFF )
    {
        pm.Aalp = (float *)aObj[ o_wi_aalp].Alloc( pm.FI, 1, F_ );
        pm.Xr0h0 = (float (*)[2])aObj[ o_wi_xr0h0].Alloc( pm.FI, 2, F_ );
    }
    else
    {
        pm.Aalp  = (float *)aObj[ o_wi_aalp ].Free();
        pm.Xr0h0 = (float (*)[2])aObj[ o_wi_xr0h0 ].Free();
    }

   if( syp->PSigm != S_OFF )
           pm.Sigw = (float *)aObj[ o_wi_sigw].Alloc( pm.FI, 1, F_ );
    else  pm.Sigw  = (float *)aObj[ o_wi_sigw ].Free();

    if( syp->PSigm != S_OFF )
        pm.Sigg = (float *)aObj[ o_wi_sigg].Alloc( pm.FI, 1, F_ );
    else  pm.Sigg  = (float *)aObj[ o_wi_sigg ].Free();

    if( pm.E )
    {
        pm.EZ = (double *)aObj[ o_wi_ez].Alloc( pm.L, 1, D_ );
        pm.Xcond = (float *)aObj[ o_wi_xcond].Alloc( pm.FI, 1, F_ );
        pm.Xeps = (float *)aObj[ o_wi_xeps].Alloc( pm.FI, 1, F_ );
    }
    else
    {
        pm.EZ    = (double *)aObj[ o_wi_ez ].Free();
        pm.Xcond = (float *)aObj[ o_wi_xcond ].Free();
        pm.Xeps  = (float *)aObj[ o_wi_xeps ].Free();
    }

    if( pm.FIat > 0 && /*pm.Lads > 0 &&Sveta 12/09/99*/ pm.FIs > 0 )
    { /* ADSORPTION AND ION EXCHANGE */
pm.SATX = (short (*)[4])aObj[ o_wi_satndx].Alloc( pm.Lads, 4, I_ );
        pm.SCM  = (char (*)[MST])aObj[ o_wi_scm].Alloc( pm.FIs, pm.FIat, A_ );
        //   memset( pm.SCM, SC_NOT_USED, pm.FIs*pm.FIat );
        pm.Nfsp = (float (*)[MST])aObj[ o_wi_nfsp].Alloc( pm.FIs, pm.FIat, F_ );
        pm.MASDT = (float (*)[MST])aObj[ o_wi_masdt].Alloc( pm.FIs, pm.FIat, F_ );
        pm.XcapA = (float (*)[MST])aObj[ o_w_xcapa].Alloc( pm.FIs, pm.FIat, F_ );
        pm.XcapB = (float (*)[MST])aObj[ o_w_xcapb].Alloc( pm.FIs, pm.FIat, F_ );
        pm.XcapD = (float (*)[MST])aObj[ o_w_xcapd].Alloc( pm.FIs, pm.FIat, F_ );
        pm.XcapF = (float (*)[MST])aObj[ o_w_xcapf].Alloc( pm.FIs, pm.FIat, F_ );
        pm.XdlA = (float (*)[MST])aObj[ o_wo_xdla].Alloc( pm.FIs, pm.FIat, F_ );
        pm.XdlB = (float (*)[MST])aObj[ o_wo_xdlb].Alloc( pm.FIs, pm.FIat, F_ );
        pm.XdlD = (float (*)[MST])aObj[ o_wo_xdld].Alloc( pm.FIs, pm.FIat, F_ );
        pm.XpsiA = (double (*)[MST])aObj[ o_wo_psia].Alloc( pm.FIs, pm.FIat, D_ );
        pm.XpsiB = (double (*)[MST])aObj[ o_wo_psib].Alloc( pm.FIs, pm.FIat, D_ );
        pm.XpsiD = (double (*)[MST])aObj[ o_wo_psid].Alloc( pm.FIs, pm.FIat, D_ );
        pm.XlamA = (float (*)[MST])aObj[ o_wi_xlama].Alloc( pm.FIs, pm.FIat, F_ );
        pm.Xetaf = (float (*)[MST])aObj[ o_wi_xetaf].Alloc( pm.FIs, pm.FIat, F_ );
        pm.XetaA = (double (*)[MST])aObj[ o_w_xetaa].Alloc( pm.FIs, pm.FIat, D_ );
        pm.XetaB = (double (*)[MST])aObj[ o_w_xetab].Alloc( pm.FIs, pm.FIat, D_ );
pm.XetaD = (double (*)[MST])aObj[ o_w_xetad].Alloc( pm.FIs, pm.FIat, D_ );
        pm.XFTS = (double (*)[MST])aObj[ o_w_xfts].Alloc( pm.FIs, pm.FIat, D_ );
pm.SATT = (char *)aObj[ o_wi_satt].Alloc( pm.Lads, 1, A_ );
pm.MASDJ = (float (*)[DFCN])aObj[ o_wo_masdj].Alloc( pm.Lads, DFCN, F_ );
//        pm.MASDJ = (float *)aObj[ o_wo_masdj].Alloc( pm.Ls, 1, F_ );
pm.D = (double (*)[MST])aObj[ o_w_d].Alloc( MST, MST, D_);
pm.lnSAC = (double (*)[4])aObj[ o_wo_lnsat].Alloc( pm.Lads, 4, D_ );
pm.SM3 = (char (*)[MAXDCNAME])aObj[ o_wd_sm3].Alloc( pm.Lads, 1, MAXDCNAME );
pm.DCC3 = (char *)aObj[ o_wi_dcc3].Alloc( pm.Lads, 1, A_ );
    }
    else
    {
        pm.SATX =(short (*)[4])aObj[ o_wi_satndx ].Free();
        pm.SCM   = (char (*)[MST])aObj[ o_wi_scm ].Free();
        pm.Nfsp  = (float (*)[MST])aObj[ o_wi_nfsp ].Free();
        pm.MASDT = (float (*)[MST])aObj[ o_wi_masdt ].Free();
        pm.XcapA = (float (*)[MST])aObj[ o_w_xcapa ].Free();
        pm.XcapB = (float (*)[MST])aObj[ o_w_xcapb ].Free();
        pm.XcapD = (float (*)[MST])aObj[ o_w_xcapd ].Free();
        pm.XcapF = (float (*)[MST])aObj[ o_w_xcapf ].Free();
        pm.XdlA  = (float (*)[MST])aObj[ o_wo_xdla ].Free();
        pm.XdlB  = (float (*)[MST])aObj[ o_wo_xdlb ].Free();
        pm.XdlD  = (float (*)[MST])aObj[ o_wo_xdld ].Free();
        pm.XpsiA = (double (*)[MST])aObj[ o_wo_psia ].Free();
        pm.XpsiB = (double (*)[MST])aObj[ o_wo_psib ].Free();
        pm.XpsiD = (double (*)[MST])aObj[ o_wo_psid ].Free();
        pm.XlamA = (float (*)[MST])aObj[ o_wi_xlama ].Free();
        pm.Xetaf = (float (*)[MST])aObj[ o_wi_xetaf ].Free();
        pm.XFTS  = (double (*)[MST])aObj[ o_w_xfts ].Free();
        pm.XetaA = (double (*)[MST])aObj[ o_w_xetaa ].Free();
        pm.XetaB = (double (*)[MST])aObj[ o_w_xetab ].Free();
pm.XetaD = (double (*)[MST])aObj[ o_w_xetad ].Free();  // added 12.09.05  KD
        pm.SATT  = (char *)aObj[ o_wi_satt ].Free();
        pm.MASDJ = (float (*)[DFCN])aObj[ o_wo_masdj ].Free();
pm.lnSAC = (double (*)[4])aObj[ o_wo_lnsat ].Free();
pm.D = (double (*)[MST])aObj[ o_w_d ].Free();
pm.SM3 = (char (*)[MAXDCNAME])aObj[ o_wd_sm3 ].Free();
pm.DCC3 = (char *)aObj[ o_wi_dcc3 ].Free();
    }

    if( pm.PG > 0 )
    {
        pm.Fug = (float *)aObj[ o_wd_fug].Alloc( pm.PG, 1, F_ );
        pm.Fug_l = (float *)aObj[ o_wd_fugl].Alloc( pm.PG, 1, F_ );
        pm.Ppg_l = (float *)aObj[ o_wd_ppgl].Alloc( pm.PG, 1, F_ );
    }
    else
    {
        pm.Fug   = (float *)aObj[ o_wd_fug ].Free();
        pm.Fug_l = (float *)aObj[ o_wd_fugl ].Free();
        pm.Ppg_l = (float *)aObj[ o_wd_ppgl ].Free();
    }

    // Part 3

    if( pm.Ls > 1 && pm.FIs > 0 )
    {
        pm.Wb = (float *)aObj[ o_wio_wb].Alloc( pm.Ls, 1, F_ );
        pm.Wabs = (float *)aObj[ o_wio_wabs].Alloc( pm.Ls, 1, F_ );
        pm.Rion = (float *)aObj[ o_wi_rion].Alloc( pm.Ls, 1, F_ );
        pm.Qp = (double *)aObj[ o_wi_qp].Alloc( pm.FIs, QPSIZE, D_ );
        pm.Qd = (double *)aObj[ o_wi_qd].Alloc( pm.FIs, QDSIZE, D_ );
    }
    else
    {
        pm.Wb    = (float *)aObj[ o_wio_wb ].Free();
        pm.Wabs  = (float *)aObj[ o_wio_wabs ].Free();
        pm.Rion  = (float *)aObj[ o_wi_rion ].Free();
        pm.Qp    = (double *)aObj[ o_wi_qp ].Free();
        pm.Qd    = (double *)aObj[ o_wi_qd ].Free();
    }

    /* pm.R = (float *)aObj[ o_w_r].Alloc( pm.N, pm.N+1, D_ ); */
    sit_dyn_new();
}


//set default information
void TMulti::set_def( int /*q*/)
{
    memcpy( &pm.PunE, "jjbC", 4 );
    memset( &pm.N, 0, 36*sizeof(short));
    memset( &pm.TC, 0, 54*sizeof(double));
    pm.sitNcat = 0;
    pm.sitNan = 0;
    pm.L1    = 0;
    pm.LsMod = 0;
    pm.LsMdc = 0;
    pm.mui   = 0;
    pm.muk   = 0;
    pm.muj   = 0;
    pm.SATX =0;
    pm.DUL   = 0;
    pm.DLL   = 0;
    pm.GEX   = 0;
    pm.PUL   = 0;
    pm.PLL   = 0;
    pm.YOF   = 0;
    pm.PMc   = 0;
    pm.DMc   = 0;
    pm.Vol   = 0;
    pm.HYM   = 0;
    pm.VL    = 0;
    pm.MM    = 0;
   pm.H0    = 0;
    pm.A0    = 0;
    pm.U0    = 0;
    pm.S0    = 0;
    pm.Cp0   = 0;
    pm.Cv0   = 0;
    pm.Pparc = 0;
    pm.Y_m   = 0;
    pm.Y_la  = 0;
    pm.Y_w   = 0;
    pm.Gamma = 0;
    pm.lnGmf = 0;
    pm.lnGmM = 0;
    pm.EZ    = 0;
    pm.Wb    = 0;
    pm.Wabs  = 0;
    pm.Rion  = 0;
    pm.Aalp  = 0;
    pm.Sigw  = 0;
    pm.Sigg  = 0;
    pm.Nfsp  = 0;
    pm.MASDT = 0;
    pm.FVOL  = 0;
    pm.FWGT  = 0;
    pm.XcapA = 0;
    pm.XcapB = 0;
    pm.XcapD = 0;
    pm.XdlA  = 0;
    pm.XdlB  = 0;
    pm.XdlD  = 0;
    pm.XpsiA = 0;
    pm.XpsiB = 0;
    pm.XpsiD = 0;
    pm.Xr0h0 = 0;
    pm.XlamA = 0;
    pm.Xetaf = 0;
    pm.Xcond = 0;
    pm.Xeps  = 0;
    pm.Awt   = 0;
    pm.A     = 0;
    pm.XFs   = 0;
        pm.Falps = 0;
        pm.Fug   = 0;
        pm.Fug_l = 0;
        pm.Ppg_l = 0;
        pm.XFTS  = 0;
        pm.MASDJ = 0;
        pm.G     = 0;
        pm.G0    = 0;
        pm.lnGam = 0;
        pm.lnGmo = 0;
//        pm.lnSAT = 0;
        pm.lnSAC = 0;
        pm.B     = 0;
        pm.U     = 0;
        pm.U_r   = 0;
        pm.C     = 0;
        pm.IC_m  = 0;
        pm.IC_lm = 0;
        pm.IC_wm = 0;
        pm.BF    = 0;
        pm.XF    = 0;
        pm.YF    = 0;
        pm.XFA   = 0;
        pm.YFA   = 0;
        pm.Falp  = 0;
        pm.XetaA = 0;
        pm.XetaB = 0;
pm.XetaD = 0;         // added 12.09.05  KD
        pm.X     = 0;
        pm.Y     = 0;
        pm.XY    = 0;
        pm.Qp    = 0;
        pm.Qd    = 0;
        pm.MU    = 0;
        pm.EMU   = 0;
        pm.NMU   = 0;
        pm.W     = 0;
        pm.Fx    = 0;
        pm.Wx    = 0;
        pm.F     = 0;
        pm.F0    = 0;
        pm.D     = 0;
     //   pm.R     = 0;
     //   pm.R1    = 0;
        pm.sMod  = 0;
        pm.SB    = 0;
        pm.SB1    = 0; // added Sveta 6/07/2001
        pm.SM    = 0;
        pm.SF    = 0;
        pm.SFs   = 0;
        pm.pbuf  = 0; //
        pm.RLC   = 0;
        pm.RSC   = 0;
        pm.RFLC  = 0;
        pm.RFSC  = 0;
        pm.ICC   = 0;
        pm.DCC   = 0;
        pm.PHC   = 0;
        pm.SCM   = 0;
        pm.SATT  = 0;
        pm.DCCW  = 0;
        pm.XcapF = 0;
        pm.SM2    = 0;
        pm.SM3    = 0;
        pm.SF2    = 0;
        pm.DCC3   = 0;
//  Added 16.11.2004 by Sveta
    pm.sitXcat = 0;
    pm.sitXan = 0;
    pm.sitE = 0;
 }

//--------------------- End of ms_multi.cpp ---------------------------
