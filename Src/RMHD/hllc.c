/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Implement the HLLC Riemann solver for relativistic MHD.

  Solve the Riemann problem for the relativistic MHD (RMHD) equations 
  using the HLLC solver of Mignone & Bodo (2006).
   
  On input, it takes left and right primitive sweep vectors 
  \c sweep->vL and \c sweep->vR at zone edge \c i+1/2;
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
    - "An HLLC Riemann solver for relativistic flows - II. 
       Magnetohydrodynamics", Mignone and Bodo, MNRAS (2006) 368, 1040.

  \authors A. Mignone (mignone@ph.unito.it)
  \date    June 10, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"
#define BX_MIN  1.e-6

/* ********************************************************************* */
void HLLC_Solver (const Sweep *sweep, int beg, int end, 
                  double *cmax, Grid *grid)
/*!
 * Solve the RMHD Riemann problem using the HLLC Riemann solver.
 *
 * \param[in,out] sweep   pointer to Sweep structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int    nv, i;
  double scrh;
  double usL[NFLX], usR[NFLX];
  double Bx, Bys, Bzs;
  double BtFBt, Bt2, FBt2;
  double a, b, c;
  double ps, vxs, vys, vzs, gammas_2, vBs;
  double vxl, vxr, alpha_l, alpha_r;

  double *uL, *uR, *SL, *SR;
  static double **Uhll, **Fhll;
  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);
  double **fL = stateL->flux, **fR = stateR->flux;
  double *a2L = stateL->a2,   *a2R = stateR->a2;
  double  *hL = stateL->h,     *hR = stateR->h;
  double  *pL = stateL->prs,   *pR = stateR->prs;

  if (Uhll == NULL){
    Uhll = ARRAY_2D(NMAX_POINT, NVAR, double);
    Fhll = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

#ifdef GLM_MHD
  GLM_Solve (sweep, beg, end, grid);
#endif

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */
   
  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

  Flux (stateL, beg, end);
  Flux (stateR, beg, end);

/* ----------------------------------------
     get max and min signal velocities
   ---------------------------------------- */

  SL = sweep->SL; SR = sweep->SR;
  HLL_Speed (stateL, stateR, SL, SR, beg, end);

/* ----------------------------------------------
           compute HLL sweep and flux
   ---------------------------------------------- */
   
  for (i = beg; i <= end; i++) {

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;

    uL = stateL->u[i];
    uR = stateR->u[i];
    scrh = 1.0/(SR[i] - SL[i]);
    for (nv = NFLX; nv--;  ){  
      Uhll[i][nv] =   SR[i]*uR[nv] - SL[i]*uL[nv] 
                    + fL[i][nv] - fR[i][nv];
      Uhll[i][nv] *= scrh;

      Fhll[i][nv]  =   SL[i]*SR[i]*(uR[nv] - uL[nv])
                     + SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
      Fhll[i][nv] *= scrh;
    }
    Uhll[i][MXn] += (pL[i] - pR[i])*scrh;
    Fhll[i][MXn] += (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
  }

/* ----------------------------------------------
           Solve Riemann problem
   ---------------------------------------------- */
  
  for (i = beg; i <= end; i++) {

    if (SL[i] >= 0.0){
      for (nv = NFLX; nv--; ){
        sweep->flux[i][nv] = fL[i][nv];
      }
      sweep->press[i] = pL[i];
    }else if (SR[i] <= 0.0){
      for (nv = NFLX; nv--; ){
        sweep->flux[i][nv] = fR[i][nv];
      }
      sweep->press[i] = pR[i];
    }else{

      uL = stateL->u[i];
      uR = stateR->u[i];
#if SHOCK_FLATTENING == MULTID
      if ((sweep->flag[i] & FLAG_HLL) || (sweep->flag[i+1] & FLAG_HLL)){
        scrh = 1.0/(SR[i] - SL[i]);
        for (nv = NFLX; nv--; ){
          sweep->flux[i][nv]  = SL[i]*SR[i]*(uR[nv] - uL[nv])
                             +  SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
          sweep->flux[i][nv] *= scrh;
        }
        sweep->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
        continue;
      }
#endif

      vxl = stateL->v[i][VXn];
      vxr = stateR->v[i][VXn];
      
      EXPAND(Bx  = Uhll[i][BXn];  , 
             Bys = Uhll[i][BXt];  ,
             Bzs = Uhll[i][BXb];)

      #if RMHD_REDUCED_ENERGY == YES
       Uhll[i][ENG] += Uhll[i][RHO];
       Fhll[i][ENG] += Fhll[i][RHO];
      #endif    

      if (fabs(Bx) < BX_MIN) {

        BtFBt = Bt2 = FBt2 = 0.0;
        a  = Fhll[i][ENG];
        b  = - (Fhll[i][MXn] + Uhll[i][ENG]);
        c  = Uhll[i][MXn];

      }else{

        BtFBt = EXPAND(0.0, + Uhll[i][BXt]*Fhll[i][BXt], 
                            + Uhll[i][BXb]*Fhll[i][BXb]);
                 
        Bt2 = EXPAND(0.0, + Uhll[i][BXt]*Uhll[i][BXt],  
                          + Uhll[i][BXb]*Uhll[i][BXb]);
                
        FBt2 = EXPAND(0.0, + Fhll[i][BXt]*Fhll[i][BXt],
                           + Fhll[i][BXb]*Fhll[i][BXb]);                
                  
        a  = Fhll[i][ENG] - BtFBt;
        b  = Bt2 + FBt2 - (Fhll[i][MXn] + Uhll[i][ENG]);
        c  = Uhll[i][MXn] - BtFBt;
      }

      /* -------------------------
          solve quadratic, get v* 
         ------------------------- */
/*
      if (fabs(a) > 1.e-12){       
        vxs = 0.5*(- b - sqrt(b*b - 4.0*a*c))/a; 
      }else{
        vxs = -c/b;
      }
*/

/*
      scrh = -0.5*(b + DSIGN(b)*sqrt(b*b - 4.0*a*c));
      vxs  = c/scrh;
*/

      scrh = 1.0 + sqrt(1.0 - 4.0*a*c/(b*b));
      vxs  = - 2.0*c/(b*scrh);

      if (fabs(Bx) < BX_MIN) {

      /* -------------------------------
           the value of vy and vz
               is irrelevant in this case  
         ------------------------------- */
 
        ps  = Fhll[i][MXn] - Fhll[i][ENG]*vxs;              

        alpha_l = (SL[i] - vxl)/(SL[i] - vxs);
        alpha_r = (SR[i] - vxr)/(SR[i] - vxs);

        usL[RHO] = uL[RHO]*alpha_l;
        usR[RHO] = uR[RHO]*alpha_r;
      
        usL[ENG] = (SL[i]*uL[ENG] - fL[i][ENG] + ps*vxs)/(SL[i] - vxs);
        usR[ENG] = (SR[i]*uR[ENG] - fR[i][ENG] + ps*vxs)/(SR[i] - vxs);

        EXPAND( usL[MXn] = (usL[ENG] + ps)*vxs; 
                usR[MXn] = (usR[ENG] + ps)*vxs;  ,
                usL[MXt] = uL[MXt]*alpha_l; 
                usR[MXt] = uR[MXt]*alpha_r;      ,  
                usL[MXb] = uL[MXb]*alpha_l; 
                usR[MXb] = uR[MXb]*alpha_r; )

        #if RMHD_REDUCED_ENERGY == YES
         usL[MXn] += usL[RHO]*vxs;
         usR[MXn] += usR[RHO]*vxs;
        #endif

        EXPAND( usL[BXn] = Bx;
                usR[BXn] = Bx;                   ,
                usL[BXt] = uL[BXt]*alpha_l;
                usR[BXt] = uR[BXt]*alpha_r;   ,  
                usL[BXb] = uL[BXb]*alpha_l;
                usR[BXb] = uR[BXb]*alpha_r; )

      }else{

        EXPAND(                                     ,       
                vys = (Bys*vxs - Fhll[i][BXt])/Bx;   ,
                vzs = (Bzs*vxs - Fhll[i][BXb])/Bx; )

        ps = (BtFBt - Fhll[i][ENG])*vxs 
            + Bx*Bx + Fhll[i][MXn] - FBt2;

        gammas_2 = EXPAND(vxs*vxs, + vys*vys, + vzs*vzs);
        gammas_2 = 1.0 - gammas_2;
        vBs = EXPAND(vxs*Bx, + vys*Bys, + vzs*Bzs);

        alpha_l = (SL[i] - vxl)/(SL[i] - vxs);
        alpha_r = (SR[i] - vxr)/(SR[i] - vxs);

        usL[RHO] = uL[RHO]*alpha_l;
        usR[RHO] = uR[RHO]*alpha_r;

        usL[ENG] = (SL[i]*uL[ENG] - fL[i][ENG] + ps*vxs - vBs*Bx)/(SL[i] - vxs);
        usR[ENG] = (SR[i]*uR[ENG] - fR[i][ENG] + ps*vxs - vBs*Bx)/(SR[i] - vxs);

        EXPAND( usL[MXn] = (usL[ENG] + ps)*vxs - vBs*Bx; 
                usR[MXn] = (usR[ENG] + ps)*vxs - vBs*Bx;                                                     ,
                usL[MXt] = (SL[i]*uL[MXt] - fL[i][MXt] - Bx*(Bys*gammas_2 + vBs*vys))/(SL[i] - vxs); 
                usR[MXt] = (SR[i]*uR[MXt] - fR[i][MXt] - Bx*(Bys*gammas_2 + vBs*vys))/(SR[i] - vxs);   ,  
                usL[MXb] = (SL[i]*uL[MXb] - fL[i][MXb] - Bx*(Bzs*gammas_2 + vBs*vzs))/(SL[i] - vxs); 
                usR[MXb] = (SR[i]*uR[MXb] - fR[i][MXb] - Bx*(Bzs*gammas_2 + vBs*vzs))/(SR[i] - vxs); )

        #if RMHD_REDUCED_ENERGY == YES
         usL[MXn] += usL[RHO]*vxs;
         usR[MXn] += usR[RHO]*vxs;
        #endif

        EXPAND( usL[BXn] = usR[BXn] = Bx;    ,  
                usL[BXt] = usR[BXt] = Bys;   ,
                usL[BXb] = usR[BXb] = Bzs; )
      }                
      
      #ifdef GLM_MHD
       usL[PSI_GLM] = usR[PSI_GLM] = uL[PSI_GLM];
      #endif

/* -------------------------------------------
       Check consistency condition 
   ------------------------------------------- */
/*
      for (nv = 0; nv < NVAR; nv++) {
        scrh  = (vxs - SL[i])*usL[nv] + (SR[i] - vxs)*usR[nv];
        scrh -= SR[i]*uR[i][nv] - SL[i]*uL[i][nv] +
                fL[i][nv] - fR[i][nv];
        if (fabs(scrh/(SR[i]-SL[i])) > 1.e-5){
          printf (" Consistency condition violated\n");
          printf ("%d %d  %12.6e\n",i,nv, scrh);
        }
      }
*/
/*  ----  Compute HLLC flux  ----  */

      if (vxs > 0.0) {
        for (nv = NFLX; nv--;   ) {
          sweep->flux[i][nv] = fL[i][nv] + SL[i]*(usL[nv] - uL[nv]);
        }
        sweep->press[i] = pL[i];
      }else {
        for (nv = NFLX; nv--; ) {
          sweep->flux[i][nv] = fR[i][nv] + SR[i]*(usR[nv] - uR[nv]);
        }
        sweep->press[i] = pR[i];
      }
    }

  }   /* -- end loop on points -- */

/* --------------------------------------------------------
              initialize souRce term
   -------------------------------------------------------- */
 
  #if DIVB_CONTROL == EIGHT_WAVES
/*
   POWELL_DIVB_SOURCE (sweep, beg, end, grid);
*/
  /* ----------------------------------------------------
       to avoid conversion problems in HLL_DIVB_SOURCE, 
       we use the HLL average provided by SR = -SL = 1 
     ---------------------------------------------------- */
/*            
   for (i = beg; i <= end; i++) {
     uL = sweep->uL[i]; uR = sweep->uR[i];
     for (nv = 0; nv < NFLX; nv++) {
       Uhll[i][nv] = 0.5*(uR[nv] + uL[nv] + fL[i][nv] - fR[i][nv]);
     }
     Uhll[i][MXn] += (pL[i] - pR[i])*0.5;
     for (nv = NFLX; nv < NVAR; nv++) Uhll[i][nv] = 0.0;
   }
*/
   HLL_DIVB_SOURCE (sweep, Uhll, beg+1, end, grid);
  #endif

}
#undef BX_MIN
