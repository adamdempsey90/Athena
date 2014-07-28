#include "copyright.h"
/*============================================================================*/
/*! \file planet-disk.c
 *  \brief Problem generator for planet embedded in a disk, using the
 *   shearing sheet approximation.
 *
 * Code must be configured using --enable-shearing-box */
/*============================================================================*/
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * PlanetPot()   - static gravitational potential of planet
 * UnstratifiedDisk() - tidal potential in shearing box
 * expr_dV2()    - computes delta(Vy)
 * expr_Lflux()	 - computes ang momentum flux rho vx vy
 * hst_*         - new history variables
 *============================================================================*/


void constant_iib(GridS *pGrid);
void constant_oib(GridS *pGrid);
static void flux_output_func(MeshS *pM, OutputS *pOut);

static Real Mp,Rsoft,Xplanet,Yplanet,Zplanet;
static Real ramp_time,insert_time;
static Real PlanetPot(const Real x1, const Real x2, const Real x3);
static Real dx2PlanetPot(const Real x1, const Real x2, const Real x3);
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static Real expr_dV2(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j,const int k);
static Real hst_rho_dVy2(const GridS *pG,const int i, const int j, const int k);
static Real expr_Lflux(const GridS *pG,const int i,const int j,const int k);
static Real BesselK(const int order,const Real x);
static Real bessi0(const Real x);
#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k);
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,BCFlag;
  Real x1,x2,x3;
  Real den = 1.0, pres = 1.0e-6;
  static int frst=1;  /* flag so new history variables enrolled only once */

#ifdef SHEARING_BOX
/* specify xy (r-phi) plane */
  ShBoxCoord = xy;
#endif

/* Read problem parameters.  Note Omega_0 set to 10^{-3} by default */
#ifdef SHEARING_BOX
  Omega_0 = par_getd_def("problem","omega",1.0e-3);
  qshear  = par_getd_def("problem","qshear",1.5);
#endif
  Mp      = par_getd_def("problem","Mplanet",0.0);
  Xplanet = par_getd_def("problem","Xplanet",0.0);
  Yplanet = par_getd_def("problem","Yplanet",0.0);
  Zplanet = par_getd_def("problem","Zplanet",0.0);
  Rsoft   = par_getd_def("problem","Rsoft",0.1);
  ramp_time = 0.0;
  insert_time = par_getd_def("problem","insert_time",0.0);

/* Compute field strength based on beta.  */
#ifdef ISOTHERMAL
  pres = Iso_csound2;
#endif

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

/* Initialize d, M, and P.  With FARGO do not initialize the background shear */

      pGrid->U[k][j][i].d  = den;
//      pGrid->U[k][j][i].d = 1.0+.5*(1-.02)*(tanh((x1-10.0)/3.5)-tanh((x1+10.0)/3.5))-.5*1.1*(tanh((x1-10.0)/15.0)-tanh((x1+10.0)/15.0));
      pGrid->U[k][j][i].M1 = 0.0;
      pGrid->U[k][j][i].M2 = 0.0;
#ifdef SHEARING_BOX
#ifndef FARGO
      pGrid->U[k][j][i].M2 -= den*(qshear*Omega_0*x1);
#endif
#endif
      pGrid->U[k][j][i].M3 = 0.0;
#ifdef ADIABATIC
      pGrid->U[k][j][i].E = pres/Gamma_1
        + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) 
             + SQR(pGrid->U[k][j][i].M3))/den;
#endif

    }
  }}

/* enroll gravitational potential of planet & shearing-box potential fns */

  StaticGravPot = PlanetPot;
  ShearingBoxPot = UnstratifiedDisk;

/* enroll new history variables, only once  */

  if (frst == 1) {
    dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
    dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
#ifdef ADIABATIC
    dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif
    frst = 0;
  }

/* With viscosity and/or resistivity, read diffusion coeffs */
#ifdef VISCOSITY
  nu_iso = par_getd_def("problem","nu_iso",0.0);
  nu_aniso = par_getd_def("problem","nu_aniso",0.0);
#endif

/* Enroll outflow BCs if perdiodic BCs NOT selected.  This assumes the root
 * level grid is specified by the <domain1> block in the input file */

  BCFlag = par_geti_def("domain1","bc_ix1",0);
  if (BCFlag != 4) {
    if (pDomain->Disp[0] == 0) bvals_mhd_fun(pDomain, left_x1,  constant_iib);
  }
  BCFlag = par_geti_def("domain1","bc_ox1",0);
  if (BCFlag != 4) {
    if (pDomain->MaxX[0] == pDomain->RootMaxX[0])
      bvals_mhd_fun(pDomain, right_x1, constant_oib);
  }

  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

/* 'problem_read_restart' must enroll gravity on restarts */

void problem_read_restart(MeshS *pM, FILE *fp)
{
  int nl,nd,BCFlag_ix1,BCFlag_ox1;
/* Read Omega, and with viscosity and/or resistivity, read eta_Ohm and nu_V */

#ifdef SHEARING_BOX
  Omega_0 = par_getd_def("problem","omega",1.0e-3);
  qshear  = par_getd_def("problem","qshear",1.5);
#endif
  Mp      = par_getd_def("problem","Mplanet",0.0);
  Xplanet = par_getd_def("problem","Xplanet",0.0);
  Yplanet = par_getd_def("problem","Yplanet",0.0);
  Zplanet = par_getd_def("problem","Zplanet",0.0);
  Rsoft   = par_getd_def("problem","Rsoft",0.1);
  ramp_time = 0.0;
  insert_time = par_getd_def("problem","insert_time",0.0);
#ifdef VISCOSITY
  nu_iso = par_getd_def("problem","nu_iso",0.0);
  nu_aniso = par_getd_def("problem","nu_aniso",0.0);
#endif

/* enroll gravitational potential of planet & shearing-box potential fns */

  StaticGravPot = PlanetPot;
  ShearingBoxPot = UnstratifiedDisk;

/* enroll new history variables */

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
  dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
#ifdef ADIABATIC
  dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif

  BCFlag_ix1 = par_geti_def("domain1","bc_ix1",0);
  BCFlag_ox1 = par_geti_def("domain1","bc_ox1",0);
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Disp[0] == 0 && BCFlag_ix1 != 4) 
        bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x1,  constant_iib);
      if (pM->Domain[nl][nd].MaxX[0] == pM->Domain[nl][nd].RootMaxX[0] 
          && BCFlag_ox1 != 4)
        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x1, constant_oib);
    }
  }

  return;
}

/* Get_user_expression computes dVy */
ConsFun_t get_usr_expr(const char *expr)
{
   if(strcmp(expr,"dVy")==0) return expr_dV2;
   return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  if(strcmp(name,"flux")==0) return flux_output_func;
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  ramp_time = pM->time;
}

void Userwork_after_loop(MeshS *pM)
{
}

/*------------------------------------------------------------------------------
 * PlanetPot:
 */
/*! \fn static Real PlanetPot(const Real x1, const Real x2, const Real x3)
 *  \brief static gravitational potential of planet */
static Real PlanetPot(const Real x1, const Real x2, const Real x3)
{
  Real rad,phi=0.0;
  rad = sqrt(SQR(x1-Xplanet) + SQR(x2-Yplanet) + SQR(x3-Zplanet));
  phi = -1.0*MIN(1.0,(ramp_time/(insert_time+0.0001)))*Mp/(sqrt(rad*rad+Rsoft*Rsoft));
	
//	double kval = 1.0;
//	kval *= 2*M_PI/40.0;

// 	 rad = sqrt(Rsoft*Rsoft+SQR(x1-Xplanet));
  
// 	phi = -1.0*MIN(1.0,(ramp_time/(insert_time+0.0001)))*Mp*cos(x1-Xplanet);
  
//   phi = -MIN(1.0,(ramp_time/(insert_time+0.0001)))*(Mp/M_PI)* 
//			BesselK(0,fabs(kval*rad)) * cos(kval*(x2-Yplanet));

//	phi = -1.0*MIN(1.0,(ramp_time/(insert_time+0.0001)))*Mp *
//			cos(Rsoft+sqrt(SQR(x2-Yplanet))) * exp(-sqrt(SQRT(x1-Xplanet)));


  return phi;
}
static Real dx2PlanetPot(const Real x1, const Real x2, const Real x3)
{
  Real rad,phi=0.0;
  double kval = 1.0;
//	kval *= 2*M_PI/40.0;
//  rad = sqrt(SQR(x1-Xplanet) + SQR(x2-Yplanet) + SQR(x3-Zplanet));
//  phi = MIN(1.0,(ramp_time/(insert_time+0.0001)))*Mp*x2/pow(rad+Rsoft,3);

   rad = sqrt(Rsoft*Rsoft+SQR(x1-Xplanet));
  
//  phi=MIN(1.0,(ramp_time/(insert_time+0.0001)))*Mp*sin(x1-Xplanet);
  
   phi = 2.0*MIN(1.0,(ramp_time/(insert_time+0.0001)))*(Mp/M_PI) * 
			BesselK(0,fabs(kval*rad)) * kval*(sqrt(SQR(x1-Xplanet))/rad)*sin(kval*(x2-Yplanet));	

//	phi = 1.0*MIN(1.0,(ramp_time/(insert_time+0.0001)))*Mp *
//			sin(Rsoft+sqrt(SQR(x2-Yplanet))) * exp(-sqrt(SQRT(x1-Xplanet)));
  return phi;
}
/*------------------------------------------------------------------------------
 *! \fn static Real UnstratifiedDisk(const Real x1, const Real x2,const Real x3)
 *  \brief tidal potential in shearing box
 */

static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifdef SHEARING_BOX
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
#endif
  return phi;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_dV2(const GridS *pG, const int i, const int j, 
 *			     const int k)
 *  \brief Computes delta(Vy) 
 */

static Real expr_dV2(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef SHEARING_BOX
#ifdef FARGO
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
#endif
}
/*! \fn static Real expr_Lflux(const GridS *pG,const int i,const int j, 
 *				   const int k)
 *  \brief Output angluar momentum flux at given point, used for tab & vtk files.*/
static Real expr_Lflux(const GridS *pG,const int i,const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef SHEARING_BOX
#ifdef FARGO
  return pG->U[k][j][i].M1*(pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  return pG->U[k][j][i].M1*
    (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
#endif
}

/*---------------------------------------------------------------------------*/
/* Hydro history variables:
 * hst_rho_Vx_dVy: Reynolds stress, added as history variable.
 * hst_rho_dVy2: KE in y-velocity fluctuations
 * hst_E_total: total energy (including tidal potential).
 */

/*! \fn static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j, 
 *				   const int k)
 *  \brief Reynolds stress, added as history variable.*/
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef SHEARING_BOX
#ifdef FARGO
  return pG->U[k][j][i].M1*(pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  return pG->U[k][j][i].M1*
    (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
#endif
}

/*! \fn static Real hst_rho_dVy2(const GridS *pG, const int i, const int j,  
 *				const int k)
 *  \brief KE in y-velocity fluctuations */
static Real hst_rho_dVy2(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,dVy;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef SHEARING_BOX
#ifdef FARGO
  dVy = (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  dVy = (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
#endif
  return pG->U[k][j][i].d*dVy*dVy;
}

#ifdef ADIABATIC
/*! \fn static Real hst_E_total(const GridS *pG, const int i, const int j, 
 *			        const int k)
 *  \brief total energy (including tidal potential). */
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,phi;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  phi = UnstratifiedDisk(x1, x2, x3);

  return pG->U[k][j][i].E + pG->U[k][j][i].d*phi;
}
#endif /* ADIABATIC */

/*---------------------------------------------------------------------------*/
/*! \fn void constant_iib(GridS *pGrid)
 *  \brief Sets boundary condition on left X boundary (iib) to constant
 * state (initial values).
 */

void constant_iib(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  Real x1,x2,x3;
  Real den = 1.0, pres = 1.0e-6;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
      cc_pos(pGrid,(is-i),j,k,&x1,&x2,&x3);

/* Initialize d, M, and P.  With FARGO do not initialize the background shear */

      pGrid->U[k][j][is-i].d  = den;
      pGrid->U[k][j][is-i].M1 = 0.0;
      pGrid->U[k][j][is-i].M2 = 0.0;
#ifdef SHEARING_BOX
#ifndef FARGO
      pGrid->U[k][j][is-i].M2 -= den*(qshear*Omega_0*x1);
#endif
#endif
      pGrid->U[k][j][is-i].M3 = 0.0;
#ifdef ADIABATIC
      pGrid->U[k][j][is-i].E = pres/Gamma_1
        + 0.5*(SQR(pGrid->U[k][j][is-i].M1) + SQR(pGrid->U[k][j][is-i].M2)
             + SQR(pGrid->U[k][j][is-i].M3))/den;
#endif
      }
    }
  } 
}

/*---------------------------------------------------------------------------*/
/*! \fn void constant_oib(GridS *pGrid)
 *  \brief  Sets boundary condition on right X boundary (oib) to constant
 * state (initial values).
 */

void constant_oib(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  Real x1,x2,x3;
  Real den = 1.0, pres = 1.0e-6;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
      cc_pos(pGrid,(ie+i),j,k,&x1,&x2,&x3);

/* Initialize d, M, and P.  With FARGO do not initialize the background shear */

      pGrid->U[k][j][ie+i].d  = den;
      pGrid->U[k][j][ie+i].M1 = 0.0;
      pGrid->U[k][j][ie+i].M2 = 0.0;
#ifdef SHEARING_BOX
#ifndef FARGO
      pGrid->U[k][j][ie+i].M2 -= den*(qshear*Omega_0*x1);
#endif
#endif
      pGrid->U[k][j][ie+i].M3 = 0.0;
#ifdef ADIABATIC
      pGrid->U[k][j][ie+i].E = pres/Gamma_1
        + 0.5*(SQR(pGrid->U[k][j][ie+i].M1) + SQR(pGrid->U[k][j][ie+i].M2)
             + SQR(pGrid->U[k][j][ie+i].M3))/den;
#endif
      }
    }
  }
}
/*---------------------------------------------------------------------------*/
/*! \fn static void flux_output_func(MeshS *pM, OutputS *pOut) 
 *  \brief  New output format which outputs y-integrated angular momentum fluxes
 *  Currently can only be used with 1 proc and 1 domain.
 */
static void flux_output_func(MeshS *pM, OutputS *pOut)
{
  GridS *pG=pM->Domain[0][0].Grid;
  int nx1,nx2,nx3, ind, i,j,k, ind_i,ind_k;
  Real lx1, lx2, lx3;
  PrimS W[7],Ws[7],We[7];
  Real x1[7],x2[7],x3[7],xs1[7],xs2[7],xs3[7],xe1[7],xe2[7],xe3[7];
  Real dmin, dmax;
  Real **Fluxx=NULL;
  Real **FluxH=NULL;
  Real **FluxNu=NULL;
  Real **Th=NULL;
  Real **outCoordsx1=NULL;
  Real **outCoordsx3=NULL;
  Real **vx=NULL;
  Real **vy=NULL;
  Real **sigvx=NULL;
  Real **osigx=NULL;
  Real **davg=NULL;
  
  
  FILE *pfile;
  char *fname;
	
  nx1 = pG->Nx[0]; nx3 = pG->Nx[2]; nx2 = pG->Nx[1];
  lx1 = pG->MaxX[0] - pG->MinX[0];
  lx2 = pG->MaxX[1] - pG->MinX[1];
  lx3 = pG->MaxX[2] - pG->MinX[2];
 printf("%d, %d, %d, %g, %g, %g\n",nx1,nx2,nx3,lx1,lx2,lx3);
#ifdef MPI_PARALLEL  
  printf("%d, %d, %d, %g, %g, %g \n", myID_Comm_world,nx1,nx3,lx1,lx2,lx3); 

  printf("IND %d: (%d,%d), (%d,%d)\n",myID_Comm_world,pG->is,pG->js,pG->ie,pG->je);

#endif

  Fluxx=(Real **)calloc_2d_array(nx3,nx1,sizeof(Real));
  if (Fluxx == NULL) return;
  FluxH=(Real **)calloc_2d_array(nx3,nx1,sizeof(Real));
  if (FluxH == NULL) return;
  FluxNu=(Real **)calloc_2d_array(nx3,nx1,sizeof(Real));
  if (FluxNu == NULL) return;
  Th=(Real **)calloc_2d_array(nx3,nx1,sizeof(Real));
  if (Th == NULL) return;
  outCoordsx1=(Real **)calloc_2d_array(nx3,nx1,sizeof(Real));
  if (outCoordsx1 == NULL) return;
  outCoordsx3=(Real **)calloc_2d_array(nx3,nx1,sizeof(Real));
  if (outCoordsx3 == NULL) return;
  vx=(Real **)calloc_2d_array(nx3,nx1,sizeof(Real));
  if (vx == NULL) return;
  vy=(Real **)calloc_2d_array(nx3,nx1,sizeof(Real));
  if (vy == NULL) return;
  sigvx=(Real **)calloc_2d_array(nx3,nx1,sizeof(Real));
  if (sigvx == NULL) return;
  osigx=(Real **)calloc_2d_array(nx3,nx1,sizeof(Real));
  if (osigx == NULL) return;
  davg=(Real **)calloc_2d_array(nx3,nx1,sizeof(Real));
  if (davg == NULL) return;
/* Open file and write header */  
  
  if((fname = ath_fname(NULL,pM->outfilename,NULL,NULL,num_digit,
            pOut->num,pOut->id,"tab")) == NULL){
          ath_error("[dump_tab]: Error constructing filename\n");
  }
  if((pfile = fopen(fname,"w")) == NULL){
  	ath_error("[dump_tab]: Unable to open ppm file %s\n",fname);
  }
  free(fname);
  
#ifdef MPI_PARALLEL
  if(myID_Comm_world == 0) 
#endif
  	fprintf(pfile,"#t=%12.8e	x,FH,Fx,Fnu,Th,vx,vy,davg,sigvx,omsigx \n", pOut->t);

/* Compute y-integrated fluxes explicitly.
 * For the derivatives use a central difference method.
 * For the integration use a composite trapezoid method.
 * Both of these make use of the ghost cells for the boundary values.
 * 1	FluxH = < d * vx1 *vx2 >
 * 2	Fluxx = < 0.5 * omega * d * x1 * vx1 >
 * 3	FluxNu = - < nu_iso * d/dx vx2 >
 * 4	Th = - < d * d/dy phi >
 * 5	vx = < vx1 >
 * 6	vy = < vx2 >
 * 7	davg = < d >
 * 8	sigvx = < d * vx1 >
 * 9	osigx = < 0.5 * sig * omega * x1 >
 *
 *         
 *		x	4   x
 *      1   0	2
 *		x	3	x
 *
 * The variables are stored up as arrays of primitives W[7].
 * W = ( W[k][j][i], W[k][j][i-1], W[k][j][i+1], W[k][j-1][i], 
 *		 W[k][j+1][i], W[k-1][j][i], W[k+1][j][i] )
 * This is needed for the integration and differentiation.
 *
 * The background shear is taken out of vy when doing computations.
*/

  for(k=pG->ks; k <=pG->ke; k++) {
  	for(i=pG->is; i<= pG->ie; i++) {
		ind_i=i-pG->is; ind_k=k-pG->ks; 
		for(ind=0; ind < 7; ind++) {
			if (ind==0) {
				cc_pos(pG,i,pG->js,k,&xs1[ind],&xs2[ind],&xs3[ind]);
				cc_pos(pG,i,pG->je,k,&xe1[ind],&xe2[ind],&xe3[ind]);
				Ws[ind] = Cons_to_Prim(&(pG->U[k][pG->js][i])); 
				We[ind] = Cons_to_Prim(&(pG->U[k][pG->je][i])); 
			}
			if (ind > 0 && ind < 3) {
				cc_pos(pG,i+2*ind-3,pG->js,k,&xs1[ind],&xs2[ind],&xs3[ind]);
				cc_pos(pG,i+2*ind-3,pG->je,k,&xe1[ind],&xe2[ind],&xe3[ind]);
				Ws[ind] = Cons_to_Prim(&(pG->U[k][pG->js][i+2*ind-3])); 
				We[ind] = Cons_to_Prim(&(pG->U[k][pG->je][i+2*ind-3]));
			}
			if (ind > 2 && ind < 5) {
				cc_pos(pG,i,pG->js+2*ind-7,k,&xs1[ind],&xs2[ind],&xs3[ind]);
				cc_pos(pG,i,pG->je+2*ind-7,k,&xe1[ind],&xe2[ind],&xe3[ind]);
				Ws[ind] = Cons_to_Prim(&(pG->U[k][pG->js+2*ind-7][i])); 
				We[ind] = Cons_to_Prim(&(pG->U[k][pG->je+2*ind-7][i]));
			}
			if (ind > 4 && ind < 7) {
				if (pG->MinX[2] == pG->MaxX[2]) {
					cc_pos(pG,i,pG->js,k,&xs1[ind],&xs2[ind],&xs3[ind]);
					cc_pos(pG,i,pG->je,k,&xe1[ind],&xe2[ind],&xe3[ind]);
					Ws[ind] = Cons_to_Prim(&(pG->U[k][pG->js][i])); 
					We[ind] = Cons_to_Prim(&(pG->U[k][pG->je][i]));
				}
				else {
					cc_pos(pG,i,pG->js,k+2*ind-11,&xs1[ind],&xs2[ind],&xs3[ind]);
					cc_pos(pG,i,pG->je,k+2*ind-11,&xe1[ind],&xe2[ind],&xe3[ind]);
					Ws[ind] = Cons_to_Prim(&(pG->U[k+2*ind-11][pG->js][i])); 
					We[ind] = Cons_to_Prim(&(pG->U[k+2*ind-11][pG->je][i]));
				}
			}
			Ws[ind].V2 = Ws[ind].V2 + qshear*Omega_0*xs1[ind];
			We[ind].V2 = We[ind].V2 + qshear*Omega_0*xe1[ind];
/* If using d' instead of d then put in
 * W[ind] = W[ind]->d - d0;
 */
		}
/* Set initial values for the integration */
		FluxH[ind_k][ind_i] = 0.5*( Ws[0].d * Ws[0].V1 * Ws[0].V2 + 
						    We[0].d * We[0].V1 * We[0].V2 );  	
		
		Fluxx[ind_k][ind_i] = 0.5*( Ws[0].d * Ws[0].V1 * xs1[0] +
							We[0].d * We[0].V1 * xe1[0] );
		
#ifdef VISCOSITY 
		FluxNu[ind_k][ind_i] = 0.5*( Ws[2].V2 - Ws[1].V2 +
							 We[2].V2 - We[1].V2 ); 
		
#else
		FluxNu[ind_k][ind_i] = 0;
#endif   	
  	
  		Th[ind_k][ind_i] = 0.5*( Ws[0].d * dx2PlanetPot(xs1[0],xs2[0],xs3[0]) +
  						 We[0].d * dx2PlanetPot(xe1[0],xe2[0],xe3[0]) );
		
		vx[ind_k][ind_i] = 0.5*( Ws[0].V1 + We[0].V1 );
		vy[ind_k][ind_i] = 0.5*( Ws[0].V2 + We[0].V2 );
		sigvx[ind_k][ind_i] = 0.5*( Ws[0].d * Ws[0].V1 + We[0].d * We[0].V1 );
		osigx[ind_k][ind_i] = 0.5*( Ws[0].d * xs1[0] + We[0].d * xe1[0] );
		davg[ind_k][ind_i] = 0.5*(Ws[0].d + We[0].d);
  	 	for(j=(pG->js)+1; j < pG->je; j++) {
  			for(ind=0; ind < 7; ind++) {
				if (ind==0) {
					cc_pos(pG,i,j,k,&x1[ind],&x2[ind],&x3[ind]);
					W[ind] = Cons_to_Prim(&(pG->U[k][j][i])); 
				}
				if (ind > 0 && ind < 3) {
					cc_pos(pG,i+2*ind-3,j,k,&x1[ind],&x2[ind],&x3[ind]);
					W[ind] = Cons_to_Prim(&(pG->U[k][j][i+2*ind-3])); 
				}
				if (ind > 2 && ind < 5) {
					cc_pos(pG,i,j+2*ind-7,k,&x1[ind],&x2[ind],&x3[ind]);
					W[ind] = Cons_to_Prim(&(pG->U[k][j+2*ind-7][i])); 
				}
				if (ind > 4 && ind < 7) {
					if (pG->MinX[2] == pG->MaxX[2]) {
						cc_pos(pG,i,j,k,&x1[ind],&x2[ind],&x3[ind]);
						W[ind] = Cons_to_Prim(&(pG->U[k][j][i])); 
					}
					else {
						cc_pos(pG,i,j,k+2*ind-11,&x1[ind],&x2[ind],&x3[ind]);
						W[ind] = Cons_to_Prim(&(pG->U[k+2*ind-11][j][i])); 
					}
				}
				W[ind].V2 = W[ind].V2 + qshear*Omega_0*x1[ind];
				
			}	 
			
			FluxH[ind_k][ind_i] += W[0].d * W[0].V1 * W[0].V2;
			Fluxx[ind_k][ind_i] += W[0].d * W[0].V1 * x1[0];
#ifdef VISCOSITY
			FluxNu[ind_k][ind_i] += W[0].d * (W[2].V2 - W[1].V1);
#endif
			Th[ind_k][ind_i] += W[0].d * dx2PlanetPot(x1[0],x2[0],x3[0]);
			vx[ind_k][ind_i] += W[0].V1;
			vy[ind_k][ind_i] += W[0].V2;
			sigvx[ind_k][ind_i] += W[0].d * W[0].V1;
			osigx[ind_k][ind_i] += W[0].d * x1[0];
			davg[ind_k][ind_i] += W[0].d;
			
  		}	
  		FluxH[ind_k][ind_i] *= (pG->dx2)/lx2;
  		Fluxx[ind_k][ind_i] *= .5*Omega_0*(pG->dx2)/lx2;
#ifdef	VISCOSITY
		FluxNu[ind_k][ind_i] *= -nu_iso*(pG->dx2)/(2*lx2*(pG->dx1));
#endif 
  		Th[ind_k][ind_i] *= -1.0*(pG->dx2)/lx2;				
  		vx[ind_k][ind_i] *= (pG->dx2)/lx2;
  		vy[ind_k][ind_i] *= (pG->dx2)/lx2;
  		sigvx[ind_k][ind_i] *= (pG->dx2)/lx2;
  		osigx[ind_k][ind_i] *= (0.5*Omega_0*(pG->dx2))/lx2;
  		davg[ind_k][ind_i] *= (pG->dx2)/lx2;
  	 	outCoordsx1[ind_k][ind_i]=x1[0]; outCoordsx3[ind_k][ind_i]=x3[0];
  	}
  }  

/* Quantities are ready to be written to output file
 * Format (Not outputting x3 at the moment: 
 * x1	(x3)	FluxH	Fluxx	Fluxnu	Th	vx	vy	davg	sigvx	osigx
*/
  
  
  for(k=pG->ks; k<=pG->ke; k++) {
  	for(i=pG->is; i<=pG->ie; i++) {
 		ind_k=k-pG->ks; ind_i=i-pG->is;
  		if (lx3==0) {
  			fprintf(pfile,"%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
  					outCoordsx1[ind_k][ind_i],FluxH[ind_k][ind_i],
  					Fluxx[ind_k][ind_i],FluxNu[ind_k][ind_i],Th[ind_k][ind_i],
  					vx[ind_k][ind_i],vy[ind_k][ind_i],davg[ind_k][ind_i],sigvx[ind_k][ind_i],
  					osigx[ind_k][ind_i]);
  		}
  		else {
  			fprintf(pfile,"%12.8e	%12.8e	%12.8e	%12.8e	%12.8e	%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
  					outCoordsx1[ind_k][ind_i],outCoordsx3[ind_k][ind_i],FluxH[ind_k][ind_i],
  					Fluxx[ind_k][ind_i],FluxNu[ind_k][ind_i],Th[ind_k][ind_i],
  					vx[ind_k][ind_i],vy[ind_k][ind_i],davg[ind_k][ind_i],sigvx[ind_k][ind_i],
  					osigx[ind_k][ind_i]);
  		}
  	}
  }
  
  fclose(pfile);
  free_2d_array(Fluxx); 
  free_2d_array(FluxH); 
  free_2d_array(FluxNu); 
  free_2d_array(Th); 
  free_2d_array(outCoordsx1);
  free_2d_array(outCoordsx3);
  free_2d_array(vx);
  free_2d_array(vy);
  free_2d_array(sigvx);
  free_2d_array(osigx);
  free_2d_array(davg);
  return;
}

static Real BesselK(const int order,const Real x)
{
	Real y,ans; 
	if (order==0) {
		if (x <= 2.0) { 
		y=x*x/4.0;
		ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420 +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2 +y*(0.10750e-3+y*0.74e-5))))));
		} 
		else { 
			y=2.0/x;
			ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1 +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2 +y*(-0.251540e-2+y*0.53208e-3))))));
		}
	}
// 	if (order == 1) {
// 		if (x <= 2.0) { 
// 			y=x*x/4.0;
// 			ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144 +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1 +y*(-0.110404e-2+y*(-0.4686e-4)))))));
// 		} 
// 		else { 
// 			y=2.0/x;
// 		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619 +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2 +y*(0.325614e-2+y*(-0.68245e-3)))))));
// 		}
// 	}
//	if (order !=1 || order !-0) ans=0;
	
	return ans; 
}

static Real bessi0(const Real x)
{
	Real ax,ans; Real y;
	if ((ax=fabs(x)) < 3.75) { 
		y=x/3.75;
		y*=y; 
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
		+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))); 
	} 
	else {
		y=3.75/ax; 
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2 +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1 +y*0.392377e-2))))))));
	}
	return ans;
}
