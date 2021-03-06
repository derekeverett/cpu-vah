/*
 * SourceTerms.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */
#include <stdlib.h> // for system.exit
#include <stdio.h> // for printf
#include <math.h> // for math functions

#include "edu/osu/rhic/trunk/hydro/SourceTerms.h"
#include "edu/osu/rhic/core/util/FiniteDifference.h"
#include "edu/osu/rhic/trunk/hydro/EnergyMomentumTensor.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

#include "edu/osu/rhic/trunk/hydro/FullyDiscreteKurganovTadmorScheme.h" // for const params

#include "edu/osu/rhic/trunk/eos/EquationOfState.h" // for bulk terms
#include "edu/osu/rhic/trunk/hydro/AnisotropicDistributionFunctions.h"

#include "edu/osu/rhic/core/muscl/FluxLimiter.h"

#define USE_APPROX_DERIVATIVE //use approximate derivative with minmod flux limiter for derivatives of primary variables and pi^munu

// paramters for the analytic parameterization of the bulk viscosity \zeta/S
#define A_1 -13.77
#define A_2 27.55
#define A_3 13.45

#define LAMBDA_1 0.9
#define LAMBDA_2 0.25
#define LAMBDA_3 0.9
#define LAMBDA_4 0.22

#define SIGMA_1 0.025
#define SIGMA_2 0.13
#define SIGMA_3 0.0025
#define SIGMA_4 0.022

inline PRECISION bulkViscosityToEntropyDensity(PRECISION T) {
	PRECISION x = T/1.01355;
	if(x > 1.05)
		return LAMBDA_1*exp(-(x-1)/SIGMA_1) + LAMBDA_2*exp(-(x-1)/SIGMA_2)+0.001;
	else if(x < 0.995)
		return LAMBDA_3*exp((x-1)/SIGMA_3)+ LAMBDA_4*exp((x-1)/SIGMA_4)+0.03;
	else
		return A_1*x*x + A_2*x - A_3;
}

void loadSourceTermsX(const PRECISION * const __restrict__ I, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s,
PRECISION d_dx
) {
	//=========================================================
	// spatial derivatives of the conserved variables \pi^{\mu\nu}
	//=========================================================

	#ifndef USE_APPROX_DERIVATIVE

	PRECISION facX = 1 / d_dx / 2;
	int ptr = 25; // 5 * n (with n = 4 corresponding to pitt)
	PRECISION dxpitt = (*(I + ptr + 3) - *(I + ptr + 1)) * facX;
	ptr += 5;
	PRECISION dxpitx = (*(I + ptr + 3) - *(I + ptr + 1)) * facX;
	ptr += 5;
	PRECISION dxpity = (*(I + ptr + 3) - *(I + ptr + 1)) * facX;
	ptr += 5;
	PRECISION dxpitn = (*(I + ptr + 3) - *(I + ptr + 1)) * facX;
	ptr += 5;
	PRECISION dxpixx = (*(I + ptr + 3) - *(I + ptr + 1)) * facX;
	ptr += 5;
	PRECISION dxpixy = (*(I + ptr + 3) - *(I + ptr + 1)) * facX;
	ptr += 5;
	PRECISION dxpixn = (*(I + ptr + 3) - *(I + ptr + 1)) * facX;
	ptr += 20;

	#else

	int ptr = 25; // 5 * n (with n = 4 corresponding to pitt)
	PRECISION dxpitt = approximateDerivative( *(I + ptr + 1), *(I + ptr + 2), *(I + ptr + 3) )/d_dx;
	ptr += 5;
	PRECISION dxpitx = approximateDerivative( *(I + ptr + 1), *(I + ptr + 2), *(I + ptr + 3) )/d_dx;
	ptr += 5;
	PRECISION dxpity = approximateDerivative( *(I + ptr + 1), *(I + ptr + 2), *(I + ptr + 3) )/d_dx;
	ptr += 5;
	PRECISION dxpitn = approximateDerivative( *(I + ptr + 1), *(I + ptr + 2), *(I + ptr + 3) )/d_dx;
	ptr += 5;
	PRECISION dxpixx = approximateDerivative( *(I + ptr + 1), *(I + ptr + 2), *(I + ptr + 3) )/d_dx;
	ptr += 5;
	PRECISION dxpixy = approximateDerivative( *(I + ptr + 1), *(I + ptr + 2), *(I + ptr + 3) )/d_dx;
	ptr += 5;
	PRECISION dxpixn = approximateDerivative( *(I + ptr + 1), *(I + ptr + 2), *(I + ptr + 3) )/d_dx;
	ptr += 20;

	#endif

	PRECISION ut = u->ut[s];
	PRECISION ux = u->ux[s];

	//=========================================================
	// set dx terms in the source terms
	//=========================================================
	PRECISION vx = ux / ut;
#ifndef PI
	S[0] = dxpitt*vx - dxpitx;
	S[1] = dxpitx*vx - dxpixx;
#else
	ptr=(NUMBER_CONSERVED_VARIABLES-1)*5;
	#ifndef USE_APPROX_DERIVATIVE
	PRECISION dxPi = (*(I + ptr + 3) - *(I + ptr + 1)) * facX;
	#else
	PRECISION dxPi = approximateDerivative( *(I + ptr + 1), *(I + ptr + 2), *(I + ptr + 3) )/d_dx;
	#endif
	dxPi *= 1.5;
	S[0] = dxpitt * vx - dxpitx - vx * dxPi;
	S[1] = dxpitx * vx - dxpixx - dxPi;
#endif
	S[2] = dxpity * vx - dxpixy;
	S[3] = dxpitn * vx - dxpixn;
}

void loadSourceTermsY(const PRECISION * const __restrict__ J, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s,
PRECISION d_dy
) {
	//=========================================================
	// spatial derivatives of the conserved variables \pi^{\mu\nu}
	//=========================================================

	#ifndef USE_APPROX_DERIVATIVE

	PRECISION facY = 1 / d_dy / 2;
	int ptr = 25; // 5 * n (with n = 4 corresponding to pitt)
	PRECISION dypitt = (*(J + ptr + 3) - *(J + ptr + 1)) * facY;
	ptr += 5;
	PRECISION dypitx = (*(J + ptr + 3) - *(J + ptr + 1)) * facY;
	ptr += 5;
	PRECISION dypity = (*(J + ptr + 3) - *(J + ptr + 1)) * facY;
	ptr += 5;
	PRECISION dypitn = (*(J + ptr + 3) - *(J + ptr + 1)) * facY;
	ptr += 10;
	PRECISION dypixy = (*(J + ptr + 3) - *(J + ptr + 1)) * facY;
	ptr += 10;
	PRECISION dypiyy = (*(J + ptr + 3) - *(J + ptr + 1)) * facY;
	ptr += 5;
	PRECISION dypiyn = (*(J + ptr + 3) - *(J + ptr + 1)) * facY;
	ptr += 10;

	#else

	int ptr = 25; // 5 * n (with n = 4 corresponding to pitt)
	PRECISION dypitt = approximateDerivative( *(J + ptr + 1), *(J + ptr + 2), *(J + ptr + 3) )/d_dy;
	ptr += 5;
	PRECISION dypitx = approximateDerivative( *(J + ptr + 1), *(J + ptr + 2), *(J + ptr + 3) )/d_dy;
	ptr += 5;
	PRECISION dypity = approximateDerivative( *(J + ptr + 1), *(J + ptr + 2), *(J + ptr + 3) )/d_dy;
	ptr += 5;
	PRECISION dypitn = approximateDerivative( *(J + ptr + 1), *(J + ptr + 2), *(J + ptr + 3) )/d_dy;
	ptr += 10;
	PRECISION dypixy = approximateDerivative( *(J + ptr + 1), *(J + ptr + 2), *(J + ptr + 3) )/d_dy;
	ptr += 10;
	PRECISION dypiyy = approximateDerivative( *(J + ptr + 1), *(J + ptr + 2), *(J + ptr + 3) )/d_dy;
	ptr += 5;
	PRECISION dypiyn = approximateDerivative( *(J + ptr + 1), *(J + ptr + 2), *(J + ptr + 3) )/d_dy;
	ptr += 10;

	#endif

	PRECISION ut = u->ut[s];
	PRECISION uy = u->uy[s];

	//=========================================================
	// set dy terms in the source terms
	//=========================================================
	PRECISION vy = uy / ut;
#ifndef PI
	S[0] = dypitt*vy - dypity;
	S[2] = dypity*vy - dypiyy;
#else
	ptr=(NUMBER_CONSERVED_VARIABLES-1)*5;
	#ifndef USE_APPROX_DERIVATIVE
	PRECISION dyPi = (*(J + ptr + 3) - *(J + ptr + 1)) * facY;
	#else
	PRECISION dyPi = approximateDerivative( *(J + ptr + 1), *(J + ptr + 2), *(J + ptr + 3) )/d_dy;
	#endif
	dyPi *= 1.5;
	S[0] = dypitt * vy - dypity - vy * dyPi;
	S[2] = dypity * vy - dypiyy - dyPi;
#endif
	S[1] = dypitx * vy - dypixy;
	S[3] = dypitn * vy - dypiyn;
}

void loadSourceTermsZ(const PRECISION * const __restrict__ K, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s, PRECISION t,
PRECISION d_dz
) {
	//=========================================================
	// spatial derivatives of the conserved variables \pi^{\mu\nu}
	//=========================================================

	#ifndef USE_APPROX_DERIVATIVE

	PRECISION facZ = 1 / d_dz / 2;
	int ptr = 25; // 5 * n (with n = 4 corresponding to pitt)
	PRECISION dnpitt = (*(K + ptr + 3) - *(K + ptr + 1)) * facZ;
	ptr += 5;
	PRECISION dnpitx = (*(K + ptr + 3) - *(K + ptr + 1)) * facZ;
	ptr += 5;
	PRECISION dnpity = (*(K + ptr + 3) - *(K + ptr + 1)) * facZ;
	ptr += 5;
	PRECISION dnpitn = (*(K + ptr + 3) - *(K + ptr + 1)) * facZ;
	ptr += 15;
	PRECISION dnpixn = (*(K + ptr + 3) - *(K + ptr + 1)) * facZ;
	ptr += 10;
	PRECISION dnpiyn = (*(K + ptr + 3) - *(K + ptr + 1)) * facZ;
	ptr += 5;
	PRECISION dnpinn = (*(K + ptr + 3) - *(K + ptr + 1)) * facZ;
	ptr += 5;

	#else

	int ptr = 25; // 5 * n (with n = 4 corresponding to pitt)
	PRECISION dnpitt = approximateDerivative( *(K + ptr + 1), *(K + ptr + 2), *(K + ptr + 3) )/d_dz;
	ptr += 5;
	PRECISION dnpitx = approximateDerivative( *(K + ptr + 1), *(K + ptr + 2), *(K + ptr + 3) )/d_dz;
	ptr += 5;
	PRECISION dnpity = approximateDerivative( *(K + ptr + 1), *(K + ptr + 2), *(K + ptr + 3) )/d_dz;
	ptr += 5;
	PRECISION dnpitn = approximateDerivative( *(K + ptr + 1), *(K + ptr + 2), *(K + ptr + 3) )/d_dz;
	ptr += 15;
	PRECISION dnpixn = approximateDerivative( *(K + ptr + 1), *(K + ptr + 2), *(K + ptr + 3) )/d_dz;
	ptr += 10;
	PRECISION dnpiyn = approximateDerivative( *(K + ptr + 1), *(K + ptr + 2), *(K + ptr + 3) )/d_dz;
	ptr += 5;
	PRECISION dnpinn = approximateDerivative( *(K + ptr + 1), *(K + ptr + 2), *(K + ptr + 3) )/d_dz;
	ptr += 5;

	#endif

	PRECISION ut = u->ut[s];
	PRECISION un = u->un[s];

	//=========================================================
	// set dn terms in the source terms
	//=========================================================
	PRECISION vn = un / ut;
#ifndef PI
	S[0] = dnpitt*vn - dnpitn;
	S[3] = dnpitn*vn - dnpinn;
#else
	ptr=(NUMBER_CONSERVED_VARIABLES-1)*5;
	#ifndef USE_APPROX_DERIVATIVE
	PRECISION dnPi = (*(K + ptr + 3) - *(K + ptr + 1)) * facZ;
	#else
	PRECISION dnPi = approximateDerivative( *(K + ptr + 1), *(K + ptr + 2), *(K + ptr + 3) )/d_dz;
	#endif
	dnPi *= 1.5;
	S[0] = dnpitt * vn - dnpitn - vn * dnPi;
	S[3] = dnpitn * vn - dnpinn - dnPi / powf(t, 2);
#endif
	S[1] = dnpitx * vn - dnpixn;
	S[2] = dnpity * vn - dnpiyn;
}

void loadSourceTerms(const PRECISION * const __restrict__ Q, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u,
PRECISION utp, PRECISION uxp, PRECISION uyp, PRECISION unp,
PRECISION t, const PRECISION * const __restrict__ evec, const PRECISION * const __restrict__ pvec,
int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_etabar, PRECISION d_dt, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz,
int i, int j, int k, double x, double y, double z,
const CONSERVED_VARIABLES * const __restrict__ currentVars
) {
	//=========================================================
	// conserved variables
	//=========================================================
	PRECISION ttt = Q[0];
	PRECISION ttx = Q[1];
	PRECISION tty = Q[2];
	PRECISION ttn = Q[3];
	PRECISION pl = Q[4];
	// \pi^{\mu\nu}_{\perp}
#ifdef PIMUNU
	PRECISION pitt = Q[5];
	PRECISION pitx = Q[6];
	PRECISION pity = Q[7];
	PRECISION pitn = Q[8];
	PRECISION pixx = Q[9];
	PRECISION pixy = Q[10];
	PRECISION pixn = Q[11];
	PRECISION piyy = Q[12];
	PRECISION piyn = Q[13];
	PRECISION pinn = Q[14];
#else
	PRECISION pitt = 0;
	PRECISION pitx = 0;
	PRECISION pity = 0;
	PRECISION pitn = 0;
	PRECISION pixx = 0;
	PRECISION pixy = 0;
	PRECISION pixn = 0;
	PRECISION piyy = 0;
	PRECISION piyn = 0;
	PRECISION pinn = 0;
#endif
	// W^{\mu}_{\perp z}
#ifdef W_TZ_MU
	PRECISION WtTz = Q[15];
	PRECISION WxTz = Q[16];
	PRECISION WyTz = Q[17];
	PRECISION WnTz = Q[18];
#else
	PRECISION WtTz = 0;
	PRECISION WxTz = 0;
	PRECISION WyTz = 0;
	PRECISION WnTz = 0;
#endif
	// \Pi
#ifdef PI
	PRECISION Pi = Q[NUMBER_CONSERVED_VARIABLES-1];
#else
	PRECISION Pi = 0;
#endif

	//=========================================================
	// primary variables
	//=========================================================
	PRECISION *utvec = u->ut;
	PRECISION *uxvec = u->ux;
	PRECISION *uyvec = u->uy;
	PRECISION *unvec = u->un;

	PRECISION e = evec[s];
	PRECISION p = pvec[s];
	PRECISION ut = utvec[s];
	PRECISION ux = uxvec[s];
	PRECISION uy = uyvec[s];
	PRECISION un = unvec[s];

	int stride = d_ncx * d_ncy;

	//=========================================================
	// spatial derivatives of primary variables
	//=========================================================
	PRECISION facX = 1/d_dx/2;
	PRECISION facY = 1/d_dy/2;
	PRECISION facZ = 1/d_dz/2;

	#ifndef USE_APPROX_DERIVATIVE
	//use central finite difference

	// dx of u^{\mu} components
	PRECISION dxut = (*(utvec + s + 1) - *(utvec + s - 1)) * facX;
	PRECISION dxux = (*(uxvec + s + 1) - *(uxvec + s - 1)) * facX;
	PRECISION dxuy = (*(uyvec + s + 1) - *(uyvec + s - 1)) * facX;
	PRECISION dxun = (*(unvec + s + 1) - *(unvec + s - 1)) * facX;

	// dy of u^{\mu} components
	PRECISION dyut = (*(utvec + s + d_ncx) - *(utvec + s - d_ncx)) * facY;
	PRECISION dyux = (*(uxvec + s + d_ncx) - *(uxvec + s - d_ncx)) * facY;
	PRECISION dyuy = (*(uyvec + s + d_ncx) - *(uyvec + s - d_ncx)) * facY;
	PRECISION dyun = (*(unvec + s + d_ncx) - *(unvec + s - d_ncx)) * facY;

	// dn of u^{\mu} components
	int stride = d_ncx * d_ncy;
	PRECISION dnut = (*(utvec + s + stride) - *(utvec + s - stride)) * facZ;
	PRECISION dnux = (*(uxvec + s + stride) - *(uxvec + s - stride)) * facZ;
	PRECISION dnuy = (*(uyvec + s + stride) - *(uyvec + s - stride)) * facZ;
	PRECISION dnun = (*(unvec + s + stride) - *(unvec + s - stride)) * facZ;

	// pressure
	PRECISION dxp = (*(pvec + s + 1) - *(pvec + s - 1)) * facX;
	PRECISION dyp = (*(pvec + s + d_ncx) - *(pvec + s - d_ncx)) * facY;
	PRECISION dnp = (*(pvec + s + stride) - *(pvec + s - stride)) * facZ;

	// energy
	PRECISION dxe = (*(evec + s + 1) - *(evec + s - 1)) * facX;
	PRECISION dye = (*(evec + s + d_ncx) - *(evec + s - d_ncx)) * facY;
	PRECISION dne = (*(evec + s + stride) - *(evec + s - stride)) * facZ;

	#else
	//use minmod flux limiter to dampen oscillations

	PRECISION e_right =  evec[s+1];
	PRECISION p_right =  pvec[s+1];
	PRECISION ut_right = utvec[s+1];
	PRECISION ux_right = uxvec[s+1];
	PRECISION uy_right = uyvec[s+1];
	PRECISION un_right = unvec[s+1];

	PRECISION e_left =  evec[s-1];
	PRECISION p_left =  pvec[s-1];
	PRECISION ut_left = utvec[s-1];
	PRECISION ux_left = uxvec[s-1];
	PRECISION uy_left = uyvec[s-1];
	PRECISION un_left = unvec[s-1];

	PRECISION e_top =  evec[s+d_ncx];
	PRECISION p_top =  pvec[s+d_ncx];
	PRECISION ut_top = utvec[s+d_ncx];
	PRECISION ux_top = uxvec[s+d_ncx];
	PRECISION uy_top = uyvec[s+d_ncx];
	PRECISION un_top = unvec[s+d_ncx];

	PRECISION e_bottom =  evec[s-d_ncx];
	PRECISION p_bottom =  pvec[s-d_ncx];
	PRECISION ut_bottom = utvec[s-d_ncx];
	PRECISION ux_bottom = uxvec[s-d_ncx];
	PRECISION uy_bottom = uyvec[s-d_ncx];
	PRECISION un_bottom = unvec[s-d_ncx];

	PRECISION e_forward =  evec[s+stride];
	PRECISION p_forward =  pvec[s+stride];
	PRECISION ut_forward = utvec[s+stride];
	PRECISION ux_forward = uxvec[s+stride];
	PRECISION uy_forward = uyvec[s+stride];
	PRECISION un_forward = unvec[s+stride];

	PRECISION e_backward =  evec[s-stride];
	PRECISION p_backward =  pvec[s-stride];
	PRECISION ut_backward = utvec[s-stride];
	PRECISION ux_backward = uxvec[s-stride];
	PRECISION uy_backward = uyvec[s-stride];
	PRECISION un_backward = unvec[s-stride];

	// dx of u^{\mu} components
	PRECISION dxut = approximateDerivative(ut_left, ut, ut_right) / d_dx;
	PRECISION dxux = approximateDerivative(ux_left, ux, ux_right) / d_dx;
	PRECISION dxuy = approximateDerivative(uy_left, uy, uy_right) / d_dx;
	PRECISION dxun = approximateDerivative(un_left, un, un_right) / d_dx;

	// dy of u^{\mu} components
	PRECISION dyut = approximateDerivative(ut_bottom, ut, ut_top) / d_dy;
	PRECISION dyux = approximateDerivative(ux_bottom, ux, ux_top) / d_dy;
	PRECISION dyuy = approximateDerivative(uy_bottom, uy, uy_top) / d_dy;
	PRECISION dyun = approximateDerivative(un_bottom, un, un_top) / d_dy;

	//dn of u^{\mu} components
	PRECISION dnut = approximateDerivative(ut_backward, ut, ut_forward) / d_dz;
	PRECISION dnux = approximateDerivative(ux_backward, ux, ux_forward) / d_dz;
	PRECISION dnuy = approximateDerivative(uy_backward, uy, uy_forward) / d_dz;
	PRECISION dnun = approximateDerivative(un_backward, un, un_forward) / d_dz;

	// pressure
	PRECISION dxp = approximateDerivative(p_left, p, p_right) / d_dx;
	PRECISION dyp = approximateDerivative(p_bottom, p, p_top) / d_dy;
	PRECISION dnp = approximateDerivative(p_backward, p, p_forward) / d_dz;

	// energy
	PRECISION dxe = approximateDerivative(e_left, e, e_right) / d_dx;
	PRECISION dye = approximateDerivative(e_bottom, e, e_top) / d_dy;
	PRECISION dne = approximateDerivative(e_backward, e, e_forward) / d_dz;

	#endif


	//=========================================================
	// Deriviatives of v
	//=========================================================
	PRECISION vx = ux/ut;
	PRECISION vy = uy/ut;
	PRECISION vn = un/ut;
	PRECISION dxvx = (dxux - vx * dxut)/ ut;
	PRECISION dyvy = (dyuy - vy * dyut)/ ut;
	PRECISION dnvn = (dnun - vn * dnut)/ ut;
	PRECISION dkvk = dxvx + dyvy + dnvn;

	PRECISION t2 = t*t;
	PRECISION t3 = t*t2;
	PRECISION un2 = un*un;
	PRECISION ut2 = ut*ut;

	/************************************************************************************\
	 * Gradient terms of the fluid velocity
	/************************************************************************************/
	// time derivatives of u
	PRECISION dtut = (ut - utp) / d_dt;
	PRECISION dtux = (ux - uxp) / d_dt;
	PRECISION dtuy = (uy - uyp) / d_dt;
	PRECISION dtun = (un - unp) / d_dt;

	// covariant derivatives
	PRECISION Dut = ut*dtut + ux*dxut + uy*dyut + un*dnut + t*un*un;
	PRECISION DuxUpper = ut*dtux + ux*dxux + uy*dyux + un*dnux;
	PRECISION Dux = -DuxUpper;
	PRECISION DuyUpper = ut*dtuy + ux*dxuy + uy*dyuy + un*dnuy;
	PRECISION Duy = -DuyUpper;
	PRECISION DunUpper = ut*dtun + ux*dxun + uy*dyun + un*dnun + 2*ut*un/t;
	PRECISION Dun = -t2*DunUpper;

	PRECISION dut = Dut -t*un*un;
	PRECISION dux = ut*dtux + ux*dxux + uy*dyux + un*dnux;
	PRECISION duy = ut*dtuy + ux*dxuy + uy*dyuy + un*dnuy;
	PRECISION dun = ut*dtun + ux*dxun + uy*dyun + un*dnun;

	// expansion rate
	PRECISION theta = ut / t + dtut + dxux + dyuy + dnun;

	// shear tensor
	PRECISION stt = -t * ut * un2 + (dtut - ut * dut) + (ut2 - 1) * theta / 3;
	PRECISION stx = -(t * un2 * ux) / 2 + (dtux - dxut) / 2 - (ux * dut + ut * dux) / 2 + ut * ux * theta / 3;
	PRECISION sty = -(t * un2 * uy) / 2 + (dtuy - dyut) / 2 - (uy * dut + ut * duy) / 2 + ut * uy * theta / 3;
	PRECISION stn = -un * (2 * ut2 + t2 * un2) / (2 * t) + (dtun - dnut / t2) / 2 - (un * dut + ut * dun) / 2 + ut * un * theta / 3;
	PRECISION sxx = -(dxux + ux * dux) + (1 + ux*ux) * theta / 3;
	PRECISION sxy = -(dxuy + dyux) / 2 - (uy * dux + ux * duy) / 2	+ ux * uy * theta / 3;
	PRECISION sxn = -ut * ux * un / t - (dxun + dnux / t2) / 2 - (un * dux + ux * dun) / 2 + ux * un * theta / 3;
	PRECISION syy = -(dyuy + uy * duy) + (1 + uy*uy) * theta / 3;
	PRECISION syn = -ut * uy * un / t - (dyun + dnuy / t2) / 2 - (un * duy + uy * dun) / 2 + uy * un * theta / 3;
	PRECISION snn = -ut * (1 + 2 * t2 * un2) / t3 - dnun / t2 - un * dun + (1 / t2 + un2) * theta / 3;

	// vorticity tensor
	PRECISION wtx = (dtux + dxut) / 2 + (ux * dut - ut * dux) / 2 + t * un2 * ux / 2;
	PRECISION wty = (dtuy + dyut) / 2 + (uy * dut - ut * duy) / 2 + t * un2 * uy / 2;
	PRECISION wtn = (t2 * dtun + 2 * t * un + dnut) / 2 + (t2 * un * dut - ut * Dun) + t3 * un*un2 / 2;
	PRECISION wxy = (dyux - dxuy) / 2 + (uy * dux - ux * duy) / 2;
	PRECISION wxn = (dnux - t2 * dxun) / 2 + (t2 * un * dux - ux * Dun) / 2;
	PRECISION wyn = (dnuy - t2 * dyun) / 2 + (t2 * un * duy - uy * Dun) / 2;
	// anti-symmetric vorticity components
	PRECISION wxt = wtx;
	PRECISION wyt = wty;
	PRECISION wnt = wtn / t2;
	PRECISION wyx = -wxy;
	PRECISION wnx = -wxn / t2;
	PRECISION wny = -wyn / t2;

	/************************************************************************************\
	 * Split transverse and longitudinal gradients of the fluid velocity
	/************************************************************************************/
	// transverse flow velocity
	PRECISION uT2 = ux*ux+uy*uy;
	PRECISION uT = sqrt(uT2);
	PRECISION uTdtuT = (ux*dtux+uy*dtuy);
	PRECISION uTdxuT = (ux*dxux+uy*dxuy);
	PRECISION uTdyuT = (ux*dyux+uy*dyuy);
	PRECISION uTdnuT = (ux*dnux+uy*dnuy);
	PRECISION uTduT = ut*uTdtuT+ux*uTdxuT+uy*uTdyuT+un*uTdnuT;

	PRECISION F = 1+uT2;
	PRECISION F2 = F*F;
	PRECISION FS = sqrt(1+uT2);

	double z0 = t*un/sqrt(F);
	double z1 = 0;
	double z2 = 0;
	double z3 = ut/t/sqrt(F);
	double z02 = z0*z0;
	double z0z3 = z0*z3;
	double z32 = z3*z3;

	// g^{\mu\nu}
	double gtt = 1;
	double gxx = -1;
	double gyy = -1;
	double gnn = -1/t2;
	// \Delta_{\perp}^{\mu\nu}
	double Xitt = gtt-ut*ut+z02;
	double Xitx = -ut*ux;
	double Xity = -ut*uy;
	double Xitn = -ut*un+z0z3;
	double Xixx = gxx-ux*ux;
	double Xixy = -ux*uy;
	double Xixn = -ux*un;
	double Xiyy = gyy-uy*un;
	double Xiyn = -uy*un;
	double Xinn = gnn-un*un+z32;

	// covariant derivatives of z
	double dz0 = (ut*un+t*dun-t*un/F*uTduT)/sqrt(F);
	double dz3 = (-ut2/t2+dut/t-ut/t/F*uTduT)/sqrt(F);
	double Dz0 = dz0+t*un*z3;
	double Dz3Upper = dz3+(ut*z3+un*z0)/t;
	double Dz3 = -t2*Dz3Upper;

	double A = z02*stt-2*t2*z0z3*stn+t2*t2*z32*snn;
	// B
	double B0 = z0*stt-t2*z3*stn;
	double B1 = z0*stx-t2*z3*sxn;
	double B2 = z0*sty-t2*z3*syn;
	double B3 = z0*stn-t2*z3*snn;
	// Bw
	double Bw0 = z3*wtn;
	double Bw1 = z0*wtx+z3*wxn;
	double Bw2 = z0*wty+z3*wyn;
	double Bw3 = z0*wtn/t2;

	// transverse and longitudinal expansion scalars
	double thetaT = 2*theta/3+A;
	double zDzu = theta/3-A;
	zDzu = theta-thetaT;

	// transverse shear tensor
	double sttT = stt+2*z0*B0+z02*A-0.5*(1-ut*ut+z02)*A;
	double stxT = stx+z0*B1-0.5*(-ut*ux)*A;
	double styT = sty+z0*B2-0.5*(-ut*uy)*A;
	double stnT = stn+z0*B3+z3*B0+z0z3*A-0.5*(-ut*un+z0z3)*A;
	double sxxT = sxx-0.5*(-1-ux*ux)*A;
	double sxyT = sxy-0.5*(-ux*uy)*A;
	double sxnT = sxn+z3*B1-0.5*(-ux*un)*A;
	double syyT = syy-0.5*(-1-uy*uy)*A;
	double synT = syn+z3*B2-0.5*(-uy*un)*A;
	double snnT = snn+2*z3*B3+z32*A-0.5*(-1/t2-un*un+z32)*A;

	/************************************************************************************\
	 * Anisotropic hydro stuff
	/************************************************************************************/
	if(e==0) e=1.e-7;
	PRECISION T = effectiveTemperature(e);
	PRECISION cs2 = speedOfSoundSquared(e);

	double a = pl/e;
	double a2 = a*a;
	double a3 = a2*a;
	double a4 = a3*a;
	double a5 = a4*a;
	double a6 = a5*a;
	double a7 = a6*a;
	double a8 = a7*a;
	double a9 = a8*a;
	double a10 = a9*a;
	double a11 = a10*a;
	double a12 = a11*a;
	double a13 = a12*a;
	double a14 = a13*a;
	double a15 = a14*a;

	PRECISION Rtilde = (-6.674731906076046e-6 + 0.004617789933500251*a + 0.7207562721999754*a2 + 9.097427250602184*a3 - 4.475814747302824*a4 - 36.37501529319408*a5 +
     46.868405146729316*a6 - 15.833867583743228*a7)/
   (0.06856675185266 + 2.9181587012768597*a + 11.951184087839218*a2 - 29.708257843442173*a3 - 2.618233802059826*a4 + 34.646239784689065*a5 -
     19.62596366454439*a6 + 2.374808442453899*a7);

	PRECISION Rhat = (0.0024792827625583747 + 1.943027171680747*a + 53.46970495217282*a2 + 19.989171951866325*a3 - 347.1285593126723*a4 + 412.2647882672885*a5 -
     140.53693383827797*a6)/(0.5061402347582388 + 29.466067530916984*a + 126.07947638942892*a2 - 334.420268508072*a3 + 86.57706367583984*a4 +
     183.53625188578846*a5 - 91.68259808111912*a6);

	PRECISION Rbar0 = (0.015267955823446243 + 7.725572805021035*a + 421.0063884634789*a2 + 3422.877939650926*a3 - 5785.670846299543*a4 - 12261.66452089229*a5 +
     31491.409484673808*a6 - 22737.05146992673*a7 + 5441.373392185447*a8)/
   (0.05470696094814806 + 14.505878005231883*a + 522.6643024173569*a2 + 2731.7776413939037*a3 - 6161.1991042880445*a4 -
     3989.4375208972588*a5 + 15008.260526258282*a6 - 10243.036679405379*a7 + 2116.74060159494*a8);

	PRECISION Rgamma = (0.0001373796340585521 - 0.6149634004722385*a + 3.1968875683514253*a2 + 246.35973783799196*a3 - 764.319750320186*a4 +
     834.160165071214*a5 - 371.64955466673234*a6 + 52.87411921963021*a7)/
   (-1.673322188488071 + 13.343782941396997*a + 561.7566534476223*a2 - 1790.2296622275915*a3 + 1896.4688704912812*a4 -
     658.7933063369629*a5 - 85.96181900698849*a6 + 65.09739194472589*a7);

	double Rbar0P = (2.473173363908116e-10 - 4.899839370307281e-6*a + 155055.91462124084*a10 - 275435.45350226434*a11 + 350689.68825705117*a12 -
     299725.38986957155*a13 + 151477.08809203724*a14 - 33196.47417939176*a15 - 0.004301975027942015*a2 - 0.14858206981041563*a3 +
     6.249255189587875*a4 - 92.79641927240235*a5 + 807.175057749925*a6 - 4760.015905266286*a7 + 20324.533122685436*a8 -
     64758.869552496515*a9)/(0.00008222793468208523 + 0.03411917870833943*a + 4.895969276094396e6*a10 - 8.84162305829353e6*a11 +
     1.1445063656613324e7*a12 - 9.918442713390596e6*a13 + 5.065882388219598e6*a14 - 1.1181016364928822e6*a15 + 0.23871740573818725*a2 -
     23.50912574236691*a3 + 417.4953123877312*a4 - 4234.215775452717*a5 + 29824.022790048104*a6 - 157419.8447785501*a7 +
     641300.6529027821*a8 - 2.0248032895288002e6*a9);

	double R0g1m2031 = (-1.0233483506421896e-7 - 0.005510394233958543*a - 0.5161308003737349*a2 - 4.115511461930346*a3 + 6.378431203946746*a4 + 3.926438723664259*a5 -
     8.465485699618803*a6 + 2.7925630611642154*a7)/
   (0.001958161993306958 + 0.22517859370360388*a + 2.883216830325076*a2 + 0.1905363935371778*a3 - 12.192584184275201*a4 + 10.729468548753893*a5 -
     0.8635431725599291*a6 - 0.9690254375998808*a7);

	// longitudinal pressure
	PRECISION dxpl = (*(currentVars->pl + s + 1) - *(currentVars->pl + s - 1)) * facX;
	PRECISION dypl = (*(currentVars->pl + s + d_ncx) - *(currentVars->pl + s - d_ncx)) * facY;
	PRECISION dnpl = (*(currentVars->pl + s + stride) - *(currentVars->pl + s - stride)) * facZ;

	PRECISION ptHat = transversePressureHat(e, p, pl);
	PRECISION pt = ptHat + 1.5*Pi;
	PRECISION DP = pl-pt;
	if(isnan(pl)) printf("found pl nan \n");
	// L functions
	PRECISION Ltt = DP*t2*un2/F;
	PRECISION Ltn = DP*ut*un/F;
	PRECISION Lnn = DP*ut2/F/t2;
	// W functions
	double Wtt = 2*WtTz*z0;
	double Wtx = WxTz*z0;
	double Wty = WyTz*z0;
	double Wtn = WtTz*z3+WnTz*z0;
	double Wxx = 0;
	double Wxy = 0;
	double Wxn = WxTz*z3;
	double Wyy = 0;
	double Wyn = WyTz*z3;
	double Wnn = 2*WnTz*z3;

	// derivatives of z0
	double dxz0 = t*(dxun/sqrt(F)-un/pow(F,1.5)*uTdxuT);
	double dyz0 = t*(dyun/sqrt(F)-un/pow(F,1.5)*uTdyuT);
	double dnz0 = t*(dnun/sqrt(F)-un/pow(F,1.5)*uTdnuT);
	// derivatives of z3
	double dxz3 = (dxut/sqrt(F)-ut/pow(F,1.5)*uTdxuT)/t;
	double dyz3 = (dyut/sqrt(F)-ut/pow(F,1.5)*uTdyuT)/t;
	double dnz3 = (dnut/sqrt(F)-ut/pow(F,1.5)*uTdnuT)/t;

	// derivative of W^{\mu}_{\perp z}
#ifdef W_TZ_MU
	// WtTz
	PRECISION dxWtTz = (*(currentVars->WtTz + s + 1) - *(currentVars->WtTz + s - 1)) * facX;
	PRECISION dyWtTz = (*(currentVars->WtTz + s + d_ncx) - *(currentVars->WtTz + s - d_ncx)) * facY;
	PRECISION dnWtTz = (*(currentVars->WtTz + s + stride) - *(currentVars->WtTz + s - stride)) * facZ;
	// WxTz
	PRECISION dxWxTz = (*(currentVars->WxTz + s + 1) - *(currentVars->WxTz + s - 1)) * facX;
	PRECISION dyWxTz = (*(currentVars->WxTz + s + d_ncx) - *(currentVars->WxTz + s - d_ncx)) * facY;
	PRECISION dnWxTz = (*(currentVars->WxTz + s + stride) - *(currentVars->WxTz + s - stride)) * facZ;
	// WyTz
	PRECISION dxWyTz = (*(currentVars->WyTz + s + 1) - *(currentVars->WyTz + s - 1)) * facX;
	PRECISION dyWyTz = (*(currentVars->WyTz + s + d_ncx) - *(currentVars->WyTz + s - d_ncx)) * facY;
	PRECISION dnWyTz = (*(currentVars->WyTz + s + stride) - *(currentVars->WyTz + s - stride)) * facZ;
	// WnTz
	PRECISION dxWnTz = (*(currentVars->WnTz + s + 1) - *(currentVars->WnTz + s - 1)) * facX;
	PRECISION dyWnTz = (*(currentVars->WnTz + s + d_ncx) - *(currentVars->WnTz + s - d_ncx)) * facY;
	PRECISION dnWnTz = (*(currentVars->WnTz + s + stride) - *(currentVars->WnTz + s - stride)) * facZ;
#else
	// WnTz
	PRECISION dxWtTz = 0;
	PRECISION dyWtTz = 0;
	PRECISION dnWtTz = 0;
	// WxTz
	PRECISION dxWxTz = 0;
	PRECISION dyWxTz = 0;
	PRECISION dnWxTz = 0;
	// WyTz
	PRECISION dxWyTz = 0;
	PRECISION dyWyTz = 0;
	PRECISION dnWyTz = 0;
	// WnTz
	PRECISION dxWnTz = 0;
	PRECISION dyWnTz = 0;
	PRECISION dnWnTz = 0;
#endif

	//=========================================================
	// derivatives of W in conservation law source terms
	//=========================================================
	// Wtt
	double dxWtt = 2*t*(WtTz*dxun-un*WtTz*uTdxuT/F+un*dxWtTz)/FS;
	double dyWtt = 2*t*(WtTz*dyun-un*WtTz*uTdyuT/F+un*dyWtTz)/FS;
	double dnWtt = 2*t*(WtTz*dnun-un*WtTz*uTdnuT/F+un*dnWtTz)/FS;
	// Wtx
	double dxWtx = t*(WxTz*dxun-un*WxTz*uTdxuT/F+un*dxWxTz)/FS;
	double dyWtx = t*(WxTz*dyun-un*WxTz*uTdyuT/F+un*dyWxTz)/FS;
	double dnWtx = t*(WxTz*dnun-un*WxTz*uTdnuT/F+un*dnWxTz)/FS;
	// Wty
	double dxWty = t*(WyTz*dxun-un*WyTz*uTdxuT/F+un*dxWyTz)/FS;
	double dyWty = t*(WyTz*dyun-un*WyTz*uTdyuT/F+un*dyWyTz)/FS;
	double dnWty = t*(WyTz*dnun-un*WyTz*uTdnuT/F+un*dnWyTz)/FS;
	// Wtn
	double dxWtn = (WtTz*dxut-ut*WtTz*uTdxuT/F+ut*dxWyTz)/FS/t+t*(WnTz*dxun-un*WnTz*uTdxuT/F+un*dxWnTz)/FS;
	double dyWtn = (WtTz*dyut-ut*WtTz*uTdyuT/F+ut*dyWyTz)/FS/t+t*(WnTz*dyun-un*WnTz*uTdyuT/F+un*dyWnTz)/FS;
	double dnWtn = (WtTz*dnut-ut*WtTz*uTdnuT/F+ut*dnWyTz)/FS/t+t*(WnTz*dnun-un*WnTz*uTdnuT/F+un*dnWnTz)/FS;
	// Wxn
	double dxWxn = (WxTz*dxut-ut*WxTz*uTdxuT/F+ut*dxWxTz)/FS/t;
	double dnWxn = (WxTz*dnut-ut*WxTz*uTdnuT/F+ut*dnWxTz)/FS/t;
	// Wyn
	double dyWyn = (WyTz*dyut-ut*WyTz*uTdyuT/F+ut*dyWyTz)/FS/t;
	double dnWyn = (WyTz*dnut-ut*WyTz*uTdnuT/F+ut*dnWyTz)/FS/t;
	// Wnn
	double dnWnn = 2*(WnTz*dnut-ut*WnTz*uTdnuT/F+ut*dnWnTz)/FS/t;

	PRECISION tnn = (e+pt)*un*un+pt/t2+Lnn+Wnn+pinn;

	// deritive of Rbar
	double dxRbar0 = Rbar0P*(dxpl-a*dxe)/e;
	double dyRbar0 = Rbar0P*(dypl-a*dye)/e;
	double dnRbar0 = Rbar0P*(dnpl-a*dne)/e;
	// derivative of transverse pressure
	PRECISION dxptHat = 0.5*(dxe-dxpl-Rbar0*(dxe-3*dxp)-(e-3*p)*dxRbar0);
	if (isnan(e)) printf("found e nan \n");
	if (isnan(p)) printf("found p nan \n");
	if (isnan(dxe)) printf("found dxe nan \n");
	if (isnan(dxp)) printf("found dxp nan \n");
	if (isnan(dxpl)) printf("found dxpl nan \n");
	if (isnan(dye)) printf("found dye nan \n");
	if (isnan(dyp)) printf("found dyp nan \n");
	if (isnan(dypl)) printf("found dypl nan \n");
	if (isnan(Rbar0)) printf("found Rbar0 nan \n");
	if (isnan(dxRbar0)) printf("found dxRbar0 nan \n");
	if (isnan(dyRbar0)) printf("found dyRbar0 nan \n");
	if (isnan(dnRbar0)) printf("found dnRbar0 nan \n");
	PRECISION dyptHat = 0.5*(dye-dypl-Rbar0*(dye-3*dyp)-(e-3*p)*dyRbar0);
	PRECISION dnptHat = 0.5*(dne-dnpl-Rbar0*(dne-3*dnp)-(e-3*p)*dnRbar0);

	// derivative of \Delta P -- NEED TO ADD DERIVATIVES OF BULK PI
	PRECISION dxDP = dxpl-dxptHat;
	PRECISION dyDP = dypl-dyptHat;
	PRECISION dnDP = dnpl-dnptHat;

	// spatial derivatives of Ltt
	PRECISION dxLtt = (2*DP*dxun*t2*un)/F + (dxDP*t2*un2)/F - (2*DP*uTdxuT*t2*un2)/F2;
	PRECISION dyLtt = (2*DP*dyun*t2*un)/F + (dyDP*t2*un2)/F - (2*DP*uTdyuT*t2*un2)/F2;
	PRECISION dnLtt = (2*DP*dnun*t2*un)/F + (dnDP*t2*un2)/F - (2*DP*uTdnuT*t2*un2)/F2;
	// spatial derivatives of Ltn
	PRECISION dxLtn = (DP*dxut*un)/F + (DP*dxun*ut)/F + (dxDP*un*ut)/F - (2*DP*uTdxuT*un*ut)/F2;
	PRECISION dyLtn = (DP*dyut*un)/F + (DP*dyun*ut)/F + (dyDP*un*ut)/F - (2*DP*uTdyuT*un*ut)/F2;
	PRECISION dnLtn = (DP*dnut*un)/F + (DP*dnun*ut)/F + (dnDP*un*ut)/F - (2*DP*uTdnuT*un*ut)/F2;
	// spatial derivatives of Lnn
	PRECISION dnLnn = (2*DP*dnut*ut)/(F*t2) + (dnDP*ut2)/(F*t2) - (2*DP*uTdnuT*ut2)/(F2*t2);

	PRECISION zeta_zz = Rtilde*e-3*pl;
	PRECISION zeta_zT = (Rtilde*e + pl)/2.;

	//==============================================
	// second-order transport coefficients
	//==============================================
	double beta_lPi=0;
	double delta_lPi=0;
	double lambda_piPi=0;
	double beta_PiPi=0;
	double delta_PiPi=0;
	double lambda_Pipi=0;
	secondOrderTransportCoefficientsZ(e, p, pl, cs2, T, &beta_lPi, &delta_lPi, &lambda_piPi, &beta_PiPi, &delta_PiPi, &lambda_Pipi);
	// Pl
	double lambda_lWu = R0g1m2031;
	double lambda_lWT = 2+lambda_lWu;
	double lambda_lpi = Rgamma;
	// W
	double delta_WW = lambda_lWu/2-1;
	double lambda_WWu = 2 + lambda_lWu;
	double lambda_WWT = 1 + delta_WW;
	double lambda_Wpiu = lambda_lpi;
	double lambda_WpiT = lambda_Wpiu-1;
	// pi^\mu\nu
	double delta_pipi = (3+Rgamma)/2;
	double tau_pipi = 4*(delta_pipi-1/2)/3;
	double lambda_pipi = Rgamma-1;
	double lambda_piWu = delta_WW-1/2;
	double lambda_piWT = lambda_piWu+2;

	PRECISION taupiInv = 0.2 * T/d_etabar;

	PRECISION psT = pitt * sttT - 2 * pitx * stxT - 2 * pity * styT + pixx * sxxT + 2 * pixy * sxyT + piyy * syyT - 2 * pitn * stnT * t2 + 2 * pixn * sxnT * t2
			+ 2 * piyn * synT * t2 + pinn * snnT * t2 * t2;

	// IL1
	double IL1 = -(WtTz*B0 - WxTz*B1 - WyTz*B2 - t2*WnTz*B3) + (WtTz*Bw0 - WxTz*Bw1 - WyTz*Bw2 - t2*WnTz*Bw3);
	// IL2
	double IL2 = (WtTz*B0 - WxTz*B1 - WyTz*B2 - t2*WnTz*B3) + (WtTz*Bw0 - WxTz*Bw1 - WyTz*Bw2 - t2*WnTz*Bw3);
	// IL3
	double IL3 = psT;
	// IL
	double IL = lambda_lWu * IL1 - lambda_lWT * IL2 - lambda_lpi * IL3;
//	IL=0;

//	PRECISION dPL = -(pl-p)*taupiInv + zeta_zz*zDzu - zeta_zT*thetaT - lambda_lpi * psT;
	PRECISION dPL = -(pl-p)*taupiInv + zeta_zz*zDzu - zeta_zT*thetaT - beta_lPi * Pi * zDzu - delta_lPi * Pi * thetaT + IL;
	S[4] = dPL / ut + dkvk * pl;

	if(isnan(S[4])) {
		printf("=======================================================================================\n");
		printf("Found Nan in S[4]:\n");
		printf("pl=%.9f;\n",pl);
		printf("Grid point = (%d, %d, %d) = (%.3f, %.3f, %.3f)\n", i, j, k, x, y, z);
		printf("Rhat=%.9f;Rtilde=%.39f;\n",Rhat,Rtilde);
		printf("pl = %.9f\n",pl);
		printf("e = %.9f\n",e);
		printf("a = %.9f\n",a);
		printf("T=%.9f;e=%.9f;p=%.9f;\n",T,e,p);
		printf("ut=%.9f;ux=%.9f;uy=%.9f\n",ut,ux,uy);
		printf("ux[s+1]=%.9f; ux[s-1]=%.9f\n",*(uxvec + s + 1),*(uxvec + s - 1));
		printf("uy[s+1]=%.9f; uy[s-1]=%.9f\n",*(uyvec + s + 1),*(uyvec + s - 1));
		printf("dx=%.3f;dy=%.3f;facX=%.3f;facY=%.3f\n",d_dx,d_dy,facX,facY);
		printf("zeta_zz=%.3f;zeta_zT=%.3f;dtu0=%.3f;zDzu=%.3f;thetaT=%.3f;taupiInv=%.3f;\n",zeta_zz,zeta_zT,dtut,zDzu,thetaT,taupiInv);
		printf("=======================================================================================\n");
		exit(-1);
	}

	/************************************************************************************\
	 * T^{\mu\nu} source terms
	/************************************************************************************/
	double Htt = Ltt+Wtt+pitt;
	S[0] = -(ttt / t + t * tnn) + dkvk * (Htt - pt) - vx * dxptHat - vy * dyptHat - vn * dnptHat
				+ vx*(dxLtt+dxWtt)-dxWtx + vy*(dyLtt+dyWtt)-dyWty + vn*(dnLtt+dnWtt)-dnLtn-dnWtn;
	S[1] = -ttx / t - dxptHat + dkvk * (Wtx + pitx) + vx*dxWtx + vy*dyWtx + vn*dnWtx - dnWxn;
	S[2] = -tty / t - dyptHat + dkvk * (Wty + pity) + vx*dxWty + vy*dyWty + vn*dnWty - dnWyn;
	S[3] = -3 * ttn / t - dnptHat / t2 + dkvk * (Ltn + Wtn + pitn)
				+ vx*(dxLtn+dxWtn)-dxWxn + vy*(dyLtn+dyWtn)-dyWyn + vn*(dnLtn+dnLtn)-dnLnn-dnWnn;

	if(isnan(S[0])) {
		printf("=======================================================================================\n");
		printf("Found Nan in S[0]:\n");
		printf("pl=%.9f;\n",pl);
		printf("Grid point = (%d, %d, %d) = (%.3f, %.3f, %.3f)\n", i, j, k, x, y, z);
		printf("Rhat=%.9f;Rtilde=%.39f;\n",Rhat,Rtilde);
		printf("pl = %.9f\n",pl);
		printf("e = %.9f\n",e);
		printf("a = %.9f\n",a);
		printf("T=%.9f;e=%.9f;p=%.9f;\n",T,e,p);
		printf("ut=%.9f;ux=%.9f;uy=%.9f\n",ut,ux,uy);
		printf("ux[s+1]=%.9f; ux[s-1]=%.9f\n",*(uxvec + s + 1),*(uxvec + s - 1));
		printf("uy[s+1]=%.9f; uy[s-1]=%.9f\n",*(uyvec + s + 1),*(uyvec + s - 1));
		printf("dx=%.3f;dy=%.3f;facX=%.3f;facY=%.3f\n",d_dx,d_dy,facX,facY);
		printf("zeta_zz=%.3f;zeta_zT=%.3f;dtu0=%.3f;zDzu=%.3f;thetaT=%.3f;taupiInv=%.3f;\n",zeta_zz,zeta_zT,dtut,zDzu,thetaT,taupiInv);
		printf("taupiInv=%.3f;\n",taupiInv);
		printf("-------------------------------------------\n");
		printf("ttt=%.9f;tnn=%.9f\n",ttt,tnn);
		printf("dkvk=%.9f;vx=%.9f;vy=%.9f;vn=%.9f\n",dkvk,vx,vy,vn);
		printf("pt=%.9f;dxpt=%.9f;dypt=%.9f;dnpt=%.9f\n",pt,dxptHat,dyptHat,dnptHat);
		printf("Ltt=%.9f;dxLtt=%.9f;dyLtt=%.9f;dnLtt=%.9f\n",Ltt,dxLtt,dyLtt,dnLtt);
		printf("Ltn=%.9f;dxLtn=%.9f;dyLtn=%.9f;dnLtn=%.9f\n",Ltn,dxLtn,dyLtn,dnLtn);
		printf("dnLnn=%.9f\n",dnLnn);
		printf("dxpt=%.9f;dypt=%.9f;dnpt=%.9f;\n",dxptHat,dyptHat,dnptHat);
		printf("dxDP=%.9f;dyDP=%.9f;dnDP=%.9f;\n",dxDP,dyDP,dnDP);
		printf("DP=%.9f;F=%.9f;uT=%.9f\n",DP,F,uT);
		printf("uTdxuT=%.9f;uTdyuT=%.9f;uTdnuT=%.9f;\n",uTdxuT,uTdyuT,uTdnuT);
		printf("dxux=%.9f;dyux=%.9f;dnux=%.9f;\n",dxux,dyux,dnux);
		printf("dxuy=%.9f;dyuy=%.9f;dnuy=%.9f;\n",dxuy,dyuy,dnuy);
		printf("=======================================================================================\n");
		exit(-1);
	}

	/************************************************************************************\
	 * \pi^{\mu\nu}_\perp source terms
	/************************************************************************************/
#ifdef PIMUNU
	double etaBar_TC = ((3-Rtilde)*e-2*pl)/8;

	// I1
	PRECISION I1Utt = 2 * ut * (pitt * Dut + pitx * Dux + pity * Duy + pitn * Dun);
	PRECISION I1Utx = (pitt * ux + pitx * ut) * Dut + (pitx * ux + pixx * ut) * Dux + (pity * ux + pixy * ut) * Duy + (pitn * ux + pixn * ut) * Dun;
	PRECISION I1Uty = (pitt * uy + pity * ut) * Dut + (pitx * uy + pixy * ut) * Dux + (pity * uy + piyy * ut) * Duy + (pitn * uy + piyn * ut) * Dun;
	PRECISION I1Utn = (pitt * un + pitn * ut) * Dut + (pitx * un + pixn * ut) * Dux + (pity * un + piyn * ut) * Duy + (pitn * un + pinn * ut) * Dun;
	PRECISION I1Uxx = 2 * ux * (pitx * Dut + pixx * Dux + pixy * Duy + pixn * Dun);
	PRECISION I1Uxy = (pitx * uy + pity * ux) * Dut + (pixx * uy + pixy * ux) * Dux + (pixy * uy + piyy * ux) * Duy + (pixn * uy + piyn * ux) * Dun;
	PRECISION I1Uxn = (pitx * un + pitn * ux) * Dut + (pixx * un + pixn * ux) * Dux + (pixy * un + piyn * ux) * Duy + (pixn * un + pinn * ux) * Dun;
	PRECISION I1Uyy = 2 * uy * (pity * Dut + pixy * Dux + piyy * Duy + piyn * Dun);
	PRECISION I1Uyn = (pity * un + pitn * uy) * Dut + (pixy * un + pixn * uy) * Dux + (piyy * un + piyn * uy) * Duy + (piyn * un + pinn * uy) * Dun;
	PRECISION I1Unn = 2 * un * (pitn * Dut + pixn * Dux + piyn * Duy + pinn * Dun);

//Dz0=0;
//Dz3=0;

	PRECISION I1Ztt = 2 * z0 * (pitt * Dz0 + pitn * Dz3);
	PRECISION I1Ztx = ( pitx * z0) * Dz0 + (pitn * ux + pixn * z0) * Dz3;
	PRECISION I1Zty = (pity * z0) * Dz0 + (pitn * uy + piyn * z0) * Dz3;
	PRECISION I1Ztn = (pitt * z3 + pitn * z0) * Dz0 + (pitn * z3 + pinn * z0) * Dz3;
	PRECISION I1Zxx = 0;
	PRECISION I1Zxy = 0;
	PRECISION I1Zxn = (pitx * z3) * Dz0 + (pixn * z3) * Dz3;
	PRECISION I1Zyy = 0;
	PRECISION I1Zyn = (pity * z3) * Dz0 + (piyn * z3) * Dz3;
	PRECISION I1Znn = 2 * z3 * (pitn * Dz0 + pinn * Dz3);

	PRECISION I1tt = I1Utt - I1Ztt;
	PRECISION I1tx = I1Utx - I1Ztx;
	PRECISION I1ty = I1Uty - I1Zty;
	PRECISION I1tn = I1Utn - I1Ztn;
	PRECISION I1xx = I1Uxx - I1Zxx;
	PRECISION I1xy = I1Uxy - I1Zxy;
	PRECISION I1xn = I1Uxn - I1Zxn;
	PRECISION I1yy = I1Uyy - I1Zyy;
	PRECISION I1yn = I1Uyn - I1Zyn;
	PRECISION I1nn = I1Unn - I1Znn;

	// I2
	PRECISION I2tt = thetaT * pitt;
	PRECISION I2tx = thetaT * pitx;
	PRECISION I2ty = thetaT * pity;
	PRECISION I2tn = thetaT * pitn;
	PRECISION I2xx = thetaT * pixx;
	PRECISION I2xy = thetaT * pixy;
	PRECISION I2xn = thetaT * pixn;
	PRECISION I2yy = thetaT * piyy;
	PRECISION I2yn = thetaT * piyn;
	PRECISION I2nn = thetaT * pinn;

	// I4
///*
	PRECISION ux2 = ux * ux;
	PRECISION uy2 = uy * uy;
	PRECISION ps = pitt * sttT - 2 * pitx * stxT - 2 * pity * styT + pixx * sxxT + 2 * pixy * sxyT + piyy * syyT - 2 * pitn * stnT * t2 + 2 * pixn * sxnT * t2
			+ 2 * piyn * synT * t2 + pinn * snnT * t2 * t2;
	PRECISION ps2 = ps / 2;
	PRECISION I4tt = (pitt * sttT - pitx * stxT - pity * styT - t2 * pitn * stnT) - (1 - ut2+z02) * ps2;
	PRECISION I4tx = (pitt * stxT + pitx * sttT) / 2 - (pitx * sxxT + pixx * stxT) / 2 - (pity * sxyT + pixy * styT) / 2 - t2 * (pitn * sxnT + pixn * stnT) / 2
			+ (ut * ux) * ps2;
	PRECISION I4ty = (pitt * styT + pity * sttT) / 2 - (pitx * sxyT + pixy * stxT) / 2 - (pity * syyT + piyy * styT) / 2 - t2 * (pitn * synT + piyn * stnT) / 2
			+ (ut * uy) * ps2;
	PRECISION I4tn = (pitt * stnT + pitn * sttT) / 2 - (pitx * sxnT + pixn * stxT) / 2 - (pity * synT + piyn * styT) / 2 - t2 * (pitn * snnT + pinn * stnT) / 2
			+ (ut * un-z0z3) * ps2;
	PRECISION I4xx = (pitx * stxT - pixx * sxxT - pixy * sxyT - t2 * pixn * sxnT) + (1 + ux2) * ps2;
	PRECISION I4xy = (pitx * styT + pity * stxT) / 2 - (pixx * sxyT + pixy * sxxT) / 2 - (pixy * syyT + piyy * sxyT) / 2 - t2 * (pixn * synT + piyn * sxnT) / 2
			+ (ux * uy) * ps2;
	PRECISION I4xn = (pitx * stnT + pitn * stxT) / 2 - (pixx * sxnT + pixn * sxxT) / 2 - (pixy * synT + piyn * sxyT) / 2 - t2 * (pixn * snnT + pinn * sxnT) / 2
			+ (ux * un) * ps2;
	PRECISION I4yy = (pity * styT - pixy * sxyT - piyy * syyT - t2 * piyn * synT) + (1 + uy2) * ps2;
	PRECISION I4yn = (pity * stnT + pitn * styT) / 2 - (pixy * sxnT + pixn * sxyT) / 2 - (piyy * synT + piyn * syyT) / 2 - t2 * (piyn * snnT + pinn * synT) / 2
			+ (uy * un) * ps2;
	PRECISION I4nn = (pitn * stnT - pixn * sxnT - piyn * synT - t2 * pinn * snnT) + (1 / t2 + un2-z32) * ps2;
//*/
/*
	 PRECISION I4tt = 0;
	 PRECISION I4tx = 0;
	 PRECISION I4ty = 0;
	 PRECISION I4tn = 0;
	 PRECISION I4xx = 0;
	 PRECISION I4xy = 0;
	 PRECISION I4xn = 0;
	 PRECISION I4yy = 0;
	 PRECISION I4yn = 0;
	 PRECISION I4nn = 0;
//*/

	// I5
	PRECISION I5tt = zDzu * pitt;
	PRECISION I5tx = zDzu * pitx;
	PRECISION I5ty = zDzu * pity;
	PRECISION I5tn = zDzu * pitn;
	PRECISION I5xx = zDzu * pixx;
	PRECISION I5xy = zDzu * pixy;
	PRECISION I5xn = zDzu * pixn;
	PRECISION I5yy = zDzu * piyy;
	PRECISION I5yn = zDzu * piyn;
	PRECISION I5nn = zDzu * pinn;

	// I6
	double L = z0*WtTz - t2*z3*WnTz;
	double WB = WtTz * B0 - WxTz * B1 - WyTz * B2 - t2 * WnTz * B3;
	// I6S
	double I6tt_S = ((WtTz*B0 + WtTz*B0) + (WtTz*z0 + z0*WtTz)*A + (z0*B0 + z0*B0)*L + 2*z0*z0*L*A)/2 - Xitt*(WB + L*A)/2;
	double I6tx_S = ((WtTz*B1 + WxTz*B0) + (WtTz*z1 + z0*WxTz)*A + (z0*B1 + z1*B0)*L + 2*z0*z1*L*A)/2 - Xitx*(WB + L*A)/2;
	double I6ty_S = ((WtTz*B2 + WyTz*B0) + (WtTz*z2 + z0*WyTz)*A + (z0*B2 + z2*B0)*L + 2*z0*z2*L*A)/2 - Xity*(WB + L*A)/2;
	double I6tn_S = ((WtTz*B3 + WnTz*B0) + (WtTz*z3 + z0*WnTz)*A + (z0*B3 + z3*B0)*L + 2*z0*z3*L*A)/2 - Xitn*(WB + L*A)/2;
	double I6xx_S = ((WxTz*B1 + WxTz*B1) + (WxTz*z1 + z1*WxTz)*A + (z1*B1 + z1*B1)*L + 2*z1*z1*L*A)/2 - Xixx*(WB + L*A)/2;
	double I6xy_S = ((WxTz*B2 + WyTz*B1) + (WxTz*z2 + z1*WyTz)*A + (z2*B1 + z1*B2)*L + 2*z1*z2*L*A)/2 - Xixy*(WB + L*A)/2;
	double I6xn_S = ((WxTz*B3 + WnTz*B1) + (WxTz*z3 + z1*WnTz)*A + (z1*B3 + z3*B1)*L + 2*z1*z3*L*A)/2 - Xixn*(WB + L*A)/2;
	double I6yy_S = ((WyTz*B2 + WyTz*B2) + (WyTz*z2 + z2*WyTz)*A + (z2*B2 + z2*B2)*L + 2*z2*z2*L*A)/2 - Xiyy*(WB + L*A)/2;
	double I6yn_S = ((WyTz*B3 + WnTz*B2) + (WyTz*z3 + z2*WnTz)*A + (z2*B3 + z3*B2)*L + 2*z2*z3*L*A)/2 - Xiyn*(WB + L*A)/2;
	double I6nn_S = ((WnTz*B3 + WnTz*B3) + (WnTz*z3 + z3*WnTz)*A + (z3*B3 + z3*B3)*L + 2*z3*z3*L*A)/2 - Xinn*(WB + L*A)/2;
	// I6W
	double I6tt_W = ((WtTz*Bw0 + WtTz*Bw0) + (z0*Bw0 + z0*Bw0)*L)/2 - Xitt*(WB)/2;
	double I6tx_W = ((WtTz*Bw1 + WxTz*Bw0) + (z0*Bw1 + z1*Bw0)*L)/2 - Xitx*(WB)/2;
	double I6ty_W = ((WtTz*Bw2 + WyTz*Bw0) + (z0*Bw2 + z2*Bw0)*L)/2 - Xity*(WB)/2;
	double I6tn_W = ((WtTz*Bw3 + WnTz*Bw0) + (z0*Bw3 + z3*Bw0)*L)/2 - Xitn*(WB)/2;
	double I6xx_W = ((WxTz*Bw1 + WxTz*Bw1) + (z1*Bw1 + z1*Bw1)*L)/2 - Xixx*(WB)/2;
	double I6xy_W = ((WxTz*Bw2 + WyTz*Bw1) + (z2*Bw1 + z1*Bw2)*L)/2 - Xixy*(WB)/2;
	double I6xn_W = ((WxTz*Bw3 + WnTz*Bw1) + (z1*Bw3 + z3*Bw1)*L)/2 - Xixn*(WB)/2;
	double I6yy_W = ((WyTz*Bw2 + WyTz*Bw2) + (z2*Bw2 + z2*Bw2)*L)/2 - Xiyy*(WB)/2;
	double I6yn_W = ((WyTz*Bw3 + WnTz*Bw2) + (z2*Bw3 + z3*Bw2)*L)/2 - Xiyn*(WB)/2;
	double I6nn_W = ((WnTz*Bw3 + WnTz*Bw3) + (z3*Bw3 + z3*Bw3)*L)/2 - Xinn*(WB)/2;
	// I6
	double I6tt = -I6tt_S + I6tt_W;
	double I6tx = -I6tx_S + I6tx_W;
	double I6ty = -I6ty_S + I6ty_W;
	double I6tn = -I6tn_S + I6tn_W;
	double I6xx = -I6xx_S + I6xx_W;
	double I6xy = -I6xy_S + I6xy_W;
	double I6xn = -I6xn_S + I6xn_W;
	double I6yy = -I6yy_S + I6yy_W;
	double I6yn = -I6yn_S + I6yn_W;
	double I6nn = -I6nn_S + I6nn_W;
/*
	I6tt = 0;
	I6tx = 0;
	I6ty = 0;
	I6tn = 0;
	I6xx = 0;
	I6xy = 0;
	I6xn = 0;
	I6yy = 0;
	I6yn = 0;
	I6nn = 0;
//*/

	// I7
	double I7tt = I6tt_S + I6tt_W;
	double I7tx = I6tx_S + I6tx_W;
	double I7ty = I6ty_S + I6ty_W;
	double I7tn = I6tn_S + I6tn_W;
	double I7xx = I6xx_S + I6xx_W;
	double I7xy = I6xy_S + I6xy_W;
	double I7xn = I6xn_S + I6xn_W;
	double I7yy = I6yy_S + I6yy_W;
	double I7yn = I6yn_S + I6yn_W;
	double I7nn = I6nn_S + I6nn_W;
/*
	I7tt = 0;
	I7tx = 0;
	I7ty = 0;
	I7tn = 0;
	I7xx = 0;
	I7xy = 0;
	I7xn = 0;
	I7yy = 0;
	I7yn = 0;
	I7nn = 0;
//*/


	PRECISION Itt = I1tt - delta_pipi * I2tt + tau_pipi * I4tt - lambda_pipi * I5tt + lambda_piWu * I6tt - lambda_piWT * I7tt - lambda_piPi * Pi * sttT;
	PRECISION Itx = I1tx - delta_pipi * I2tx + tau_pipi * I4tx - lambda_pipi * I5tx + lambda_piWu * I6tx - lambda_piWT * I7tx - lambda_piPi * Pi * sttT;
	PRECISION Ity = I1ty - delta_pipi * I2ty + tau_pipi * I4ty - lambda_pipi * I5ty + lambda_piWu * I6ty - lambda_piWT * I7ty - lambda_piPi * Pi * sttT;
	PRECISION Itn = I1tn - delta_pipi * I2tn + tau_pipi * I4tn - lambda_pipi * I5tn + lambda_piWu * I6tn - lambda_piWT * I7tn - lambda_piPi * Pi * sttT;
	PRECISION Ixx = I1xx - delta_pipi * I2xx + tau_pipi * I4xx - lambda_pipi * I5xx + lambda_piWu * I6xx - lambda_piWT * I7xx - lambda_piPi * Pi * sttT;
	PRECISION Ixy = I1xy - delta_pipi * I2xy + tau_pipi * I4xy - lambda_pipi * I5xy + lambda_piWu * I6xy - lambda_piWT * I7xy - lambda_piPi * Pi * sttT;
	PRECISION Ixn = I1xn - delta_pipi * I2xn + tau_pipi * I4xn - lambda_pipi * I5xn + lambda_piWu * I6xn - lambda_piWT * I7xn - lambda_piPi * Pi * sttT;
	PRECISION Iyy = I1yy - delta_pipi * I2yy + tau_pipi * I4yy - lambda_pipi * I5yy + lambda_piWu * I6yy - lambda_piWT * I7yy - lambda_piPi * Pi * sttT;
	PRECISION Iyn = I1yn - delta_pipi * I2yn + tau_pipi * I4yn - lambda_pipi * I5yn + lambda_piWu * I6yn - lambda_piWT * I7yn - lambda_piPi * Pi * sttT;
	PRECISION Inn = I1nn - delta_pipi * I2nn + tau_pipi * I4nn - lambda_pipi * I5nn + lambda_piWu * I6nn - lambda_piWT * I7nn - lambda_piPi * Pi * sttT;

	PRECISION dpitt = 2 * etaBar_TC * sttT - pitt * taupiInv - Itt - 2 * un * t * pitn;
	PRECISION dpitx = 2 * etaBar_TC * stxT - pitx * taupiInv - Itx - un * t * pixn;
	PRECISION dpity = 2 * etaBar_TC * styT - pity * taupiInv - Ity - un * t * piyn;
	PRECISION dpitn = 2 * etaBar_TC * stnT - pitn * taupiInv - Itn - un * t * pinn - (ut * pitn + un * pitt) / t;
	PRECISION dpixx = 2 * etaBar_TC * sxxT - pixx * taupiInv - Ixx;
	PRECISION dpixy = 2 * etaBar_TC * sxyT - pixy * taupiInv - Ixy;
	PRECISION dpixn = 2 * etaBar_TC * sxnT - pixn * taupiInv - Ixn - (ut * pixn + un * pitx) / t;
	PRECISION dpiyy = 2 * etaBar_TC * syyT - piyy * taupiInv - Iyy;
	PRECISION dpiyn = 2 * etaBar_TC * synT - piyn * taupiInv - Iyn - (ut * piyn + un * pity) / t;
	PRECISION dpinn = 2 * etaBar_TC * snnT - pinn * taupiInv - Inn - 2 * (ut * pinn + un * pitn) / t;

	S[5] = dpitt / ut + pitt * dkvk;
	S[6] = dpitx / ut + pitx * dkvk;
	S[7] = dpity / ut + pity * dkvk;
	S[8] = dpitn / ut + pitn * dkvk;
	S[9] = dpixx / ut + pixx * dkvk;
	S[10] = dpixy / ut + pixy * dkvk;
	S[11] = dpixn / ut + pixn * dkvk;
	S[12] = dpiyy / ut + piyy * dkvk;
	S[13] = dpiyn / ut + piyn * dkvk;
	S[14] = dpinn / ut + pinn * dkvk;
#endif
	/************************************************************************************\
	 * W^{\mu}_{\perp z} source terms
	/************************************************************************************/
#ifdef W_TZ_MU
	double etaWu = zeta_zT/2.;
	double etaWT = ((Rtilde+1)*e-2*pl)/4;

	// I1
	double WDU = WtTz*Dut + WxTz*Dux + WyTz*Duy + WnTz*Dun;
	double I1t = ut*WDU + (pitt-WtTz*z0-1.5*Pi*Xitt)*Dz0 + (pitn-WnTz*z0-1.5*Pi*Xitn)*Dz3;
	double I1x = ux*WDU + (pitx-1.5*Pi*Xitx)*Dz0 + (pixn-1.5*Pi*Xixn)*Dz3;
	double I1y = uy*WDU + (pity-1.5*Pi*Xity)*Dz0 + (piyn-1.5*Pi*Xiyn)*Dz3;
	double I1n = un*WDU + (pitn-WtTz*z3-1.5*Pi*Xitn)*Dz0 + (pinn-WnTz*z3-1.5*Pi*Xinn)*Dz3;

	// I2
	double I2t = WtTz*thetaT;
	double I2x = WxTz*thetaT;
	double I2y = WyTz*thetaT;
	double I2n = WnTz*thetaT;

	// I3
	double I3t = WtTz * zDzu;
	double I3x = WxTz * zDzu;
	double I3y = WyTz * zDzu;
	double I3n = WnTz * zDzu;

	// I4
	double I4t = WtTz * sttT - WxTz * stxT - WyTz * styT - t2 * WnTz * stnT;
	double I4x = WtTz * stxT - WxTz * sxxT - WyTz * sxyT - t2 * WnTz * sxnT;
	double I4y = WtTz * styT - WxTz * sxyT - WyTz * syyT - t2 * WnTz * synT;
	double I4n = WtTz * stnT - WxTz * sxnT - WyTz * synT - t2 * WnTz * snnT;

	// I5
	double I5t = 0;
	double I5x = 0;
	double I5y = 0;
	double I5n = 0;

	// I6
	double BB0 = Bw0-B0;
	double BB1 = Bw1-B1;
	double BB2 = Bw2-B2;
	double BB3 = Bw3-B3;
	double I6t = pitt*BB0 - pitx*BB1 - pity*BB2 - t2*pitn*BB3;
	double I6x = pitx*BB0 - pixx*BB1 - pixy*BB2 - t2*pixn*BB3;
	double I6y = pity*BB0 - pixy*BB1 - piyy*BB2 - t2*piyn*BB3;
	double I6n = pitn*BB0 - pixn*BB1 - piyn*BB2 - t2*pinn*BB3;

	// I7
	double St = B0 + z0*A + Bw0;
	double Sx = B1 + Bw1;
	double Sy = B2 + Bw2;
	double Sn = B3 + z3*A + Bw3;
	double I7t = pitt * St - pitx * Sx - pity * Sy - t2 * pitn * Sn;
	double I7x = pitx * St - pixx * Sx - pixy * Sy - t2 * pixn * Sn;
	double I7y = pity * St - pixy * Sx - piyy * Sy - t2 * piyn * Sn;
	double I7n = pitn * St - pixn * Sx - piyn * Sy - t2 * pinn * Sn;

	// J1
	double J1t = -B0-z0*A+Bw0;
	double J1x = -B1+Bw1;
	double J1y = -B2+Bw2;
	double J1n = -B3-z3*A+Bw3;
	// J2
	double J2t = B0+z0*A+Bw0;
	double J2x = B1+Bw1;
	double J2y = B2+Bw2;
	double J2n = B3+z3*A+Bw3;

	double Jt = 2*(etaWu * J1t - etaWT * J2t);
	double Jx = 2*(etaWu * J1x - etaWT * J2x);
	double Jy = 2*(etaWu * J1y - etaWT * J2y);
	double Jn = 2*(etaWu * J1n - etaWT * J2n);

	// I
	double IWt = I1t - delta_WW * I2t + lambda_WWu * I3t - lambda_WWT * I4t - I5t - lambda_Wpiu * I6t + lambda_WpiT * I7t;
	double IWx = I1x - delta_WW * I2x + lambda_WWu * I3x - lambda_WWT * I4x - I5x - lambda_Wpiu * I6x + lambda_WpiT * I7x;
	double IWy = I1y - delta_WW * I2y + lambda_WWu * I3y - lambda_WWT * I4y - I5y - lambda_Wpiu * I6y + lambda_WpiT * I7y;
	double IWn = I1n - delta_WW * I2n + lambda_WWu * I3n - lambda_WWT * I4n - I5n - lambda_Wpiu * I6n + lambda_WpiT * I7n;

	// G^{\mu}_{W} -- Christofel symbols from covariant derivative
	double GWt = t*un*WnTz;
	double GWn = (ut*WnTz+un*WtTz)/t;

	PRECISION dWtTz = Jt - WtTz*taupiInv - IWt - GWt;
	PRECISION dWxTz = Jx - WxTz*taupiInv - IWx;
	PRECISION dWyTz = Jy - WyTz*taupiInv - IWy;
	PRECISION dWnTz = Jn - WnTz*taupiInv - IWn - GWn;

	S[15] = dWtTz / ut + WtTz * dkvk;
	S[16] = dWxTz / ut + WxTz * dkvk;
	S[17] = dWyTz / ut + WyTz * dkvk;
	S[18] = dWnTz / ut + WnTz * dkvk;
#endif
	/************************************************************************************\
	 * \Pi source terms
	/************************************************************************************/
#ifdef PI
	PRECISION b = 1.0/3.0 - cs2;
	PRECISION b2 = b*b;

   PRECISION zetabar = bulkViscosityToEntropyDensity(T);
	PRECISION tauPiInv = 15*b2*T/zetabar;

	double zeta_z = (Rhat-Rbar0)*(e-3*p)/3;
	double zeta_T = -(Rbar0+Rhat)*(e-3*p)/6;

//	PRECISION dPi = -((1-Rbar0)*(e-3*p)/3+Pi)*tauPiInv - zeta_z*zDzu -zeta_T*thetaT;

	PRECISION dPi;

	#ifndef PIMUNU
	dPi = -((1-Rbar0)*(e-3*p)/3+Pi)*tauPiInv - zeta_z*zDzu -zeta_T*thetaT - beta_PiPi * Pi * zDzu - delta_PiPi * Pi * thetaT;
	#else
	dPi = -((1-Rbar0)*(e-3*p)/3+Pi)*tauPiInv - zeta_z*zDzu -zeta_T*thetaT - beta_PiPi * Pi * zDzu - delta_PiPi * Pi * thetaT + lambda_Pipi * ps;
	#endif
	S[NUMBER_CONSERVED_VARIABLES-1] = dPi / ut + Pi * dkvk;
#endif
}
