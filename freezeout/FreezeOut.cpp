#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

//return a 4 dimensional linear interpolation inside the hypercube, given the values
//on the corners (a0000 through a1111) and edge lengths x0 through x3
double linearInterp4D(double x0, double x1, double x2, double x3,
                      double a0000, double a1000, double a0100, double a0010, double a0001,
                      double a1100, double a1010, double a1001,
                      double a0110, double a0101, double a0011,
                      double a1110, double a1101, double a0111, double a1011, double a1111)
{
  double result = 0;
  result = ((1-x0) * (1-x1) * (1-x2) * (1-x3) * a0000)
            + ((x0) * (1-x1) * (1-x2) * (1-x3) * a1000)
            + ((1-x0) * (x1) * (1-x2) * (1-x3) * a0100)
            + ((1-x0) * (1-x1) * (x2) * (1-x3) * a0010)
            + ((1-x0) * (1-x1) * (1-x2) * (x3) * a0001)
            + ((x0) * (x1) * (1-x2) * (1-x3) * a1100)
            + ((x0) * (1-x1) * (x2) * (1-x3) * a1010)
            + ((x0) * (1-x1) * (1-x2) * (x3) * a1001)
            + ((1-x0) * (x1) * (x2) * (1-x3) * a0110)
            + ((1-x0) * (x1) * (1-x2) * (x3) * a0101)
            + ((1-x0) * (1-x1) * (x2) * (x3) * a0011)
            + ((x0) * (x1) * (x2) * (1-x3) * a1110)
            + ((x0) * (x1) * (1-x2) * (x3) * a1101)
            + ((x0) * (1-x1) * (x2) * (x3) * a1011)
            + ((1-x0) * (x1) * (x2) * (x3) * a0111)
            + ((x0) * (x1) * (x2) * (x3) * a1111);

  return result;
}

double linearInterp3D(double x0, double x1, double x2,
                      double a000, double a100, double a010, double a001,
                      double a110, double a101, double a011, double a111)
{
  double result = 0;
  result = ((1-x0) * (1-x1) * (1-x2) * a000)
            + ((x0) * (1-x1) * (1-x2) * a100)
            + ((1-x0) * (x1) * (1-x2) * a010)
            + ((1-x0) * (1-x1) * (x2) * a001)
            + ((x0) * (x1) * (1-x2) * a110)
            + ((x0) * (1-x1) * (x2) * a101)
            + ((1-x0) * (x1) * (x2) * a011)
            + ((x0) * (x1) * (x2)  * a111);

  return result;
}

void swapAndSetHydroVariables(double ****energy_density_evoution, double *****hydrodynamic_evolution,
                              CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e,
                              FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz, int FOFREQ)
{
#pragma omp parallel for collapse(3)
  for (int ix = 2; ix < nx+2; ix++)
  {
    for (int iy = 2; iy < ny+2; iy++)
    {
      for (int iz = 2; iz < nz+2; iz++)
      {
        int s = columnMajorLinearIndex(ix, iy, iz, nx+4, ny+4);
        //previous hydro variable values written to zeroth index
        energy_density_evoution[0][ix-2][iy-2][iz-2] = energy_density_evoution[FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[0][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[0][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[1][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[1][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[2][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[2][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[3][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[3][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[4][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[4][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[5][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[5][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[6][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[6][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[7][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[7][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[8][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[8][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[9][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[9][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[10][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[10][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[11][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[11][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[12][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[12][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[13][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[13][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[14][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[14][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[15][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[15][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[16][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[16][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[17][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[17][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[18][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[18][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[19][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[19][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evolution[20][0][ix-2][iy-2][iz-2] = hydrodynamic_evolution[20][FOFREQ][ix-2][iy-2][iz-2];

        //current hydro variable values written to first index
        energy_density_evoution[1][ix-2][iy-2][iz-2] = (double)e[s];
        hydrodynamic_evolution[0][1][ix-2][iy-2][iz-2] = (double)(u->ut[s]);
        hydrodynamic_evolution[1][1][ix-2][iy-2][iz-2] = (double)(u->ux[s]);
        hydrodynamic_evolution[2][1][ix-2][iy-2][iz-2] = (double)(u->uy[s]);
        hydrodynamic_evolution[3][1][ix-2][iy-2][iz-2] = (double)(u->un[s]);
        hydrodynamic_evolution[4][1][ix-2][iy-2][iz-2] = (double)(e[s]);
        hydrodynamic_evolution[5][1][ix-2][iy-2][iz-2] = (double)(q->pl[s]);
	#ifdef PIMUNU
        hydrodynamic_evolution[6][1][ix-2][iy-2][iz-2] = (double)(q->pitt[s]);
        hydrodynamic_evolution[7][1][ix-2][iy-2][iz-2] = (double)(q->pitx[s]);
        hydrodynamic_evolution[8][1][ix-2][iy-2][iz-2] = (double)(q->pity[s]);
        hydrodynamic_evolution[9][1][ix-2][iy-2][iz-2] = (double)(q->pitn[s]);
        hydrodynamic_evolution[10][1][ix-2][iy-2][iz-2] = (double)(q->pixx[s]);
        hydrodynamic_evolution[11][1][ix-2][iy-2][iz-2] = (double)(q->pixy[s]);
        hydrodynamic_evolution[12][1][ix-2][iy-2][iz-2] = (double)(q->pixn[s]);
        hydrodynamic_evolution[13][1][ix-2][iy-2][iz-2] = (double)(q->piyy[s]);
        hydrodynamic_evolution[14][1][ix-2][iy-2][iz-2] = (double)(q->piyn[s]);
        hydrodynamic_evolution[15][1][ix-2][iy-2][iz-2] = (double)(q->pinn[s]);
	#endif
  #ifdef W_TZ_MU
        hydrodynamic_evolution[16][1][ix-2][iy-2][iz-2] = (double)(q->WtTz[s]);
        hydrodynamic_evolution[17][1][ix-2][iy-2][iz-2] = (double)(q->WxTz[s]);
        hydrodynamic_evolution[18][1][ix-2][iy-2][iz-2] = (double)(q->WyTz[s]);
        hydrodynamic_evolution[19][1][ix-2][iy-2][iz-2] = (double)(q->WnTz[s]);
	#endif
	#ifdef PI
        hydrodynamic_evolution[20][1][ix-2][iy-2][iz-2] = (double)(q->Pi[s]);
	#endif
      }
    }
  }
}

void setHydroVariables(double ****energy_density_evoution, double *****hydrodynamic_evolution,
                              CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e,
                              FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz, int FOFREQ, int n)
{
  int nFO = n % FOFREQ;
#pragma omp parallel for collapse(3)
  for (int ix = 2; ix < nx+2; ix++)
  {
    for (int iy = 2; iy < ny+2; iy++)
    {
      for (int iz = 2; iz < nz+2; iz++)
      {
        int s = columnMajorLinearIndex(ix, iy, iz, nx+4, ny+4);
        energy_density_evoution[nFO+1][ix-2][iy-2][iz-2] = (double)e[s];
        hydrodynamic_evolution[0][nFO+1][ix-2][iy-2][iz-2] = (double)(u->ut[s]);
        hydrodynamic_evolution[1][nFO+1][ix-2][iy-2][iz-2] = (double)(u->ux[s]);
        hydrodynamic_evolution[2][nFO+1][ix-2][iy-2][iz-2] = (double)(u->uy[s]);
        hydrodynamic_evolution[3][nFO+1][ix-2][iy-2][iz-2] = (double)(u->un[s]);
        hydrodynamic_evolution[4][nFO+1][ix-2][iy-2][iz-2] = (double)(e[s]);
        hydrodynamic_evolution[5][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pl[s]);
	#ifdef PIMUNU
        hydrodynamic_evolution[6][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pitt[s]);
        hydrodynamic_evolution[7][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pitx[s]);
        hydrodynamic_evolution[8][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pity[s]);
        hydrodynamic_evolution[9][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pitn[s]);
        hydrodynamic_evolution[10][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pixx[s]);
        hydrodynamic_evolution[11][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pixy[s]);
        hydrodynamic_evolution[12][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pixn[s]);
        hydrodynamic_evolution[13][nFO+1][ix-2][iy-2][iz-2] = (double)(q->piyy[s]);
        hydrodynamic_evolution[14][nFO+1][ix-2][iy-2][iz-2] = (double)(q->piyn[s]);
        hydrodynamic_evolution[15][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pinn[s]);
	#endif
  #ifdef W_TZ_MU
        hydrodynamic_evolution[16][nFO+1][ix-2][iy-2][iz-2] = (double)(q->WtTz[s]);
        hydrodynamic_evolution[17][nFO+1][ix-2][iy-2][iz-2] = (double)(q->WxTz[s]);
        hydrodynamic_evolution[18][nFO+1][ix-2][iy-2][iz-2] = (double)(q->WyTz[s]);
        hydrodynamic_evolution[19][nFO+1][ix-2][iy-2][iz-2] = (double)(q->WnTz[s]);
  #endif
	#ifdef PI
        hydrodynamic_evolution[20][nFO+1][ix-2][iy-2][iz-2] = (double)(q->Pi[s]);
	#endif
      }
    }
  }
}
void writeEnergyDensityToHypercube4D(double ****hyperCube, double ****energy_density_evoution, int it, int ix, int iy, int iz)
{
  hyperCube[0][0][0][0] = energy_density_evoution[it][ix][iy][iz];
  hyperCube[1][0][0][0] = energy_density_evoution[it+1][ix][iy][iz];
  hyperCube[0][1][0][0] = energy_density_evoution[it][ix+1][iy][iz];
  hyperCube[0][0][1][0] = energy_density_evoution[it][ix][iy+1][iz];
  hyperCube[0][0][0][1] = energy_density_evoution[it][ix][iy][iz+1];
  hyperCube[1][1][0][0] = energy_density_evoution[it+1][ix+1][iy][iz];
  hyperCube[1][0][1][0] = energy_density_evoution[it+1][ix][iy+1][iz];
  hyperCube[1][0][0][1] = energy_density_evoution[it+1][ix][iy][iz+1];
  hyperCube[0][1][1][0] = energy_density_evoution[it][ix+1][iy+1][iz];
  hyperCube[0][1][0][1] = energy_density_evoution[it][ix+1][iy][iz+1];
  hyperCube[0][0][1][1] = energy_density_evoution[it][ix][iy+1][iz+1];
  hyperCube[1][1][1][0] = energy_density_evoution[it+1][ix+1][iy+1][iz];
  hyperCube[1][1][0][1] = energy_density_evoution[it+1][ix+1][iy][iz+1];
  hyperCube[1][0][1][1] = energy_density_evoution[it+1][ix][iy+1][iz+1];
  hyperCube[0][1][1][1] = energy_density_evoution[it][ix+1][iy+1][iz+1];
  hyperCube[1][1][1][1] = energy_density_evoution[it+1][ix+1][iy+1][iz+1];
}
void writeEnergyDensityToHypercube3D(double ***hyperCube, double ****energy_density_evoution, int it, int ix, int iy)
{
  hyperCube[0][0][0] = energy_density_evoution[it][ix][iy][0];
  hyperCube[1][0][0] = energy_density_evoution[it+1][ix][iy][0];
  hyperCube[0][1][0] = energy_density_evoution[it][ix+1][iy][0];
  hyperCube[0][0][1] = energy_density_evoution[it][ix][iy+1][0];
  hyperCube[1][1][0] = energy_density_evoution[it+1][ix+1][iy][0];
  hyperCube[1][0][1] = energy_density_evoution[it+1][ix][iy+1][0];
  hyperCube[0][1][1] = energy_density_evoution[it][ix+1][iy+1][0];
  hyperCube[1][1][1] = energy_density_evoution[it+1][ix+1][iy+1][0];
}
double interpolateVariable4D(double *****hydrodynamic_evolution, int ivar, int it, int ix, int iy, int iz, double tau_frac, double x_frac, double y_frac, double z_frac)
{
  double result = linearInterp4D(tau_frac, x_frac, y_frac, z_frac,
    hydrodynamic_evolution[ivar][it][ix][iy][iz], hydrodynamic_evolution[ivar][it+1][ix][iy][iz], hydrodynamic_evolution[ivar][it][ix+1][iy][iz], hydrodynamic_evolution[ivar][it][ix][iy+1][iz], hydrodynamic_evolution[ivar][it][ix][iy][iz+1],
    hydrodynamic_evolution[ivar][it+1][ix+1][iy][iz], hydrodynamic_evolution[ivar][it+1][ix][iy+1][iz], hydrodynamic_evolution[ivar][it+1][ix][iy][iz+1],
    hydrodynamic_evolution[ivar][it][ix+1][iy+1][iz], hydrodynamic_evolution[ivar][it][ix+1][iy][iz+1], hydrodynamic_evolution[ivar][it][ix][iy+1][iz+1],
    hydrodynamic_evolution[ivar][it+1][ix+1][iy+1][iz], hydrodynamic_evolution[ivar][it+1][ix+1][iy][iz+1], hydrodynamic_evolution[ivar][it][ix+1][iy+1][iz+1], hydrodynamic_evolution[ivar][it+1][ix][iy+1][iz+1], hydrodynamic_evolution[ivar][it+1][ix+1][iy+1][iz+1]);
    return result;
}

double interpolateVariable3D(double *****hydrodynamic_evolution, int ivar, int it, int ix, int iy, double tau_frac, double x_frac, double y_frac)
{
  double result = linearInterp3D(tau_frac, x_frac, y_frac,
    hydrodynamic_evolution[ivar][it][ix][iy][0], hydrodynamic_evolution[ivar][it+1][ix][iy][0], hydrodynamic_evolution[ivar][it][ix+1][iy][0], hydrodynamic_evolution[ivar][it][ix][iy+1][0],
    hydrodynamic_evolution[ivar][it+1][ix+1][iy][0], hydrodynamic_evolution[ivar][it+1][ix][iy+1][0], hydrodynamic_evolution[ivar][it][ix+1][iy+1][0], hydrodynamic_evolution[ivar][it+1][ix+1][iy+1][0]);
    return result;
}
