/*
 * InitialConditions.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include <math.h> // for math functions
#include <stdio.h> // for printf
#include <stdlib.h> //TEMP

#include "edu/osu/rhic/trunk/ic/InitialConditions.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"
#include "edu/osu/rhic/core/ic/GlauberModel.h"
#include "edu/osu/rhic/core/ic/MonteCarloGlauberModel.h"
#include "edu/osu/rhic/harness/hydro/HydroParameters.h"
#include "edu/osu/rhic/trunk/eos/EquationOfState.h"

#define THETA_FUNCTION(X) ((double)X < (double)0 ? (double)0 : (double)1)

//*********************************************************************************************************\
//* Read in all initial profiles from a single or seperate file
//*********************************************************************************************************/

//this reads all hydro variables from a single file; this way we do not need to fetch the coordinates many times
//note that the file must contain values for all dissipative currents, even if they are zero !!!
void setInitialTmunuFromFile(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;

    float x, y, z, e_in, p_in, ut_in, ux_in, uy_in, un_in;
    //#ifdef PIMUNU
    float pitt_in, pitx_in, pity_in, pitn_in, pixx_in, pixy_in, pixn_in, piyy_in, piyn_in, pinn_in;
    //#endif
    //#ifdef PI
    float Pi_in;
    //#endif
    FILE *fileIn;
    char fname[255];

    sprintf(fname, "%s/%s", rootDirectory, "/input/Tmunu.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open Tmunu.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &x, &y, &z, &e_in, &p_in, &ut_in, &ux_in, &uy_in, &un_in, &pitt_in, &pitx_in, &pity_in, &pitn_in, &pixx_in, &pixy_in, &pixn_in, &piyy_in, &piyn_in, &pinn_in, &Pi_in);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    e[s] =  (PRECISION) e_in + (PRECISION)1.0e-3;
                    //ep[s] = (PRECISION) e_in; //set previous step to same value
                    p[s] = p_in;
                    u->ut[s] = ut_in;
                    u->ux[s] = ux_in;
                    u->uy[s] = uy_in;
                    u->un[s] = un_in;

		    //There is no host variable up unlike cpu-vh. This is accomplished in copyHostToDeviceMemory()
		    /*
                    up->ut[s] = ut_in; //set previous step to same value
                    up->ux[s] = ux_in; //...
                    up->uy[s] = uy_in;
                    up->un[s] = un_in;
		    */
#ifdef PIMUNU
                    q->pitt[s] = pitt_in;
                    q->pitx[s] = pitx_in;
                    q->pity[s] = pity_in;
                    q->pitn[s] = pitn_in;
                    q->pixx[s] = pixx_in;
                    q->pixy[s] = pixy_in;
                    q->pixn[s] = pixn_in;
                    q->piyy[s] = piyy_in;
                    q->piyn[s] = piyn_in;
                    q->pinn[s] = pinn_in;
#endif
#ifdef PI
                    q->Pi[s] = Pi_in;
#endif
                }
            }
        }
    }
    fclose(fileIn);
}

//this function reads a separate file for every hydrodynamic variable
void setInitialTmunuFromFiles(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;

    float x, y, z, value;
    FILE *fileIn;
    char fname[255];

    //energy density
    sprintf(fname, "%s/%s", rootDirectory, "/input/e.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open e.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    e[s] =  (PRECISION) value + (PRECISION)1.0e-3;

                    //set the pressure based on EOS rather than reading from file
                    p[s] = equilibriumPressure(e[s]);
                    //ep[s] = (PRECISION) value;
                    //printf("e [ %d ] = %f\n", s, e[s]);
                }
            }
        }
    }
    fclose(fileIn);

		//longitidunal pressure
    sprintf(fname, "%s/%s", rootDirectory, "/input/pL.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pL.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pl[s] =  (PRECISION) value + (PRECISION)1.0e-3;
                }
            }
        }
    }
    fclose(fileIn);

    //pressure
    /*
    sprintf(fname, "%s/%s", rootDirectory, "/input/p.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open p.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    p[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    */
    //ut
    sprintf(fname, "%s/%s", rootDirectory, "/input/ut.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open ut.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    u->ut[s] =  (PRECISION) value;
                    up->ut[s] = (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //ux
    sprintf(fname, "%s/%s", rootDirectory, "/input/ux.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open ux.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    u->ux[s] =  (PRECISION) value;
                    up->ux[s] = (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //uy
    sprintf(fname, "%s/%s", rootDirectory, "/input/uy.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open uy.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    u->uy[s] =  (PRECISION) value;
                    up->uy[s] = (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //un
    sprintf(fname, "%s/%s", rootDirectory, "/input/un.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open un.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    u->un[s] =  (PRECISION) value;
                    up->un[s] = (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

#ifdef PIMUNU
    //pitt
    sprintf(fname, "%s/%s", rootDirectory, "/input/residual_pitt.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open residual_pitt.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pitt[s] =  (PRECISION) value;
                    //q->pitt[s] =  (PRECISION) 0.0;
                }
            }
        }
    }
    fclose(fileIn);

    //pitx
    sprintf(fname, "%s/%s", rootDirectory, "/input/residual_pitx.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open residual_pitx.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pitx[s] =  (PRECISION) value;
                    //q->pitx[s] =  (PRECISION) 0.0;
                }
            }
        }
    }
    fclose(fileIn);

    //pity
    sprintf(fname, "%s/%s", rootDirectory, "/input/residual_pity.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open residual_pity.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pity[s] =  (PRECISION) value;
                    //q->pity[s] =  (PRECISION) 0.0;
                }
            }
        }
    }
    fclose(fileIn);

    //pitn
    sprintf(fname, "%s/%s", rootDirectory, "/input/residual_pitn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open residual_pitn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pitn[s] =  (PRECISION) value;
                    //q->pitn[s] =  (PRECISION) 0.0;
                }
            }
        }
    }
    fclose(fileIn);

    //pixx
    sprintf(fname, "%s/%s", rootDirectory, "/input/residual_pixx.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open residual_pixx.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pixx[s] =  (PRECISION) value;
                    //q->pixx[s] =  (PRECISION) 0.0;
                }
            }
        }
    }
    fclose(fileIn);

    //pixy
    sprintf(fname, "%s/%s", rootDirectory, "/input/residual_pixy.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open residual_pixy.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pixy[s] =  (PRECISION) value;
                    //q->pixy[s] =  (PRECISION) 0.0;
                }
            }
        }
    }
    fclose(fileIn);

    //pixn
    sprintf(fname, "%s/%s", rootDirectory, "/input/residual_pixn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open residual_pixn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pixn[s] =  (PRECISION) value;
                    //q->pitn[s] =  (PRECISION) 0.0;
                }
            }
        }
    }
    fclose(fileIn);

    //piyy
    sprintf(fname, "%s/%s", rootDirectory, "/input/residual_piyy.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open residual_piyy.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->piyy[s] =  (PRECISION) value;
                    //q->piyy[s] =  (PRECISION) 0.0;
                }
            }
        }
    }
    fclose(fileIn);

    //piyn
    sprintf(fname, "%s/%s", rootDirectory, "/input/residual_piyn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open residual_piyn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->piyn[s] =  (PRECISION) value;
                    //q->piyn[s] =  (PRECISION) 0.0;
                }
            }
        }
    }
    fclose(fileIn);

    //pinn
    sprintf(fname, "%s/%s", rootDirectory, "/input/residual_pinn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open residual_pinn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pinn[s] =  (PRECISION) value;
                    //q->pinn[s] =  (PRECISION) 0.0;
                }
            }
        }
    }
    fclose(fileIn);
#endif
#ifdef PI
    //bulk
    sprintf(fname, "%s/%s", rootDirectory, "/input/residualBulk.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open residualBulk.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->Pi[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
#endif
#ifdef W_TZ_MU
    //longitudinal momentum diffusion vector
    sprintf(fname, "%s/%s", rootDirectory, "/input/W_perp_t.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open W_perp_t.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->WtTz[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
		sprintf(fname, "%s/%s", rootDirectory, "/input/W_perp_x.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open W_perp_x.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->WxTz[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
		sprintf(fname, "%s/%s", rootDirectory, "/input/W_perp_y.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open W_perp_y.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->WyTz[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
		sprintf(fname, "%s/%s", rootDirectory, "/input/W_perp_n.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open W_perp_n.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->WnTz[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
#endif
}

/*********************************************************************************************************\
 * Set initial flow profile
 *		- u^\mu = (1, 0, 0, 0)
 * 	- No transverse flow (ux = uy = 0)
 *		- Longitudinal scaling flow (u_z = z/t, i.e. un = 0)
/*********************************************************************************************************/
void setFluidVelocityInitialCondition(void * latticeParams, void * hydroParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double t0 = hydro->initialProperTimePoint;

	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				PRECISION ux = 0;
				PRECISION uy = 0;
				PRECISION un = 0;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = sqrt(1+ux*ux+uy*uy+t0*t0*un*un);
			}
		}
	}
}

/*********************************************************************************************************\
 * Set initial shear-stress tensor \pi^\mu\nu
 *		- Navier-Stokes value, i.e. \pi^\mu\nu = 2 * (\epsilon + P) / T * \eta/S * \sigma^\mu\nu
 * 	- No initial pressure anisotropies (\pi^\mu\nu = 0)
/*********************************************************************************************************\
void setPimunuNavierStokesInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	PRECISION dx = (PRECISION)(lattice->latticeSpacingX);
	PRECISION dz = (PRECISION)(lattice->latticeSpacingRapidity);

	PRECISION etabar = (PRECISION)(hydro->shearViscosityToEntropyDensity);
	PRECISION t = hydro->initialProperTimePoint;

	PRECISION e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
//				double T = pow(e[s]/e0, 0.25);
				PRECISION T = effectiveTemperature(e[s]);
				if (T == 0) T = 1.e-3;
				PRECISION pinn = -2/(3*t*t*t)*etabar*(e[s]+p[s])/T;
#ifdef PIMUNU
				q->pitt[s] = 0;
				q->pitx[s] = 0;
				q->pity[s] = 0;
				q->pitn[s] = 0;
				q->pixx[s] = -t*t*pinn/2;
				q->pixy[s] = 0;
				q->pixn[s] = 0;
				q->piyy[s] = -t*t*pinn/2;
				q->piyn[s] = 0;
				q->pinn[s] = pinn;
#endif
#ifdef PI
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
				PRECISION x = T/1.01355;
				PRECISION zetabar = A_1*x*x + A_2*x - A_3;
				if(x > 1.05)
					zetabar = LAMBDA_1*exp(-(x-1)/SIGMA_1) + LAMBDA_2*exp(-(x-1)/SIGMA_2)+0.001;
				else if(x < 0.995)
					zetabar = LAMBDA_3*exp((x-1)/SIGMA_3)+ LAMBDA_4*exp((x-1)/SIGMA_4)+0.03;
				q->Pi[s] = -zetabar*(e[s]+p[s])/T/t;
#endif
			}
		}
	}
}
/*********************************************************************************************************/
void setPimunuInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
	int initializePimunuNavierStokes = hydro->initializePimunuNavierStokes;
	if (initializePimunuNavierStokes==1) {
		printf("Initialize \\pi^\\mu\\nu to its asymptotic Navier-Stokes value.\n");
#ifdef PI
		printf("Initialize \\Pi to its asymptotic Navier-Stokes value.\n");
#endif
//		setPimunuNavierStokesInitialCondition(latticeParams, initCondParams, hydroParams);
		return;
	}
	else {
		printf("Initialize \\pi^\\mu\\nu to zero.\n");
		struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
		int nx = lattice->numLatticePointsX;
		int ny = lattice->numLatticePointsY;
		int nz = lattice->numLatticePointsRapidity;
		for(int i = 2; i < nx+2; ++i) {
			for(int j = 2; j < ny+2; ++j) {
				for(int k = 2; k < nz+2; ++k) {
					int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
#ifdef PIMUNU
			  		q->pitt[s] = 0;
			  		q->pitx[s] = 0;
			  		q->pity[s] = 0;
			  		q->pitn[s] = 0;
			  		q->pixx[s] = 0;
			  		q->pixy[s] = 0;
			  		q->pixn[s] = 0;
			  		q->piyy[s] = 0;
			  		q->piyn[s] = 0;
			  		q->pinn[s] = 0;
#endif
#ifdef PI
			  		q->Pi[s] = 0;
#endif
				}
			}
		}
		return;
	}
}

/*********************************************************************************************************\
 * Constant initial energy density distribution
/*********************************************************************************************************/
void setConstantEnergyDensityInitialCondition(void * latticeParams, void * initCondParams) {
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	double initialEnergyDensity = initCond->initialEnergyDensity;

	double T0 = 3.05;
	double ed = equilibriumEnergyDensity(T0);

	// temp
//	ed = 93.2104;

	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				e[s] = (PRECISION) ed;
				p[s] = equilibriumPressure(e[s]);
			}
		}
	}
}

/*********************************************************************************************************\
 * Longitudinal initial energy density distribution
/*********************************************************************************************************/
void longitudinalEnergyDensityDistribution(double * const __restrict__ eL, void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nz = lattice->numLatticePointsRapidity;

	double dz = lattice->latticeSpacingRapidity;

	double etaFlat = initCond->rapidityMean;
	double etaVariance = initCond->rapidityVariance;
	etaFlat=5.9;
	etaVariance=0.16;

	for(int k = 0; k < nz; ++k) {
		double eta = (k - (nz-1)/2)*dz;
		double etaScaled = fabs(eta) - etaFlat/2;
		double arg = -etaScaled * etaScaled / etaVariance / 2 * THETA_FUNCTION(etaScaled);
		eL[k] = exp(arg);
	}
}

/*********************************************************************************************************\
 * Continuous optical glauber Glauber initial energy density distribution
/*********************************************************************************************************/
void setGlauberInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;
	double T0 = 3.05;
//	e0 *= pow(T0,4);
	e0 = (double) equilibriumEnergyDensity(T0);

	double eT[nx*ny], eL[nz];
	energyDensityTransverseProfileAA(eT, nx, ny, dx, dy, initCondParams);
	longitudinalEnergyDensityDistribution(eL, latticeParams, initCondParams);

	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			double energyDensityTransverse = e0 * eT[i-2+(j-2)*nx];
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				double energyDensityLongitudinal = eL[k-2];
				double ed = (energyDensityTransverse * energyDensityLongitudinal) + 1.e-3;
				e[s] = (PRECISION) ed;
				p[s] = equilibriumPressure(e[s]);
			}
		}
	}
}

/*********************************************************************************************************\
 * Monte carlo Glauber initial energy density distribution
/*********************************************************************************************************/
void setMCGlauberInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;
//	double T0 = 3.05;
//	double T0 = 2.03;
//	e0 *= pow(T0,4);
//	e0 = (double) equilibriumEnergyDensity(T0);

	double eT[nx*ny], eL[nz];
	monteCarloGlauberEnergyDensityTransverseProfile(eT, nx, ny, dx, dy, initCondParams);
	longitudinalEnergyDensityDistribution(eL, latticeParams, initCondParams);

	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			double energyDensityTransverse = e0 * eT[i-2 + nx*(j-2)];
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				double energyDensityLongitudinal = eL[k-2];
				double ed = (energyDensityTransverse * energyDensityLongitudinal) + 1.e-3;
				e[s] = (PRECISION) ed;
				p[s] = equilibriumPressure(e[s]);
			}
		}
	}
}

void setMCGlauberInitialCondition_FromMS_1d(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;
//	double T0 = 3.05;
//	double T0 = 2.03;
//	e0 *= pow(T0,4);
//	e0 = (double) equilibriumEnergyDensity(T0);

double T [201] = {0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,
   0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,
   0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,
   0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,
   0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,
   0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,
   0.000599999983433382,0.0005999998354381784,0.000599999786106444,0.000599999884769913,0.00060000018076032,
   0.0005999941011173602,0.0005999899730378173,0.0005999895014264355,0.0005999943911879592,0.0006000063472271326,
   0.000599877946588576,0.0005997973040021889,0.0005998034063377464,0.0005999352404650242,0.0006002317932537975,
   0.000598531998290537,0.0005976249090491488,0.0005980995257202346,0.0006005448484943958,0.0006055498775622337,
   0.0005929250312914075,0.0005892325371511969,0.0006002560407879391,0.0006317791878479711,0.00068958562397763,
   0.0006586336190203106,0.0006957385383760284,0.0008368903716418574,0.0011180791084148677,0.0015752947382921326,
   0.0018502728711497886,0.002471821471236077,0.003574494123078321,0.005292844411203807,0.00776142592013985,0.010610427580255245,
   0.014604858793775426,0.020005364308767412,0.02707258887329805,0.03606717723543432,0.04794849148960996,0.06210377970093336,
   0.07861900728088003,0.09758013964092493,0.11907314219254352,0.14592531989125102,0.17479596371847328,0.20508570419967448,
   0.23619517186032044,0.26752499722587636,0.2996607604657861,0.3305219050505423,0.3592128240946144,0.3848379107124737,
   0.4065015580185907,0.4175742803214973,0.4243278192430877,0.42730003759931706,0.42702879820614076,0.42405196387951427,
   0.41619954971895473,0.40739422818596527,0.39885082402561056,0.39178416198295535,0.38740906680306453,0.3980978976384968,
   0.4111185612249494,0.4248964987056141,0.43785715122368213,0.448425959922345,0.4492893254608111,0.44604648958725057,
   0.43855765356585064,0.42668301866079844,0.4102827861362813,0.3847054548080138,0.3554508540007741,0.3235071105908668,
   0.2898623514545971,0.25550470346827014,0.22331697029189684,0.19191893282215056,0.1618250487394093,0.13354977572405166,
   0.10760757145645614,0.08717845484153713,0.06944493202900376,0.054255070393100166,0.041456937308070894,0.03089860014816049,
   0.022867599188902847,0.016662660677930922,0.012021983764166723,0.008683767596532416,0.006386211323950163,0.004303866307926745,
   0.0028794914316535936,0.0019921977909066707,0.0015210964814619799,0.0013452985990955251,0.0009800641972685757,
   0.0007593181746505598,0.0006531343875961552,0.0006315866924600479,0.0006647489455969237,0.0006209835593325001,
   0.0005975036950576741,0.0005898110701343721,0.000593407401924522,0.0006037944077900512,0.0006002026085996673,
   0.000598472717331823,0.0005981742504717511,0.0005988767245046843,0.000600149655915855,0.0005999480834887839,
   0.0005998596208358438,0.0005998574038676954,0.0005999145684949993,0.0006000042506284163,0.0005999965223188895,
   0.0005999933493017264,0.0005999936334525172,0.0005999962766468517,0.00060000018076032,0.0005999999854066514,
   0.0005999999202887618,0.0005999999528477066,0.0006000000505245409,0.00060000018076032,0.00060000018076032,0.00060000018076032,
   0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,
   0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,
   0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,
   0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,
   0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,0.00060000018076032,
   0.00060000018076032,0.00060000018076032,0.00060000018076032};

	int ind = 0;
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				double ed = (double) equilibriumEnergyDensity(T[ind]/0.197326938);
++ind;
				e[s] = (PRECISION) ed;
				p[s] = equilibriumPressure(e[s]);
			}
		}
	}
}

/*********************************************************************************************************\
 * Initial conditions for the Gubser ideal hydro test
 *		- set energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny
/*********************************************************************************************************/
void setIdealGubserInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			double T = 1.9048812623618392/pow(1 + pow(1 - pow(x,2) - pow(y,2),2) + 2*(1 + pow(x,2) + pow(y,2)),0.3333333333333333);
			double r = sqrt(x*x+y*y);
			double phi = atanh(2*1*r/(1+1+x*x+y*y));

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);

				e[s] = (PRECISION) (e0 * pow(T,4));
				p[s] = e[s]/3;
				u->ux[s] = (PRECISION) (sinh(phi)*x/r);
				u->uy[s] = (PRECISION) (sinh(phi)*y/r);
				u->un[s] = 0;
				u->ut[s] = sqrt(1 + u->ux[s]*u->ux[s] + u->uy[s]*u->uy[s]);
#ifdef PIMUNU
				q->pixx[s] = 0;
				q->pixy[s] = 0;
				q->piyy[s] = 0;
#endif
			}
		}
	}
}

/*********************************************************************************************************\
 * Initial conditions to use.
 *	Set the energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny.
 * 	0 - constant energy density
 *		1 - Isreal-Stewart hydrodynamic Gubser flow test
 *		2 - Continous optical Glauber
 *		3 - Ideal hydrodynamic Gubser flow test
 *		4 - Monte carlo Glauber
 *		5 - Relativistic Sod shock-tube test
/*********************************************************************************************************/
void setInitialConditions(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	int initialConditionType = initCond->initialConditionType;
	printf("Setting initial conditions: ");
	switch (initialConditionType) {
		case 0: {
			printf("constant energy density.\n");
			setConstantEnergyDensityInitialCondition(latticeParams, initCondParams);
			setFluidVelocityInitialCondition(latticeParams, hydroParams);
			setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
			return;
		}
		case 2: {
			printf("Continous optical Glauber.\n");
			setGlauberInitialCondition(latticeParams, initCondParams);
			setFluidVelocityInitialCondition(latticeParams, hydroParams);
			setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
			return;
		}
		case 3: {
			printf("Ideal hydrodynamic Gubser flow test.\n");
			setIdealGubserInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 4: {
			printf("Monte carlo Glauber.\n");
			setMCGlauberInitialCondition(latticeParams, initCondParams);
			setFluidVelocityInitialCondition(latticeParams, hydroParams);
			setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
			return;
		}
		case 10: {
      printf("Reading initial T ^mu nu from input/e.dat , input/p.dat , etc... \n");
      setInitialTmunuFromFiles(latticeParams, initCondParams, hydroParams, rootDirectory);
      return;
    }
    case 11: {
      printf("Reading initial T ^mu nu from /input/Tmunu.dat \n");
      setInitialTmunuFromFile(latticeParams, initCondParams, hydroParams, rootDirectory);
      return;
    }
		default: {
			printf("Initial condition type not defined. Exiting ...\n");
			exit(-1);
		}
	}
}
