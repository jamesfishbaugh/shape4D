#include "gridoptimize.h"

#include <stdlib.h>			// For malloc and free
#include <fftw3.h>			// For FFT

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EXGridOptimize
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
GridOptimize::GridOptimize(int, double*, int ntaux, double*, int ncy, double*,  int*, double, double*, double*)
{
    int Ndim, Nout;
    Ndim = ntaux;
    Nout = ncy;

    //this->_tauy = (double*) malloc(Nout*Ndim*sizeof(double));

    this->_tauy = new double[Nout*Ndim];
}

//----------------------------------------------------------------
// EXGrid
//----------------------------------------------------------------
// Inputs: 
//   grid - EXGrid object to copy
//
// Outputs:
//----------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------
GridOptimize::GridOptimize(const GridOptimize&)
{
    printf("In EXGridOptimize copy constructor...\n");
}


//----------------------------------------------------------------
// EXGrid
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
GridOptimize::~GridOptimize()
{
    //free(this->_tauy);
    delete [] this->_tauy;
}


//----------------------------------------------------------------
// GridOptim
//----------------------------------------------------------------
// Inputs: 
//   ncx - the number of points in cx
//   cx - array of points
//   ntaux - the dimension of points in taux
//   taux - array of weights
//   ncy - the number of points in cy
//   cy - array of points
//   grid_long - array of grid dimensions
//   grid_pas - accuracy of the grid 
//   grid_origin - origin of the grid
//   fft3k_R - FFT of the grid positions
//   
// Outputs: 
//   return - array of kernel computations
//----------------------------------------------------------------
// Computes an approximation to the kernel computation
// taux*K(cx,cy) with the corresponding grid information
//----------------------------------------------------------------
double* GridOptimize::GridOptim(int ncx, double* cx, int ntaux, double* taux, int ncy, double* cy,  int* grid_long, double grid_pas, double* grid_origin, double* fft3k_R)
{
    int Ndim, N, i, k, nx, ny, nz, nnx, Nin, Nout, c0X, c0Y, c0Z;
    double grid_pas3, rho000, rho100, rho010, rho001, rho110, rho101, rho011, rho111, deltaX, deltaY, deltaZ;
    double *alpha;
    fftw_complex *aux;
    fftw_plan forward, backward;
    nx = grid_long[0];
    ny = grid_long[1];
    nz = grid_long[2];
    nnx = (int) (nx/2 + 1);
    N = nx*ny*nz;
    grid_pas3 = grid_pas*grid_pas*grid_pas;
    Nin = ncx;
    Ndim = ntaux;
    Nout = ncy;

    //fftw_mpi_init();

    // Allocation
    //alpha = (double*) malloc(nx*ny*nz*Ndim*sizeof(double));
    alpha = new double[nx*ny*nz*Ndim];
    //this->_tauy = (double*) malloc(Nout*Ndim*sizeof(double));
    aux = (fftw_complex *) fftw_malloc(nnx*ny*nz*sizeof(fftw_complex));

    for (i=0; i<Nout*Ndim; i++)
    {
        this->_tauy[i] = 0.0;
    }

    // Tri-linear projection of the normals into the grid
    for (k=0; k<(nx*ny*nz*Ndim); k++)
    {
        alpha[k] = 0;
    }

    //printf("   In GridOptim...\n");

    for (k=0; k<Nin; k++)
    {
        c0X = (int) ((cx[3*k]   - grid_origin[0])/grid_pas);
        c0Y = (int) ((cx[1+3*k] - grid_origin[1])/grid_pas);
        c0Z = (int) ((cx[2+3*k] - grid_origin[2])/grid_pas);
        deltaX = cx[3*k]   - (grid_origin[0] + c0X*grid_pas);
        deltaY = cx[1+3*k] - (grid_origin[1] + c0Y*grid_pas);
        deltaZ = cx[2+3*k] - (grid_origin[2] + c0Z*grid_pas);
        rho000 = (grid_pas-deltaX) * (grid_pas-deltaY)  * (grid_pas-deltaZ)/grid_pas3;
        rho100 = deltaX            * (grid_pas-deltaY)  * (grid_pas-deltaZ)/grid_pas3;
        rho010 = (grid_pas-deltaX) * deltaY             * (grid_pas-deltaZ)/grid_pas3;
        rho001 = (grid_pas-deltaX) * (grid_pas-deltaY)  * deltaZ           /grid_pas3;
        rho110 = deltaX            * deltaY             * (grid_pas-deltaZ)/grid_pas3;
        rho011 = (grid_pas-deltaX) * deltaY             * deltaZ           /grid_pas3;
        rho101 = deltaX            * (grid_pas-deltaY)  * deltaZ           /grid_pas3;
        rho111 = deltaX            * deltaY             * deltaZ           /grid_pas3;

        for (i=0; i<Ndim; i++)
        {
            if (((c0X+1+(c0Y+1)*nx+(c0Z+1)*nx*ny+i*N)>=(nx*ny*nz*Ndim))||(c0X+c0Y*nx+c0Z*nx*ny+ i*N<0))
            {
                printf("ERREUR PROJECTION\n");
                printf("GRILLE: %f %d %d %d %f %f %f\n",grid_pas,nx,ny,nz,grid_origin[0],grid_origin[1],grid_origin[2]);
                printf("c0: %d %d %d\n",c0X,c0Y,c0Z);
                printf("%d >= %d\n",c0X+1+(c0Y+1)*nx+(c0Z+1)*nx*ny+i*N,nx*ny*nz*Ndim);
                return 0;
            }

            alpha[c0X   + c0Y*nx     + c0Z*nx*ny     + i*N] += rho000*taux[Ndim*k + i];
            alpha[c0X+1 + c0Y*nx     + c0Z*nx*ny     + i*N] += rho100*taux[Ndim*k + i];
            alpha[c0X   + (c0Y+1)*nx + c0Z*nx*ny     + i*N] += rho010*taux[Ndim*k + i];
            alpha[c0X   + c0Y*nx     + (c0Z+1)*nx*ny + i*N] += rho001*taux[Ndim*k + i];
            alpha[c0X+1 + (c0Y+1)*nx + c0Z*nx*ny     + i*N] += rho110*taux[Ndim*k + i];
            alpha[c0X   + (c0Y+1)*nx + (c0Z+1)*nx*ny + i*N] += rho011*taux[Ndim*k + i];
            alpha[c0X+1 + c0Y*nx     + (c0Z+1)*nx*ny + i*N] += rho101*taux[Ndim*k + i];
            alpha[c0X+1 + (c0Y+1)*nx + (c0Z+1)*nx*ny + i*N] += rho111*taux[Ndim*k + i];
        }
    }

    //printf("Beginning convolution...\n");
    // Convolution with Fast Fourier Transform
    for (i=0; i<Ndim; i++)
    {
        forward = fftw_plan_dft_r2c_3d(nz,ny,nx,alpha + i*N, aux, FFTW_ESTIMATE);
        //printf("executing forward fftw...\n");
        fftw_execute(forward);
        // ATTENTION : fft3k should be a real matrix
        for (k=0; k<(nnx*ny*nz); k++)
        {
            aux[k][0] *= fft3k_R[k]/N;
            aux[k][1] *= fft3k_R[k]/N;
        }
        backward = fftw_plan_dft_c2r_3d(nz,ny,nx, aux, alpha+i*N, FFTW_ESTIMATE);
        //printf("executing backwards fftw...\n");
        fftw_execute(backward);

        fftw_destroy_plan(forward);
        fftw_destroy_plan(backward);
    }
    //printf("Beginning interpolation on points...\n");
    // Interpolation on the points cy
    for (k=0; k<Nout; k++)
    {
        c0X = (int) ((cy[3*k]   - grid_origin[0])/grid_pas);
        c0Y = (int) ((cy[1+3*k] - grid_origin[1])/grid_pas);
        c0Z = (int) ((cy[2+3*k] - grid_origin[2])/grid_pas);
        deltaX = cy[3*k]   - (grid_origin[0] + c0X*grid_pas);
        deltaY = cy[1+3*k] - (grid_origin[1] + c0Y*grid_pas);
        deltaZ = cy[2+3*k] - (grid_origin[2] + c0Z*grid_pas);
        rho000 = (grid_pas-deltaX) * (grid_pas-deltaY)  * (grid_pas-deltaZ);
        rho100 = deltaX            * (grid_pas-deltaY)  * (grid_pas-deltaZ);
        rho010 = (grid_pas-deltaX) * deltaY             * (grid_pas-deltaZ);
        rho001 = (grid_pas-deltaX) * (grid_pas-deltaY)  * deltaZ;
        rho110 = deltaX            * deltaY             * (grid_pas-deltaZ);
        rho011 = (grid_pas-deltaX) * deltaY             * deltaZ;
        rho101 = deltaX            * (grid_pas-deltaY)  * deltaZ;
        rho111 = deltaX            * deltaY             * deltaZ;

        for (i=0; i<Ndim; i++)
        {
            if (((c0X+1+(c0Y+1)*nx+(c0Z+1)*nx*ny+i*N)>=(nx*ny*nz*Ndim))||(c0X+c0Y*nx+c0Z*nx*ny+ i*N<0))
            {
                printf("ERREUR INTERPOLATION\n");
                //return 0;
            }

            this->_tauy[i + k*Ndim] = (alpha[c0X   + c0Y*nx     + c0Z*nx*ny     + i*N]*rho000+
                    alpha[c0X+1 + c0Y*nx     + c0Z*nx*ny     + i*N]*rho100+
                    alpha[c0X   + (c0Y+1)*nx + c0Z*nx*ny     + i*N]*rho010+
                    alpha[c0X   + c0Y*nx     + (c0Z+1)*nx*ny + i*N]*rho001+
                    alpha[c0X+1 + (c0Y+1)*nx + c0Z*nx*ny     + i*N]*rho110+
                    alpha[c0X   + (c0Y+1)*nx + (c0Z+1)*nx*ny + i*N]*rho011+
                    alpha[c0X+1 + c0Y*nx     + (c0Z+1)*nx*ny + i*N]*rho101+
                    alpha[c0X+1 + (c0Y+1)*nx + (c0Z+1)*nx*ny + i*N]*rho111)/grid_pas3;
        }
    }

    //printf("Beginning memory cleanup...\n");

    // Clean up memory
    fftw_free(aux);
    delete [] alpha;

    fftw_cleanup();

    return this->_tauy;
}

