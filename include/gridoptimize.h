//--------------------------------------------------------------------------
// Erid: Static class for fast kernel computations
//--------------------------------------------------------------------------

#ifndef GRIDOPTIMIZE_H
#define GRIDOPTIMIZE_H

class GridOptimize
{

    private:

        double *_tauy;

    public:

        GridOptimize(int ncx, double* cx, int ntaux, double* taux, int ncy, double* cy,  int* grid_long, double grid_pas, double* grid_origin, double* fft3k_R);
        GridOptimize(const GridOptimize& gridOptim);
        ~GridOptimize();

        double* GridOptim(int ncx, double* cx, int ntaux, double* taux, int ncy, double* cy,  int* grid_long, double grid_pas, double* grid_origin, double* fft3k_R);

};

#endif // GRIDOPTIMIZE_H
