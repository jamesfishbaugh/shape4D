//----------------------------------------------------------------
// EXSurfaceCurrent:  3D surface represented as a current
//----------------------------------------------------------------

#include "surfacecurrent.h"

#include <math.h>				// For pow
#include <stdlib.h>				// For free

#include "gridoptimize.h"		// To make fast kernel computations

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EXSurfaceCurrent
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
SurfaceCurrent::SurfaceCurrent()
{
    this->_ny = 0;
    this->_nvy = 0;
    this->_nvx = 0;
    this->_sigmaW = 0;
    this->_timept = 0.0;
    this->_timeIndex = -1;
    this->_weight = 0.0;
    this->_grid = 0;
}

//----------------------------------------------------------------
// EXSurfaceCurrent
//----------------------------------------------------------------
// Inputs: 
//   y - a 2D array (dim, num_pts) of target points
//   vy - a 2D array (dim, num_pts) of target triangles
//   vx - a 2D array (dim, num_pts) of source trianges which
//        correspond to this target
//   sigmaW - the scale of the norm on currents
//   timept - the time corresponding to this target in whatever
//            time units used in the experiment
//   timeIndex - the index into the time discretization array 
//
// Outputs: 
//----------------------------------------------------------------
// Init constructor
//----------------------------------------------------------------
SurfaceCurrent::SurfaceCurrent(const Array2D<double>& y, const Array2D<int>& vy, const Array2D<int>& vx, double sigmaW, double timept, int timeIndex)
{
    this->_ny = y.GetWidth();
    this->_y = y;
    this->_nvy = vy.GetWidth();
    this->_vy = vy;
    this->_nvx = vx.GetWidth();
    this->_vx = vx;
    this->_sigmaW = sigmaW;
    this->_timept = timept;
    this->_timeIndex = timeIndex;
    this->_weight = 1.0f;

    this->_grid = new Grid();

    Array2D<double> tempCenters(3,this->_nvy);
    Array2D<double> tempNormals(3,this->_nvy);
    this->ComputeCentersAndNormals(y, vy, tempCenters, tempNormals);

    this->_centers = tempCenters;
    this->_normals = tempNormals;
}

//----------------------------------------------------------------
// EXSurfaceCurrent
//----------------------------------------------------------------
// Inputs: 
//   y - a 2D array (dim, num_pts) of target points
//   vy - a 2D array (dim, num_pts) of target triangles
//   vx - a 2D array (dim, num_pts) of source trianges which
//        correspond to this target
//   sigmaW - the scale of the norm on currents
//   timept - the time corresponding to this target in whatever
//            time units used in the experiment
//   timeIndex - the index into the time discretization array 
//   weight - the importance of this target in the regression
//
// Outputs: 
//----------------------------------------------------------------
// Init constructor - Recommended constructor
//----------------------------------------------------------------
SurfaceCurrent::SurfaceCurrent(const Array2D<double>& y, const Array2D<int>& vy, const Array2D<int>& vx, double sigmaW, double timept, int timeIndex, double weight)
{
    this->_ny = y.GetWidth();
    this->_y = y;
    this->_nvy = vy.GetWidth();
    this->_vy = vy;
    this->_nvx = vx.GetWidth();
    this->_vx = vx;
    this->_sigmaW = sigmaW;
    this->_timept = timept;
    this->_timeIndex = timeIndex;
    this->_weight = weight;

    this->_grid = new Grid();

    Array2D<double> tempCenters(3,this->_nvy);
    Array2D<double> tempNormals(3,this->_nvy);
    this->ComputeCentersAndNormals(y, vy, tempCenters, tempNormals);

    this->_centers = tempCenters;
    this->_normals = tempNormals;
}

//----------------------------------------------------------------
// EXSurfaceCurrent
//----------------------------------------------------------------
// Inputs: 
//   shape - EXSurfaceCurrent object to copy
//
// Outputs:
//----------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------
SurfaceCurrent::SurfaceCurrent(const SurfaceCurrent& shape)
{
    this->_ny = shape._ny;
    this->_y = shape._y;
    this->_nvy = shape._nvy;
    this->_vy = shape._vy;
    this->_nvx = shape._nvx;
    this->_vx = shape._vx;
    this->_sigmaW = shape._sigmaW;
    this->_timept = shape._timept;
    this->_timeIndex = shape._timeIndex;
    this->_centers = shape._centers;
    this->_normals = shape._normals;
    this->_grid = shape._grid;
}

//----------------------------------------------------------------
// ~EXSurfaceCurrent
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
SurfaceCurrent::~SurfaceCurrent()
{
    delete this->_grid;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// HELPER FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// ComputeCentersAndNormals
//----------------------------------------------------------------
// Inputs: 
//   pts - 2D array (dim, num_pts) of data points
//   tris - 2D array (dim, num_pts) of data triangles
//   centers - 2D array (dim, num_pts) for face centers 
//             return value
//   normals - 2D array (dim, num_pts) for face normals 
//             return value
//
// Outputs:
//----------------------------------------------------------------
// Computes face centers and normals given points and face 
// connectivity.  Returns centers and normals through reference 
// parameters.
//----------------------------------------------------------------
void SurfaceCurrent::ComputeCentersAndNormals(const Array2D<double>& pts, const Array2D<int>& tris, Array2D<double> &centers, Array2D<double> &normals)
{
    int numTris = tris.GetWidth();

    // Loop over all the faces
    //#pragma omp parallel for
    for (int f=0; f<numTris; f++)
    {
        double v[9];
        for (int k=0; k<3; k++)
        {
            for (int j=0; j<3; j++)
            {
                v[k+3*j] = pts(k,tris(j,f));
            }
        }

        // Calculate the center of each face c = (v1+v2+v3)/3;
        centers(0,f) = ((double)(v[0]+v[3]+v[6]))/3.0f;
        centers(1,f) = ((double)(v[1]+v[4]+v[7]))/3.0f;
        centers(2,f) = ((double)(v[2]+v[5]+v[8]))/3.0f;

        // Calcuate the normal of each face N = [(v2-v1)a(v3-v1)]/2;
        normals(0,f) = ((double)(((v[4]-v[1])*(v[8]-v[2])-(v[5]-v[2])*(v[7]-v[1]))))/2.0f;
        normals(1,f) = ((double)(((v[5]-v[2])*(v[6]-v[0])-(v[3]-v[0])*(v[8]-v[2]))))/2.0f;
        normals(2,f) = ((double)(((v[3]-v[0])*(v[7]-v[1])-(v[4]-v[1])*(v[6]-v[0]))))/2.0f;
    }
}

//----------------------------------------------------------------
// SpecProd
//----------------------------------------------------------------
// Inputs: 
//   A - 2D array (dim, num_pts)
//   B - 2D array (dim, num_pts) 
//   m - some integer for reshaping arrays
//
// Outputs:
//   returns - a 2D array (dim, num_pts)
//----------------------------------------------------------------
// Does some reshaping of arrays, then performs some matrix
// multiplications then reshapes the output.
//----------------------------------------------------------------
Array2D<double> SurfaceCurrent::SpecProd(const Array2D<double>& A, const Array2D<double>& B, int m)
{
    // Get dimensions of input arrays
    int n = A.GetWidth();
    int a = A.GetLength();
    int b = B.GetLength();

    // Storage for reshaped input arrays
    Array3D<double> AM(a/m,m,n);
    Array3D<double> BM(m,b/m,n);

    // Reshape A
    for (int i=0; i<a/m; i++)
    {
        for (int j=0; j<m; j++)
        {
            for (int k=0; k<n; k++)
            {
                int index = j*(a/m)+i;
                AM(i,j,k) = A(index,k);
            }
        }
    }

    // Reshape B
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<b/m; j++)
        {
            for (int k=0; k<n; k++)
            {
                int index = j*m+i;
                BM(i,j,k) = B(index,k);
            }
        }
    }

    // Storage for matrix multiplication
    Array3D<double> C(a/m,b/m,n);
    C.FillArray(0.0f);

    // Do matrix multiplication
    //#pragma omp parallel for collapse(4)
    for (int p=0; p<n; p++)
    {
        for (int i=0; i<a/m; i++)
        {
            for (int j=0; j<b/m; j++)
            {
                for (int k=0; k<m; k++)
                {
                    C(i,j,p) += AM(i,k,p)*BM(k,j,p);
                }
            }
        }
    }

    // Storage for reshaped output array
    Array2D<double> output((a*b/m/m),n);

    // Reshape the output
    for (int i=0; i<a/m; i++)
    {
        for (int j=0; j<b/m; j++)
        {
            for (int k=0; k<n; k++)
            {
                int index = j*(a/m) + i;
                output(index,k) = C(i,j,k);
            }
        }
    }

    return output;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// COMPUTE DATA MATCHING METRIC AND GRADIENT OF DATA MATCHING
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// RegScalarProduct
//----------------------------------------------------------------
// Inputs: 
//   centA - 2D array (dim, num_pts) of centers 
//   normA - 2D array (dim, num_pts) of normals
//   centB - 2D array (dim, num_pts) of centers 
//   normB - 2D array (dim, num_pts) of normals
//   sigmaW2 - squared scale of norm on currents
//
// Outputs:
//   returns - the value of data matching between A and B
//----------------------------------------------------------------
// Computes the norm on currents between A and B using exact
// kernel computations
//----------------------------------------------------------------
double SurfaceCurrent::RegScalarProduct(const Array2D<double>& centA, const Array2D<double>& normA,
                                        const Array2D<double>& centB, const Array2D<double>& normB, double sigmaW2)
{
    int numFacesA = centA.GetWidth();
    int numFacesB = centB.GetWidth();
    double scalarProd = 0.0f;

    // Loop over all faces for shape A
    //#pragma omp parallel for collapse(2) reduction(+:scalarProd)
    for (int i=0; i<numFacesA; i++)
    {
        // Loop over all faces for shape B
        for (int j=0; j<numFacesB; j++)
        {
            // The value to be used in the exp function -[(ca-cb)^2/sigmaW^2]
            double argin = -(pow(centA(0,i)-centB(0,j),2) +
                             pow(centA(1,i)-centB(1,j),2) +
                             pow(centA(2,i)-centB(2,j),2)) / sigmaW2;

            // Compute the exp
            double argout = exp(argin);

            // Update the sum of scalar products
            scalarProd += argout * (normA(0,i)*normB(0,j) +
                                    normA(1,i)*normB(1,j) +
                                    normA(2,i)*normB(2,j));
        }
    }

    return scalarProd;
}

//----------------------------------------------------------------
// GridScalarProduct
//----------------------------------------------------------------
// Inputs: 
//   centA - 2D array (dim, num_pts) of centers 
//   normA - 2D array (dim, num_pts) of normals
//   centB - 2D array (dim, num_pts) of centers 
//   normB - 2D array (dim, num_pts) of normals
//   sigmaW2 - squared scale of norm on currents
//
// Outputs:
//   returns - the value of data matching between A and B
//----------------------------------------------------------------
// Computes the norm on currents between A and B using 
// approximate kernel computations using a grid based method
//----------------------------------------------------------------
double SurfaceCurrent::GridScalarProduct(const Array2D<double>& centA, const Array2D<double>& normA,
                                         const Array2D<double>& centB, const Array2D<double>& normB, double)
{
    int numFacesA = centA.GetWidth();
    int numFacesB = centB.GetWidth();

    // Storage for column major linearized data for the grid optimizer
    double* ca = new double[3*numFacesA];
    double* na = new double[3*numFacesA];
    double* cb = new double[3*numFacesB];

    // Store linearly in column major order
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<numFacesA; j++)
        {
            int index = j*3 + i;
            ca[index] = centA(i,j);
            na[index] = normA(i,j);
        }
    }

    // Store linearly in column major order
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<numFacesB; j++)
        {
            int index = j*3 + i;
            cb[index] = centB(i,j);
        }
    }

    // Get grid information
    int* gridlong = this->GetGrid()->GetLong();
    double gridpas = this->GetGrid()->GetPas();
    double* gridorigin = this->GetGrid()->GetOrigin();
    double* gridfft3kd = this->GetGrid()->GetFFT3KD();
    // Storage for column major linearied grid fft
    double* colfft3kd = new double[((gridlong[0]/2)+1)*gridlong[1]*gridlong[2]];

    // Store linearly in column major order
    for (int i=0; i<((gridlong[0]/2)+1); i++)
    {
        for (int j=0; j<gridlong[1]; j++)
        {
            for (int k=0; k<gridlong[2]; k++)
            {
                int indexRM = i*gridlong[2]*gridlong[1] + j*gridlong[2] + k;
                int indexCM = k*((gridlong[0]/2)+1)*gridlong[1] + j*((gridlong[0]/2)+1) + i;

                colfft3kd[indexCM] = gridfft3kd[indexRM];
            }
        }
    }

    // Perform the kernel computations
    GridOptimize gridOptim(numFacesA, ca, 3, na, numFacesB, cb, gridlong, gridpas, gridorigin, colfft3kd);

    double* gridOut = gridOptim.GridOptim(numFacesA, ca, 3, na, numFacesB, cb, gridlong, gridpas, gridorigin, colfft3kd);

    //delete gridOptim;

    // The final value of the data matching metric is the dot product between normB and the kernel computations
    double sum = 0.0;
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<numFacesB; j++)
        {
            int index = j*3 + i;
            sum += normB(i,j)*gridOut[index];
        }
    }

    // Clean up memory
    //free(gridOut);
    //delete [] gridOut;
    delete [] ca;
    delete [] na;
    delete [] cb;
    delete [] colfft3kd;

    return sum;
}

//----------------------------------------------------------------
// RegGradNorm
//----------------------------------------------------------------
// Inputs: 
//   x -  2D array (dim, num_pts) of points to compare with the
//        target points 
//   centX - 2D array (dim, num_pts) of centers of faces to
//           compare with target faces
//   normX - 2D array (dim, num_pts) of normals of faces to
//           compare with target faces
//   sigmaW2 - squared scale of norm on currents
//
// Outputs:
//   returns - the gradient of data matching between A and B
//----------------------------------------------------------------
// Computes the gradient of the norm on currents using exact
// kernel computations
//----------------------------------------------------------------
Array2D<double> SurfaceCurrent::RegGradNorm(const Array2D<double>& x, const Array2D<double>& centX, const Array2D<double>& normX, double sigmaW2)
{
    int numPtsX = x.GetWidth();
    int numFacesX = centX.GetWidth();

    // The final gradient value grad_x_i for each point x_i
    Array2D<double> eta(3,numPtsX);
    eta.FillArray(0.0f);

    // Stores the final calculation for each point for the current face
    double etaloc[9];
    memset(etaloc, 0.0f, 9*sizeof(double));
    // This is the value of the first summation in equation 4.3.9 in Durrleman thesis
    double etafW[3];
    memset(etafW, 0.0f, 3*sizeof(double));
    // This is the value of the third summation in equation 4.3.9 in Durrleman thesis
    double etafDW[3];
    memset(etafDW, 0.0f, 3*sizeof(double));
    // All points (3 points, 9 values) for the current face
    double v[9];
    memset(v, 0.0f, 9*sizeof(double));

    // Loop over the source faces
    //#pragma omp critical
    for (int e=0; e<numFacesX; e++)
    {
        // Zero out the summation arrays
        double tempEtafWx = 0.0, tempEtafWy = 0.0, tempEtafWz = 0.0;
        double tempEtafDWx = 0.0, tempEtafDWy = 0.0, tempEtafDWz = 0.0;
        //memset(etafW, 0.0f, 3*sizeof(double));
        //memset(etafDW, 0.0f, 3*sizeof(double));

        // Loop over the source faces (compute the 1st and 3rd sums in equation 4.3.9 from Durrleman thesis)
        //#pragma omp parallel for reduction(+:tempEtafWx, tempEtafWy, tempEtafWz, tempEtafDWx, tempEtafDWy, tempEtafDWz)
        for (int f=0; f<numFacesX; f++)
        {
            // Calculates kernel portion for source centers from equation 4.3.9 in Durrleman thesis
            // Sum_e [exp(-(ca_e - ca_f)^2/sigmaW^2))]
            double argin = -( pow(centX(0,f)-centX(0,e),2) +
                              pow(centX(1,f)-centX(1,e),2) +
                              pow(centX(2,f)-centX(2,e),2) ) / sigmaW2;

            double argout = exp(argin);

            // Multiply argout by the normal (K(ck,cf) n_k) from equation 4.3.9 in Durrleman thesis
            tempEtafWx += argout*normX(0,f);
            tempEtafWy += argout*normX(1,f);
            tempEtafWz += argout*normX(2,f);
            //etafW[0] += argout*normX(0,f);
            //etafW[1] += argout*normX(1,f);
            //etafW[2] += argout*normX(2,f);

            // Second evaluation seems unnecessary
            argout = exp(argin);

            // The dot product of source normals (n_f^t n_k) from the 3rd summation in the Durrleman thesis
            double dotN = normX(0,f)*normX(0,e) + normX(1,f)*normX(1,e) + normX(2,f)*normX(2,e);

            // Multiply by dot product of normals and gradient_1
            tempEtafDWx += argout*dotN*(centX(0,f)-centX(0,e));
            tempEtafDWy += argout*dotN*(centX(1,f)-centX(1,e));
            tempEtafDWz += argout*dotN*(centX(2,f)-centX(2,e));
            //etafDW[0] += argout*dotN*(centX(0,f)-centX(0,e));
            //etafDW[1] += argout*dotN*(centX(1,f)-centX(1,e));
            //etafDW[2] += argout*dotN*(centX(2,f)-centX(2,e));
        }

        // Loop over the target faces (compute the 2nd and 4th sums in equation 4.3.9 from Durrleman thesis)
        //#pragma omp parallel for reduction(+:tempEtafWx, tempEtafWy, tempEtafWz, tempEtafDWx, tempEtafDWy, tempEtafDWz)
        for (int f=0; f<this->_nvy; f++)
        {
            // Calculates kernel portion bettween source/target centers from equation 4.3.9 in Durrleman thesis
            // Sum_e [exp(-(cb_e - ca_f)^2/sigmaW^2))]
            double argin = -( pow(this->_centers(0,f)-centX(0,e),2) +
                              pow(this->_centers(1,f)-centX(1,e),2) +
                              pow(this->_centers(2,f)-centX(2,e),2) ) / sigmaW2;

            double argout = exp(argin);

            // Multiply argout by the normal (K(ck,cf) n_k) from equation 4.3.9 in Durrleman thesis
            tempEtafWx -= argout*this->_normals(0,f);
            tempEtafWy -= argout*this->_normals(1,f);
            tempEtafWz -= argout*this->_normals(2,f);
            //etafW[0] -= argout*this->_normals(0,f);
            //etafW[1] -= argout*this->_normals(1,f);
            //etafW[2] -= argout*this->_normals(2,f);

            // Second evaluation seems unnecessary
            argout = exp(argin);

            // The dot product of source normals (n_f^t n'_k) from the 4th summation in the Durrleman thesis
            double dotN = this->_normals(0,f)*normX(0,e) + this->_normals(1,f)*normX(1,e) + this->_normals(2,f)*normX(2,e);

            // Multiply by dot product of normals and gradient_1
            tempEtafDWx -= argout*dotN*(this->_centers(0,f)-centX(0,e));
            tempEtafDWy -= argout*dotN*(this->_centers(1,f)-centX(1,e));
            tempEtafDWz -= argout*dotN*(this->_centers(2,f)-centX(2,e));
            //etafDW[0] -= argout*dotN*(this->_centers(0,f)-centX(0,e));
            //etafDW[1] -= argout*dotN*(this->_centers(1,f)-centX(1,e));
            //etafDW[2] -= argout*dotN*(this->_centers(2,f)-centX(2,e));
        }

        etafW[0] = tempEtafWx;
        etafW[1] = tempEtafWy;
        etafW[2] = tempEtafWz;
        etafDW[0] = tempEtafDWx;
        etafDW[1] = tempEtafDWy;
        etafDW[2] = tempEtafDWz;

        // 4/3 because we've combined the 2nd and 4th summations (not sure about the factor of 1/sigmaW^2 though...)
        etafDW[0] = 4.0f*etafDW[0]/sigmaW2/3.0f;
        etafDW[1] = 4.0f*etafDW[1]/sigmaW2/3.0f;
        etafDW[2] = 4.0f*etafDW[2]/sigmaW2/3.0f;

        // Find all points of this face
        for (int k=0; k<3; k++)
        {
            for (int j=0; j<3; j++)
            {
                v[k+3*j] = x(k,this->_vx(j,e));
            }
        }

        // Finish the gradient calculation for this point by computing the cross product of edges with values already calculated
        etaloc[0] = (v[4]-v[7])*etafW[2]-(v[5]-v[8])*etafW[1] + etafDW[0];
        etaloc[1] = (v[5]-v[8])*etafW[0]-(v[3]-v[6])*etafW[2] + etafDW[1];
        etaloc[2] = (v[3]-v[6])*etafW[1]-(v[4]-v[7])*etafW[0] + etafDW[2];
        etaloc[3] = (v[7]-v[1])*etafW[2]-(v[8]-v[2])*etafW[1] + etafDW[0];
        etaloc[4] = (v[8]-v[2])*etafW[0]-(v[6]-v[0])*etafW[2] + etafDW[1];
        etaloc[5] = (v[6]-v[0])*etafW[1]-(v[7]-v[1])*etafW[0] + etafDW[2];
        etaloc[6] = (v[1]-v[4])*etafW[2]-(v[2]-v[5])*etafW[1] + etafDW[0];
        etaloc[7] = (v[2]-v[5])*etafW[0]-(v[0]-v[3])*etafW[2] + etafDW[1];
        etaloc[8] = (v[0]-v[3])*etafW[1]-(v[1]-v[4])*etafW[0] + etafDW[2];

        // Find out the indicies of the points we just computed so we can update the gradient array in the correct locations
        int loc1 = this->_vx(0,e);
        int loc2 = this->_vx(1,e);
        int loc3 = this->_vx(2,e);

        // Keep a running sum at the points we just calculated, since the gradient at a point x_i
        // is the sum of the gradient of every face that has x_i as a vertex
        eta(0,loc1) += etaloc[0];
        eta(1,loc1) += etaloc[1];
        eta(2,loc1) += etaloc[2];
        eta(0,loc2) += etaloc[3];
        eta(1,loc2) += etaloc[4];
        eta(2,loc2) += etaloc[5];
        eta(0,loc3) += etaloc[6];
        eta(1,loc3) += etaloc[7];
        eta(2,loc3) += etaloc[8];

    }

    return eta;
}

//----------------------------------------------------------------
// GridGradNorm
//----------------------------------------------------------------
// Inputs: 
//   x -  2D array (dim, num_pts) of points to compare with the
//        target points 
//   centX - 2D array (dim, num_pts) of centers of faces to
//           compare with target faces
//   normX - 2D array (dim, num_pts) of normals of faces to
//           compare with target faces
//   sigmaW2 - squared scale of norm on currents
//
// Outputs:
//   returns - the gradient of data matching between A and B
//----------------------------------------------------------------
// Computes the gradient of the norm on currents using exact
// approximate kernel computations using a grid based method
//----------------------------------------------------------------
Array2D<double> SurfaceCurrent::GridGradNorm(const Array2D<double>& x, const Array2D<double>& centX, const Array2D<double>& normX, double sigmaW2)
{
    int numPtsX = x.GetWidth();
    int numFacesX = centX.GetWidth();

    // Concatinate centers of x and y
    Array2D<double> centers(3,numFacesX+this->_nvy);
    for (int j=0; j<3; j++)
    {
        for (int i=0; i<numFacesX+this->_nvy; i++)
        {
            if (i<numFacesX)
            {
                centers(j,i) = centX(j,i);
            }
            else
            {
                centers(j,i) = this->_centers(j,i-numFacesX);
            }
        }
    }

    // We must concatinate the normals and flip the direction of the target
    Array2D<double> normals(3,numFacesX+this->_nvy);
    for (int j=0; j<3; j++)
    {
        for (int i=0; i<numFacesX+this->_nvy; i++)
        {
            if (i<this->_nvx)
            {
                normals(j,i) = normX(j,i);
            }
            else
            {
                normals(j,i) = -1.0f*this->_normals(j,i-numFacesX);
            }
        }
    }

    // Do some matrix multiplication and reshaping
    Array2D<double> sp = this->SpecProd(centers, normals, 1);

    // Now we must concatinate the centers and the spec prod
    Array2D<double> concatSPNorms(12,numFacesX+this->_nvy);
    for (int i=0; i<12; i++)
    {
        for (int j=0; j<numFacesX+this->_nvy;j++)
        {
            if (i<3)
            {
                concatSPNorms(i,j) = normals(i,j);
            }
            else
            {
                concatSPNorms(i,j) = sp(i-3,j);
            }
        }
    }

    // Get grid information
    int* gridlong = this->GetGrid()->GetLong();
    double gridpas = this->GetGrid()->GetPas();
    double* gridorigin = this->GetGrid()->GetOrigin();
    double* gridfft3kd = this->GetGrid()->GetFFT3KD();

    // Linearize fft3kd in column major order
    double* colfft3kd = new double[((gridlong[0]/2)+1)*gridlong[1]*gridlong[2]];
    for (int i=0; i<((gridlong[0]/2)+1); i++)
    {
        for (int j=0; j<gridlong[1]; j++)
        {
            for (int k=0; k<gridlong[2]; k++)
            {
                int indexRM = i*gridlong[2]*gridlong[1] + j*gridlong[2] + k;
                int indexCM = k*((gridlong[0]/2)+1)*gridlong[1] + j*((gridlong[0]/2)+1) + i;

                colfft3kd[indexCM] = gridfft3kd[indexRM];
            }
        }
    }

    // Linearize centers in column major order
    double* czlin = new double[3*(numFacesX+this->_nvy)];
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<(numFacesX+this->_nvy); j++)
        {
            int index = j*3 + i;
            czlin[index] = centers(i,j);
        }
    }

    // Linearize X centers in column major order
    double* cphilin = new double[3*numFacesX];
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<numFacesX; j++)
        {
            int index = j*3 + i;
            cphilin[index] = centX(i,j);
        }
    }

    // Linearize the middle term in column major order
    double* middleTerm = new double[12*(numFacesX+this->_nvy)];
    for (int i=0; i<12; i++)
    {
        for (int j=0; j<numFacesX+this->_nvy; j++)
        {
            int index = j*12 + i;
            middleTerm[index] = concatSPNorms(i,j);
        }
    }

    GridOptimize gridOptim(numFacesX+this->_nvy, czlin, 12, middleTerm, numFacesX, cphilin, gridlong, gridpas, gridorigin, colfft3kd);

    double* gridOut = gridOptim.GridOptim(numFacesX+this->_nvy, czlin, 12, middleTerm, numFacesX, cphilin, gridlong, gridpas, gridorigin, colfft3kd);

    //delete gridOptim;

    Array2D<double> go0to3(3,numFacesX);
    Array2D<double> go4toEnd(9,numFacesX);
    for (int i=0; i<12; i++)
    {
        for (int j=0; j<numFacesX; j++)
        {
            int index = j*12+ i;

            if (i<3)
            {
                go0to3(i,j) = gridOut[index];
            }
            else
            {
                go4toEnd(i-3,j) = gridOut[index];
            }
        }
    }

    // Compute several specprods
    Array2D<double> spgo4ToEndXNorms = this->SpecProd(go4toEnd, normX, 3);
    Array2D<double> spXCentsXNorms = this->SpecProd(centX, normX, 1);
    Array2D<double> spspgo0to3 = this->SpecProd(spXCentsXNorms, go0to3, 3);

    Array2D<double> eta_dw(3,numFacesX);

    for (int i=0; i<3; i++)
    {
        for (int j=0; j<numFacesX; j++)
        {
            eta_dw(i,j) = (4.0f/3.0f/sigmaW2) * (spgo4ToEndXNorms(i,j) - spspgo0to3(i,j));
        }
    }

    Array2D<double> eta(3,numPtsX);
    eta.FillArray(0.0f);

    // Finish up the computation of the gradient
    for (int m=0; m<numFacesX; m++)
    {
        // Get the index to the current vertex
        int vm1 = this->_vx(0,m);
        int vm2 = this->_vx(1,m);
        int vm3 = this->_vx(2,m);

        double e1x = x(0,vm2) - x(0,vm3);
        double e1y = x(1,vm2) - x(1,vm3);
        double e1z = x(2,vm2) - x(2,vm3);

        double e2x = x(0,vm3) - x(0,vm1);
        double e2y = x(1,vm3) - x(1,vm1);
        double e2z = x(2,vm3) - x(2,vm1);

        double e3x = x(0,vm1) - x(0,vm2);
        double e3y = x(1,vm1) - x(1,vm2);
        double e3z = x(2,vm1) - x(2,vm2);

        eta(0,vm1) += eta_dw(0,m) + (e1y*go0to3(2,m) - e1z*go0to3(1,m));
        eta(1,vm1) += eta_dw(1,m) + (e1z*go0to3(0,m) - e1x*go0to3(2,m));
        eta(2,vm1) += eta_dw(2,m) + (e1x*go0to3(1,m) - e1y*go0to3(0,m));

        eta(0,vm2) += eta_dw(0,m) + (e2y*go0to3(2,m) - e2z*go0to3(1,m));
        eta(1,vm2) += eta_dw(1,m) + (e2z*go0to3(0,m) - e2x*go0to3(2,m));
        eta(2,vm2) += eta_dw(2,m) + (e2x*go0to3(1,m) - e2y*go0to3(0,m));

        eta(0,vm3) += eta_dw(0,m) + (e3y*go0to3(2,m) - e3z*go0to3(1,m));
        eta(1,vm3) += eta_dw(1,m) + (e3z*go0to3(0,m) - e3x*go0to3(2,m));
        eta(2,vm3) += eta_dw(2,m) + (e3x*go0to3(1,m) - e3y*go0to3(0,m));
    }

    // Clean up memory
    //free(gridOut);
    //delete [] gridOut;
    delete [] czlin;
    delete [] middleTerm;
    delete [] cphilin;
    delete [] colfft3kd;

    return eta;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// INTERFACE TO COMPUTE DATA MATCHING METRIC AND GRADIENT
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// Matching
//----------------------------------------------------------------
// Inputs: 
//   x - a 3D array (dim, num_pts, time) of points to compare 
//       with the target points
//   t - the time index into x to use
//
// Outputs:
//----------------------------------------------------------------
// Computes the data matching metric between sets of surfaces
// represented as currents
//----------------------------------------------------------------
double SurfaceCurrent::Matching(const Array3D<double>& x, int t)
{
    int nx = x.GetWidth();
    double match = 0.0f;

    // Extract the portion of x we are interested in
    Array2D<double> xt(3,nx);
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<nx; j++)
        {
            xt(i,j) = x(i,j,t);
        }
    }

    // Compute centers and normals for input data x
    Array2D<double> centers(3,this->_nvx);
    Array2D<double> normals(3,this->_nvx);
    this->ComputeCentersAndNormals(xt, this->_vx, centers, normals);

    // We must concatinate the centers
    Array2D<double> matchCenters(3,this->_nvx+this->_nvy);
    for (int j=0; j<3; j++)
    {
        for (int i=0; i<this->_nvx+this->_nvy; i++)
        {
            if (i<this->_nvx)
            {
                matchCenters(j,i) = centers(j,i);
            }
            else
            {
                matchCenters(j,i) = this->_centers(j,i-this->_nvx);
            }
        }
    }

    // We must concatinate the normals and flip the direction of the target
    Array2D<double> matchNormals(3,this->_nvx+this->_nvy);
    for (int j=0; j<3; j++)
    {
        for (int i=0; i<this->_nvx+this->_nvy; i++)
        {
            if (i<this->_nvx)
            {
                matchNormals(j,i) = normals(j,i);
            }
            else
            {
                matchNormals(j,i) = -1.0f*this->_normals(j,i-this->_nvx);
            }
        }
    }

    // Compute the scalar product in the space of currents
    double sigmaW2 = pow(this->GetSigmaW(),2);

    if (strcmp(this->GetKernelType(), "p3m") == 0)
    {
        match = this->GridScalarProduct(matchCenters, matchNormals, matchCenters, matchNormals, sigmaW2);
    }
    else
    {
        match = this->RegScalarProduct(matchCenters, matchNormals, matchCenters, matchNormals, sigmaW2);
    }

    return match;
}

//----------------------------------------------------------------
// GradMatching
//----------------------------------------------------------------
// Inputs: 
//   x - a 3D array (dim, num_pts, time) of points to compare 
//       with the target points
//   t - the time index into x to use
//
// Outputs:
//----------------------------------------------------------------
// Computes the gradient of the data matching metric between 
// surfaces represented as currents
//----------------------------------------------------------------
Array2D<double> SurfaceCurrent::GradMatching(const Array3D<double>& x, int t)
{
    int nx = x.GetWidth();
    Array2D<double> g(3,nx);

    // Extract the portion of x we are interested in
    Array2D<double> xt(3, nx);
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<nx; j++)
        {
            xt(i,j) = x(i,j,t);
        }
    }

    // Compute centers and normals for input data x
    Array2D<double> centers(3,this->_nvx);
    Array2D<double> normals(3,this->_nvx);
    this->ComputeCentersAndNormals(xt, this->_vx, centers, normals);

    // Compute the gradient of the norm in the space of currents
    double sigmaW2 = pow(this->GetSigmaW(),2);

    if (strcmp(this->GetKernelType(), "p3m") == 0)
    {
        g = this->GridGradNorm(xt, centers, normals, sigmaW2);

    }
    else
    {
        g = this->RegGradNorm(xt, centers, normals, sigmaW2);
    }

    return g;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// OVERLOADED OPERATORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// operator = 
//----------------------------------------------------------------
// Inputs: 
//   shape - EXSurfaceCurrent object to copy
//
// Outputs:
//   return - copied EXSurfaceCurrent object 
//----------------------------------------------------------------
// Overloaded assignment operator
//----------------------------------------------------------------
SurfaceCurrent& SurfaceCurrent::operator = (const SurfaceCurrent& shape)
{
    if (this != &shape)
    {
        this->_ny = shape._ny;
        this->_y = shape._y;
        this->_nvy = shape._nvy;
        this->_vy = shape._vy;
        this->_nvx = shape._nvx;
        this->_vx = shape._vx;
        this->_sigmaW = shape._sigmaW;
        this->_timept = shape._timept;
        this->_timeIndex = shape._timeIndex;
        this->_centers = shape._centers;
        this->_normals = shape._normals;
        this->_grid = shape._grid;
    }

    return *this;
}
