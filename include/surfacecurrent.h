//----------------------------------------------------------------
// SurfaceCurrent:  3D surface represented as a current
//----------------------------------------------------------------

#ifndef SURFACECURRENT_H
#define SURFACECURRENT_H

#include "targetdata.h"

class SurfaceCurrent : public TargetData
{

    private:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------

        Array2D<double> _centers;		// Centers of target faces
        Array2D<double> _normals;		// Normals of target faces

        //--------------------------------------------------------------------------
        // Compute data matching metric and gradient of data matching
        //--------------------------------------------------------------------------

        double RegScalarProduct(const Array2D<double>& centA, const Array2D<double>& normA,
                                const Array2D<double>& centB, const Array2D<double>& normB, double sigmaW2);
        double GridScalarProduct(const Array2D<double>& centA, const Array2D<double>& normA,
                                 const Array2D<double>& centB, const Array2D<double>& normB, double sigmaW2);
        Array2D<double> RegGradNorm(const Array2D<double>& x, const Array2D<double>& centX, const Array2D<double>& normX, double sigmaW2);
        Array2D<double> GridGradNorm(const Array2D<double>& x, const Array2D<double>& centX, const Array2D<double>& normX, double sigmaW2);

        //--------------------------------------------------------------------------
        // Helper functions
        //--------------------------------------------------------------------------

        void ComputeCentersAndNormals(const Array2D<double>& pts, const Array2D<int>& tris, Array2D<double> &centers, Array2D<double> &normals);
        Array2D<double> SpecProd(const Array2D<double>& A, const Array2D<double>& B, int m);

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        SurfaceCurrent();
        SurfaceCurrent(const Array2D<double>& y, const Array2D<int>& vy, const Array2D<int>& vx, double sigmaW, double timept, int timeIndex);
        SurfaceCurrent(const Array2D<double>& y, const Array2D<int>& vy, const Array2D<int>& vx, double sigmaW, double timept, int timeIndex, double weight);
        SurfaceCurrent(const SurfaceCurrent& shape);
        virtual ~SurfaceCurrent();

        //--------------------------------------------------------------------------
        // Interface to compute data matching metric and gradient
        //--------------------------------------------------------------------------

        virtual double Matching(const Array3D<double>& x, int t);
        virtual Array2D<double> GradMatching(const Array3D<double>& x, int t);

        //--------------------------------------------------------------------------
        // Overloaded operators
        //--------------------------------------------------------------------------
        SurfaceCurrent& operator = (const SurfaceCurrent& shape);

};

#endif // SURFACECURRENT_H
