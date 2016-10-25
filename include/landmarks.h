//----------------------------------------------------------------
// Landmarks: Target data represented as landmarks
//----------------------------------------------------------------

#ifndef LANDMARKS_H
#define LANDMARKS_H

#include "targetdata.h"

class Landmarks : public TargetData
{

    public:
        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------

        Array2D<double> _centers;		// Centers of target faces
        Array2D<double> _normals;		// Normals of target faces

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        Landmarks();
        //EXLandmarks(const EX2DArray<double>& y, const EX2DArray<int>& vx, double timept, int timeIndex);
        Landmarks(const Array2D<double>& y, const Array2D<int>& vy, const Array2D<int>& vx, double sigmaW, double timept, int timeIndex, double weight);
        Landmarks(const Landmarks& shape);
        ~Landmarks();

        //--------------------------------------------------------------------------
        // Interface to compute data matching metric and gradient
        //--------------------------------------------------------------------------

        virtual double Matching(const Array3D<double>& x, int t);
        virtual Array2D<double> GradMatching(const Array3D<double>& x, int t);

        //--------------------------------------------------------------------------
        // Overloaded operators
        //--------------------------------------------------------------------------
        Landmarks& operator = (const Landmarks& shape);

};

#endif // LANDMARKS_H
