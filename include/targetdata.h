//--------------------------------------------------------------------------
// TargetData:  Virtual base class representing target data for
// regression
//--------------------------------------------------------------------------

#ifndef TARGETDATA_H
#define TARGETDATA_H

// Constants that define the specific data representation
#define LANDMARKS    0
#define SURFACE      1

#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "grid.h"

class TargetData
{

    protected:

        //--------------------------------------------------------------------------
        // Protected member variables
        //--------------------------------------------------------------------------

        int _ny;					// Number of points of target
        Array2D<double> _y;		// Points of the target
        int _nvy;					// Number of triangles of target
        Array2D<int> _vy;			// Triangles of the target
        int _nvx;					// Number of triangles of source
        Array2D<int> _vx;			// Triangles of the source
        double _sigmaW;				// Size of the kernel
        double _timept;				// Time associated with the target
        int _timeIndex; 			// Time index associated with the target
        double _weight;				// Relative importance of the target
        char* _kernelType; 			// Kernel type for convolution

        Grid* _grid;				// A grid for this target

    public:

        //--------------------------------------------------------------------------
        // Destructor
        //--------------------------------------------------------------------------

        virtual ~TargetData() = 0;

        //--------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------

        void SetY(const Array2D<double>& y);
        void SetVy(const Array2D<int>& vy);
        void SetVx(const Array2D<int>& vx);
        void SetSigmaW(double sigmaW);
        void SetTimept(double timept);
        void SetTimeIndex(int timeIndex);
        void SetWeight(double weight);
        void SetKernelType(char* kernelType);

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        int GetNy() const;
        const Array2D<double> GetY() const;
        int GetNvy() const;
        const Array2D<int> GetVy() const;
        int GetNvx() const;
        const Array2D<int> GetVx() const;
        double GetSigmaW() const;
        double GetTimept() const;
        double GetTimeIndex() const;
        double GetWeight() const;
        const char* GetKernelType() const;
        Grid* GetGrid();
        const Grid* GetGrid() const;

        //--------------------------------------------------------------------------
        // Interface to compute data matching metric and gradient
        //--------------------------------------------------------------------------

        virtual double Matching(const Array3D<double>& x, int t) = 0;
        virtual Array2D<double> GradMatching(const Array3D<double>& x, int t) = 0;

};

#endif // TARGETDATA_H
