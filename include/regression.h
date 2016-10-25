//--------------------------------------------------------------------------
// Regression:  Virutal base class that implements common functionality
// for regression
//--------------------------------------------------------------------------

#ifndef REGRESSION_H
#define REGRESSION_H

#include "algorithm.h"
#include "regressionparams.h"
#include "targetdata.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"

#include "multiobjectcomplex.h"

using namespace std;

class Regression : public Algorithm
{

    protected:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------

        MultiObjectComplex _theSource;
        vector<MultiObjectComplex> _theTarget;

        RegressionParams _source;		// Source structure
        int _numTargets;				// Number of targets
        TargetData** _targets;          // Array of targets
        Array3D<double> _X;             // Trajectory of source points over time
        double _tau;					// Step size for integration schemes

        //--------------------------------------------------------------------------
        // Grid methods
        //--------------------------------------------------------------------------

        void SetSourceGrids();
        void UpdateSourceGrids();
        void SetTargetGrids();
        void UpdateTargetGrids();

        //--------------------------------------------------------------------------
        // Initialization methods
        //--------------------------------------------------------------------------

        void InitX();

        //--------------------------------------------------------------------------
        // Write methods
        //--------------------------------------------------------------------------
        virtual void WriteSelf() = 0;

        //--------------------------------------------------------------------------
        // Data matching and gradient of data matching
        //--------------------------------------------------------------------------

        double ComputeDataMatching(const Array3D<double>& X);
        Array3D<double> ComputeDataMatchingGradient(const Array3D<double>& X);

        //--------------------------------------------------------------------------
        // Kernel methods
        //--------------------------------------------------------------------------

        Array2D<double> RegKernelSum(int t, const Array3D<double>& weights);
        Array2D<double> GridKernelSum(int t, const Array3D<double>& weights);
        Array2D<double> RegCompDeta(int t, const Array3D<double>& eta, const Array3D<double>& momenta);
        Array2D<double> GridCompDeta(int t, const Array3D<double>& eta, const Array3D<double>& momenta);

    public:

        //--------------------------------------------------------------------------
        // Destructors
        //--------------------------------------------------------------------------

        virtual ~Regression() = 0;

        //--------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------

        void SetTheSource(const MultiObjectComplex theSource);
        void SetTargets(const vector<MultiObjectComplex>& theTargets);
        void AddTarget(const MultiObjectComplex& target);
        void SetSource(const RegressionParams& source);
        void SetNumTargets(int numTargets);
        void SetTargets(TargetData** targets);
        void SetTargets(int numTargets, TargetData** targets);

        //--------------------------------------------------------------------------
        // Methods for the optimizer
        //--------------------------------------------------------------------------

        virtual double ComputeFunctional(const Array3D<double>& parameter, double &dataMatching, double &regularity) = 0;
        virtual Array3D<double> ComputeGradient(const Array3D<double>& parameter) = 0;

        //--------------------------------------------------------------------------
        // Main algorithm method
        //--------------------------------------------------------------------------

        virtual Array3D<double> Run() = 0;
};

#endif // REGRESSION_H
