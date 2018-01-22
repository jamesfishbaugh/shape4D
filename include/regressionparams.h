//--------------------------------------------------------------------------
// RegressionParams:  Shape regression params
//--------------------------------------------------------------------------

#ifndef REGRESSIONPARAMS_H
#define REGRESSIONPARAMS_H

#include "array1d.h"
#include "array2d.h"
#include "grid.h"

class RegressionParams
{

    private:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------

        int _nx;						    // The number of source points
        Array2D<double> _x;                 // Source points array 3 x nx
        double _sigmaV;				        // Scale of deformation
        double _gamma;					    // Weight on regularity
        double _stdV;					    //
        int _T;							    // Number of time steps for discretization
        Array1D<double> _xtime;             // Time discretization array
        Array2D<double> _initV0;			// Initial velocity
        bool _useInitV0	;				    // Should we only use the initV0?
        bool _estimateBaseline;				// Should we estimate the baseline shape or leave it fixed?
        double _v0Weight;				    // The importance of velocity in optimization
        char* _kernelType;					// The kernel type
        int _maxIters;                      // Max iterations for gradient descent
        double _breakRatio;                 // Break ratio for optimzation
        bool _useFista;                     // Should we use FISTA
        double _baselineSmoothing;          // Factor for baseline shape smoothing

        Grid** _grids;                      // An array of grids at each time point

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        RegressionParams();
        RegressionParams(const Array2D<double>& x, double sigmaV, double gamma, const Array1D<double>& xtime);
        RegressionParams(const Array2D<double>& x, double sigmaV, double gamma, const Array1D<double>& xtime, const Array2D<double>& initV0);
        RegressionParams(const RegressionParams& source);
        ~RegressionParams();

        //--------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------

        void SetX(const Array2D<double>& x);
        void SetSigmaV(double sigmaV);
        void SetGamma(double gamma);
        void SetStdV(double stdV);
        void SetXTime(const Array1D<double>& xtime);
        void SetInitV0(const Array2D<double>& initV0);
        void SetV0Weight(double v0Weight);
        void SetShouldUseInitV0(bool yesNo);
        void SetShouldEstimateBaseline(bool yesNo);
        void SetKernelType(char* kernelType);
        void SetMaxIters(int maxIters);
        void SetBreakRatio(double breakRatio);
        void SetShouldUseFista(bool yesNo);
        void SetBaselineSmoothing(double baselineSmoothing);

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        int GetNx() const;
        const Array2D<double> GetX() const;
        double GetSigmaV() const;
        double GetGamma() const;
        double GetStdV() const;
        int GetT() const;
        const Array1D<double> GetXTime() const;
        bool ShouldUseInitV0() const;
        bool ShouldEstimateBaseline() const;
        const char* GetKernelType() const;
        int GetMaxIters() const;
        double GetBreakRatio() const;
        const Array2D<double> GetInitV0() const;
        double GetV0Weight() const;
        Grid* GetGridAt(int index);
        const Grid* GetGridAt(int index) const;
        bool ShouldUseFista() const;
        double GetBaselineSmoothing() const;

        //--------------------------------------------------------------------------
        // Overloaded operators
        //--------------------------------------------------------------------------

        RegressionParams& operator = (const RegressionParams& source);

};

#endif // REGRESSIONPARAMS_H
