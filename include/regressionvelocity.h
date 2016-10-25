//--------------------------------------------------------------------------
// RegressionVelocity:  Shape regression parameterized by velocity
//--------------------------------------------------------------------------

#ifndef REGRESSIONVELOCITY_H
#define REGRESSIONVELOCITY_H

#include "regression.h"

class RegressionVelocity : public Regression
{

    private:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------

        Array3D<double> _dX;            // Velocity vectors over time
        Array3D<double> _momenta;       // Momenta vectors over time

        int _iteration;
        bool _continueRegression;

        //--------------------------------------------------------------------------
        // Write methods
        //--------------------------------------------------------------------------
        void WriteSelf();

        //--------------------------------------------------------------------------
        // Computation of regression functional
        //--------------------------------------------------------------------------

        double ComputeRegularity(const Array3D<double>& momenta);

        //--------------------------------------------------------------------------
        // Helper methods
        //--------------------------------------------------------------------------

        void ComputeTrajectories(const Array3D<double>& momenta);

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        RegressionVelocity();
        RegressionVelocity(const RegressionParams& source);
        RegressionVelocity(const RegressionParams& source, int numTargets, TargetData** targets);
        virtual ~RegressionVelocity();

        //--------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------

        void SetMomenta(Array3D<double>& momenta);

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        Array2D<double> GetV0();
        Array3D<double> GetVelocity();
        Array3D<double> GetMomenta();

        //--------------------------------------------------------------------------
        // Methods for the optimizer
        //--------------------------------------------------------------------------

        virtual double ComputeFunctional(const Array3D<double>& momenta, double &dataMatching, double &regularity);
        virtual Array3D<double> ComputeGradient(const Array3D<double>& momenta);

        //--------------------------------------------------------------------------
        // Main algorithm method
        //--------------------------------------------------------------------------

        virtual Array3D<double> Run();

};

#endif // REGRESSIONVELOCITY_H
