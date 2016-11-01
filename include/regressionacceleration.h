//--------------------------------------------------------------------------
// RegressionAcceleration:  Shape regression parameterized by
// acceleration.
//--------------------------------------------------------------------------

#ifndef REGRESSIONACCELERATION_H
#define REGRESSIONACCELERATION_H

#include "regression.h"

class RegressionAcceleration : public Regression
{

    private:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------

        Array2D<double> _initV0;		// Initial initial velocity
        Array3D<double> _dX;		    // Velocity vectors
        Array3D<double> _accel;         // Acceleration vectors
        Array3D<double> _impulse;       // Final impulse vectors
        Array2D<double> _sobolevGrad;   // For smooth baseline shape estimation

        int _iteration;
        bool _continueRegression;

        //--------------------------------------------------------------------------
        // Write methods
        //--------------------------------------------------------------------------
        void WriteSelf();

        //--------------------------------------------------------------------------
        // Initialization methods
        //--------------------------------------------------------------------------

        void InitInitV0();
        void InitVelocity();

        //--------------------------------------------------------------------------
        // Computation of regression functional
        //--------------------------------------------------------------------------

        double ComputeRegularity(const Array3D<double>& impulse);
        double ComputeFunctional(const Array3D<double>& impulse, const Array2D<double> X0, const Array2D<double> V0, double &dataMatching, double &regularity);
        void ComputeGradient(const Array3D<double>& impulse, const Array2D<double>& X0, const Array2D<double>& V0,
                             Array3D<double>& gradImpulse, Array2D<double>& gradX0, Array2D<double>& gradV0);

        //--------------------------------------------------------------------------
        // Helper methods
        //--------------------------------------------------------------------------

        void SplitImpulseAndV0(const Array3D<double>& impulseAndV0, Array3D<double>&impulse, Array2D<double>& v0);
        void SplitGradient(const Array3D<double>& G, Array3D<double>& gradImpulse, Array2D<double>& gradX0, Array2D<double>& gradV0);
        void ComputeTrajectories(const Array3D<double>& impulse, const Array2D<double>& v0);
        void ComputeTrajectories(const Array3D<double>& impulse, const Array2D<double>& x0, const Array2D<double>& v0);

        //--------------------------------------------------------------------------
        // Own optimization methods
        //--------------------------------------------------------------------------
        void FISTA(Array3D<double>& impulse, Array2D<double>& X0, Array2D<double>& V0);
        void GradientDescent(Array3D<double>& impulse, Array2D<double>& X0, Array2D<double>& V0);
        void GradientDescentStep(Array3D<double>& impulseTest, Array2D<double>& X0Test, Array2D<double>& V0Test,
                                 const Array3D<double>& impulse, const Array2D<double>& X0, const Array2D<double>& V0,
                                 const Array3D<double>& gradImpulse, const Array2D<double>& gradX0, const Array2D<double>& gradV0,
                                 double stepImpulse, double stepX0andV0);

        double QDiffTerm(const Array3D<double>& impulseTest, const Array2D<double> X0Test, const Array2D<double> V0Test,
                         const Array3D<double>& impulse, const Array2D<double> X0, const Array2D<double> V0,
                         const Array3D<double>& gradImpulse, const Array2D<double>& gradX0, const Array2D<double>& gradV0,
                         double stepImpulse, double stepX0andV0);
    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        RegressionAcceleration();
        RegressionAcceleration(const RegressionParams& source);
        RegressionAcceleration(const RegressionParams& source, int numTargets, TargetData** targets);
        virtual ~RegressionAcceleration();

        //--------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------

        void SetImpulse(Array3D<double>& impulse);

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        Array3D<double> GetX();
        Array3D<double> GetVelocity();
        Array3D<double> GetAcceleration();
        Array3D<double> GetImpulse();

        //--------------------------------------------------------------------------
        // Methods for the optimizer
        //--------------------------------------------------------------------------

        virtual double ComputeFunctional(const Array3D<double>& impulseAndV0, double &dataMatching, double &regularity);
        virtual Array3D<double> ComputeGradient(const Array3D<double>& impulseAndV0);

        //--------------------------------------------------------------------------
        // Main algorithm method
        //--------------------------------------------------------------------------

        virtual Array3D<double> Run();

};

#endif // REGRESSIONACCELERATION_H
