//--------------------------------------------------------------------------
// RunExperiment:  Driver for running shape regression experiments
//--------------------------------------------------------------------------

#ifndef RUNEXPERIMENT_H
#define RUNEXPERIMENT_H

#include <vector>
#include "array3d.h"

using namespace std;

class RunExperiment
{

    private:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------

        bool _continueExp;						// Are we continuing a previous experiment?
        bool _skipRegressionVelocity;           // When continuing, we skip the RegressionVelocity if it has already finished
        Array2D<double> *_initV0;				// The initial velocity for resuming experiments
        Array3D<double> *_impulse;              // The full time impulse for resuming experiments
        Array3D<double> *_momenta;              // The full time momenta for resuming experiments

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        RunExperiment();
        ~RunExperiment();

        //--------------------------------------------------------------------------
        // Experiment drivers
        //--------------------------------------------------------------------------

        void StartExperiment(char* pathToFile);
        void ContinueExperiment(char* pathToFile);

};

#endif // RUNEXPERIMENT_H
