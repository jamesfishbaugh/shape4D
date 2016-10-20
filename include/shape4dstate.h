//--------------------------------------------------------------------------
// Shape4DState:  Save the current state of the regression.  This is
// a static class shared among all classes that need it.
//--------------------------------------------------------------------------

#ifndef SHAPE4DSTATE_H
#define SHAPE4DSTATE_H

#include <vector>

#include "array2d.h"
#include "array3d.h"

using namespace std;

class Shape4DState
{

    private:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------

        static char* _expName;								// The name of the experiment
        static char* _driverXML;							// Path to the xml driver file for the experiment

        static int _saveProgressIters;						// Save progress every n iterations

        static Array3D<double> _X;                          // The source points
        static vector<char *> _vectorNames;                 // A vector of names for vector information
        static vector< Array3D<double> > _vectors;          // A vector of EX3DArrays of vector information

        static char* _outputDir;							// The directory for output
        static vector<char *> _outputNames;                 // A vector of output names for regression shapes

        static int _T;										// The number of time points
        static int _numSourceShapes;						// The number of source shapes
        static int* _numSourcePtsArray;                     // The number of points for each source shape
        static int* _numSourceTrisArray;					// The number of tris for each source shape
        static vector< Array2D<int> > _sourceTris;          // A vector of EX2DArrays of triangles for each source shape

        static int _iteration;								// The number of iterations
        static vector<double> _dataMatching;				// An array of data matching values
        static vector<double> _regularity;                  // An array of regularity values
        static vector<double> _stepsize;					// An array of stepsize values

        static bool _save;									// Should we save anything during optimization?
        static bool _writeShapesAndVectors;                 // Is it time to save shapes and vectors?

        static bool _hasConverged;							// Has the algorithm converged?

    public:

        //--------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------

        static void SetExperimentName(char* expName);
        static void SetXMLDriver(char* driverXML);

        static void SetSaveProgressEveryNIterations(int saveProgressIters);

        static void SetX(const Array3D<double> & X);
        static void SetVectors(vector<char *> names, vector< Array3D<double> > vectors);

        static void SetOutputDir(char* outputDir);
        static void SetOutputNames(vector<char *> outputNames);

        static void SetT(int T);
        static void SetNumSourceShapes(int numSourceShapes);
        static void SetNumSourcePtsArray(int* numSourcePtsArray);
        static void SetNumSourceTrisArray(int* numSourceTrisArray);
        static void SetSourceTris(vector< Array2D<int> > sourceTris);

        static void SetDataMatchingValues(vector<double> dataMatchings);
        static void SetRegularityValues(vector<double> regularities);
        static void SetStepsizeValues(vector<double> stepsizes);

        static void SetIterations(int iterations);

        static void SetShouldSave(bool yesNo);
        static void SetShouldWriteShapesAndVectors(bool yesNo);

        static void SetHasConverged(bool yesNo);

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        static int GetIteration();
        static double GetCurrentStepsize();

        static int GetSaveProgressEveryNIterations();

        static vector<double> GetDataMatchingValues();
        static vector<double> GetRegularityValues();
        static vector<double> GetStepsizeValues();

        static bool GetShouldSave();
        static bool GetShouldWriteShapesAndVectors();

        static bool GetHasConverged();

        //--------------------------------------------------------------------------
        // Updaters
        //--------------------------------------------------------------------------

        static void ResetOptimizationState();
        static void UpdateIteration(int iterationNum, double dataMatching, double regularity, double stepsize);

        //--------------------------------------------------------------------------
        // Methods for saving
        //--------------------------------------------------------------------------

        static bool SaveStateAccel();
        static bool SaveStateVelocity();
        static bool SaveShapesAndVectors();

};

#endif // SHAPE4DSTATE_H
