#include "shape4dstate.h"
#include "polydatawriter.h"

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// STATIC MEMBER VARIABLES
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

char*Shape4DState:: _expName = (char*)"";				// The name of the experiment
char* Shape4DState::_driverXML = (char*)"";				// Path to the xml driver file for the experiment

int Shape4DState::_saveProgressIters = -1;				// Save progress every n iterations

Array3D<double> Shape4DState::_X;						// The source points
vector<char *> Shape4DState::_vectorNames;				// A vector of names for vector information
vector< Array3D<double> > Shape4DState::_vectors;		// A vector of EX2DArrays of vector information

char* Shape4DState::_outputDir = (char*)"";				// The directory for output
vector<char *> Shape4DState::_outputNames;				// A vector of output names for regression shapes

int Shape4DState::_T = 0;								// The number of time points
int Shape4DState::_numSourceShapes = 0;					// The number of source shapes
int* Shape4DState::_numSourcePtsArray = 0;				// The number of points for each source shape
int* Shape4DState::_numSourceTrisArray = 0;				// The number of tris for each source shape
vector< Array2D<int> > Shape4DState::_sourceTris;		// A vector of EX2DArrays of triangles for each source shape

int Shape4DState::_iteration = 0;						// The number of iterations
vector<double> Shape4DState::_dataMatching;				// An array of data matching values
vector<double> Shape4DState::_regularity;				// An array of regularity values
vector<double> Shape4DState::_stepsize;					// An array of stepsize values

bool Shape4DState::_save = false;						// Should we save anything during optimization?
bool Shape4DState::_writeShapesAndVectors = false;		// Is it time to save shapes and vectors?

bool Shape4DState::_hasConverged = false;				// Has the algorithm converged?

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// SETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// SetExperimentName
//----------------------------------------------------------------
// Inputs: 
//   expName - the name for the experiment
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the name of the experiment, for saving purposes.  This
// determines the name of the .exo progress file.
//----------------------------------------------------------------
void Shape4DState::SetExperimentName(char* expName)
{
    _expName = expName;
}

//----------------------------------------------------------------
// SetXMLDriver
//----------------------------------------------------------------
// Inputs: 
//   driverXML - the full path to the XML driver file
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the the full path of the XML driver
//----------------------------------------------------------------
void Shape4DState::SetXMLDriver(char* driverXML)
{
    _driverXML = driverXML;
}

//----------------------------------------------------------------
// SetSaveProgressEveryNIterations
//----------------------------------------------------------------
// Inputs: 
//   saveProgressIters - the number of iterations between saves
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the number of iterations between saves
//----------------------------------------------------------------
void Shape4DState::SetSaveProgressEveryNIterations(int saveProgressIters)
{
    _saveProgressIters = saveProgressIters;
}

//----------------------------------------------------------------
// SetX
//----------------------------------------------------------------
// Inputs: 
//   X - a 3D array of source points over time
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the source points over time
//----------------------------------------------------------------
void Shape4DState::SetX(const Array3D<double> & X)
{
    _X = X;
}

//----------------------------------------------------------------
// SetVectors
//----------------------------------------------------------------
// Inputs: 
//   names - a vector of names for corresponding vector fields
//   vectors - a vector of vector fields (velocity, accel, ...)
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the vector fields (and their names) for saving
//----------------------------------------------------------------
void Shape4DState::SetVectors(vector<char *> names, vector< Array3D<double> > vectors)
{
    _vectorNames = names;
    _vectors = vectors;
}

//----------------------------------------------------------------
// SetOutputDir
//----------------------------------------------------------------
// Inputs: 
//   outputDir - full path for output
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the full path of the directory to store output
//----------------------------------------------------------------
void Shape4DState::SetOutputDir(char* outputDir)
{
    _outputDir = outputDir;
}

//----------------------------------------------------------------
// SetOutputNames
//----------------------------------------------------------------
// Inputs: 
//   outputNames - a vector of prefixes for each saved shape
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the prefixes to use for each shape to save
//----------------------------------------------------------------
void Shape4DState::SetOutputNames(vector<char *> outputNames)
{
    _outputNames = outputNames;
}

//----------------------------------------------------------------
// SetT
//----------------------------------------------------------------
// Inputs: 
//   T - the number of time points
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the number of time points in the experiment
//----------------------------------------------------------------
void Shape4DState::SetT(int T)
{
    _T = T;
}

//----------------------------------------------------------------
// SetNumSourceShapes
//----------------------------------------------------------------
// Inputs: 
//   numSourceShapes - the number of source shapes used
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the number of source shapes used in the experiment
//----------------------------------------------------------------
void Shape4DState::SetNumSourceShapes(int numSourceShapes)
{
    _numSourceShapes = numSourceShapes;
    _numSourcePtsArray = new int[numSourceShapes];
    _numSourceTrisArray = new int[numSourceShapes];
}

//----------------------------------------------------------------
// SetNumSourcePtsArray
//----------------------------------------------------------------
// Inputs: 
//   numSourcePtsArray - an array where each element denotes the
//                       number of points corresponding to that
//                       source shape
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the number of points for each source shape
//----------------------------------------------------------------
void Shape4DState::SetNumSourcePtsArray(int* numSourcePtsArray)
{
    for (int i=0; i<_numSourceShapes; i++)
    {
        _numSourcePtsArray[i] = numSourcePtsArray[i];
    }
}

//----------------------------------------------------------------
// SetNumSourceTrisArray
//----------------------------------------------------------------
// Inputs: 
//   numSourceTrisArray - an array where each element denotes the
//                       number of tris corresponding to that
//                       source shape
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the number of triangles for each source shape
//----------------------------------------------------------------
void Shape4DState::SetNumSourceTrisArray(int* numSourceTrisArray)
{
    for (int i=0; i<_numSourceShapes; i++)
    {
        _numSourceTrisArray[i] = numSourceTrisArray[i];
    }
}

//----------------------------------------------------------------
// SetSourceTris
//----------------------------------------------------------------
// Inputs: 
//   sourceTris - a vector of triangle connectivity for each 
//				  source shape
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the connectivity (triangles) for each source shape
//----------------------------------------------------------------
void Shape4DState::SetSourceTris(vector< Array2D<int> > sourceTris)
{
    _sourceTris = sourceTris;
}

//----------------------------------------------------------------
// SetDataMatchingValues
//----------------------------------------------------------------
// Inputs: 
//   dataMatchings - a vector of data matching values obtained 
//					 from optimzation
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the data matching values obtained from optimization
//----------------------------------------------------------------
void Shape4DState::SetDataMatchingValues(vector<double> dataMatchings)
{
    _dataMatching = dataMatchings;
}

//----------------------------------------------------------------
// SetRegularityValues
//----------------------------------------------------------------
// Inputs: 
//   regularities - a vector of regularity values obtained 
//					from optimzation
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the regularity values obtained from optimization
//----------------------------------------------------------------
void Shape4DState::SetRegularityValues(vector<double> regularities)
{
    _regularity = regularities;
}

//----------------------------------------------------------------
// SetStepsizeValues
//----------------------------------------------------------------
// Inputs: 
//   stepsizes - a vector of step sizes obtained from optimzation
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the step sizes obtained from optimization
//----------------------------------------------------------------
void Shape4DState::SetStepsizeValues(vector<double> stepsizes)
{
    _stepsize = stepsizes;
}

//----------------------------------------------------------------
// SetIterations
//----------------------------------------------------------------
// Inputs: 
//   iterations - the number of finished optimization iterations
//
// Outputs: 
//    
//----------------------------------------------------------------
// Sets the number of completed optimization iterations
//----------------------------------------------------------------
void Shape4DState::SetIterations(int iterations)
{
    _iteration = iterations;
}

//----------------------------------------------------------------
// SetShouldSave
//----------------------------------------------------------------
// Inputs: 
//   yesNo - true if progress should be saved, false otherwise 
//
// Outputs: 
//    
//----------------------------------------------------------------
// Should the progress of the experiment be saved?
//----------------------------------------------------------------
void Shape4DState::SetShouldSave(bool yesNo)
{
    _save = yesNo;
}

//----------------------------------------------------------------
// SetShouldWriteShapesAndVectors
//----------------------------------------------------------------
// Inputs: 
//   yesNo - true temporary shapes/vectors should be saved, 
//           false otherwise 
//
// Outputs: 
//    
//----------------------------------------------------------------
// Should temporary shapes/vectors of the experiment be saved?
//----------------------------------------------------------------
void Shape4DState::SetShouldWriteShapesAndVectors(bool yesNo)
{
    _writeShapesAndVectors = yesNo;
}

//----------------------------------------------------------------
// SetHasConverged
//----------------------------------------------------------------
// Inputs: 
//   yesNo - true if optimization has converged, false otherwise 
//
// Outputs: 
//    
//----------------------------------------------------------------
// Has the optimization converged?
//----------------------------------------------------------------
void Shape4DState::SetHasConverged(bool yesNo)
{
    _hasConverged = yesNo;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// GETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// GetIteration
//----------------------------------------------------------------
// Inputs: 
//
// Outputs: 
//    return - the number of completed iterations
//----------------------------------------------------------------
// Returns the number of completed iterations from optimization
//----------------------------------------------------------------
int Shape4DState::GetIteration()
{
    return _iteration;
}

//----------------------------------------------------------------
// GetCurrentStepsize
//----------------------------------------------------------------
// Inputs: 
//
// Outputs: 
//    return - the current step size from optimization
//----------------------------------------------------------------
// Returns the current step size from optimization
//----------------------------------------------------------------
double Shape4DState::GetCurrentStepsize()
{
    if (_stepsize.size() == 0)
    {
        return 0.0;
    }
    else
    {
        return _stepsize.back();
    }
}	

//----------------------------------------------------------------
// GetSaveProgressEveryNIterations
//----------------------------------------------------------------
// Inputs: 
//
// Outputs: 
//    return - the number of iterations between saving progress
//----------------------------------------------------------------
// Returns the number of iterations between saving progress
//----------------------------------------------------------------
int Shape4DState::GetSaveProgressEveryNIterations()
{
    return _saveProgressIters;
}

//----------------------------------------------------------------
// GetDataMatchingValues
//----------------------------------------------------------------
// Inputs: 
//
// Outputs: 
//    return - a vector of data matching values
//----------------------------------------------------------------
// Returns the data matching values
//----------------------------------------------------------------
vector<double> Shape4DState::GetDataMatchingValues()
{
    return _dataMatching;
}

//----------------------------------------------------------------
// GetRegularityValues
//----------------------------------------------------------------
// Inputs: 
//
// Outputs: 
//    return - a vector of regularity values
//----------------------------------------------------------------
// Returns the regularity values
//----------------------------------------------------------------
vector<double> Shape4DState::GetRegularityValues()
{
    return _regularity;
}

//----------------------------------------------------------------
// GetStepsizeValues
//----------------------------------------------------------------
// Inputs: 
//
// Outputs: 
//    return - a vector of step sizes
//----------------------------------------------------------------
// Returns the step sizes
//----------------------------------------------------------------
vector<double> Shape4DState::GetStepsizeValues()
{
    return _stepsize;
}

//----------------------------------------------------------------
// GetShouldSave
//----------------------------------------------------------------
// Inputs: 
//
// Outputs: 
//    return - true if temporary progress should be saved, false
//             otherwise
//----------------------------------------------------------------
// Should the progress of the experiment be saved?
//----------------------------------------------------------------
bool Shape4DState::GetShouldSave()
{
    return _save;
}

//----------------------------------------------------------------
// GetShouldWriteShapesAndVectors
//----------------------------------------------------------------
// Inputs: 
//
// Outputs: 
//    return - true if temporary shapes/vectors should be saved, 
//			   false otherwise
//----------------------------------------------------------------
// Should temporary shapes/vectors of the experiment be saved?
//----------------------------------------------------------------
bool Shape4DState::GetShouldWriteShapesAndVectors()
{
    return _writeShapesAndVectors;
}

//----------------------------------------------------------------
// GetHasConverged
//----------------------------------------------------------------
// Inputs: 
//
// Outputs: 
//    return - true if optimization has converged, false 
//             otherwise
//----------------------------------------------------------------
// Has optimization converged?
//----------------------------------------------------------------
bool Shape4DState::GetHasConverged()
{
    return _hasConverged;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// UPDATERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// ResetOptimizationState
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//
//----------------------------------------------------------------
// Clears the functional values and sets iterations to zero
//----------------------------------------------------------------
void Shape4DState::ResetOptimizationState()
{
    _iteration = 0;
    _dataMatching.clear();
    _regularity.clear();
}

//----------------------------------------------------------------
// UpdateIteration
//----------------------------------------------------------------
// Inputs: 
//   iterationNum - the iteration number
//   dataMatching - the data matching value
//   regularity - the regularity value
//   stepsize - the step size of the optimization
//
// Outputs: 
//    
//----------------------------------------------------------------
// Update number of iterations, data matching, regularity, and
// stepsize values.  Adds each to a vector storing the values.
//----------------------------------------------------------------
void Shape4DState::UpdateIteration(int iterationNum, double dataMatching, double regularity, double stepsize)
{
    _iteration = iterationNum;
    _dataMatching.push_back(dataMatching);
    _regularity.push_back(regularity);
    _stepsize.push_back(stepsize);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// METHODS FOR SAVING
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// SaveStateAccel
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//   return - true if successful, false otherwise
//----------------------------------------------------------------
// Saves the state of the regression (acceleration)
//----------------------------------------------------------------
bool Shape4DState::SaveStateAccel()
{	
    // Set the filename
    char filename[500];
    sprintf(filename, "%s%s.exo", _outputDir, _expName);

    // Open the file and set precision
    ofstream exoshapeStateFile(filename);
    exoshapeStateFile.precision(20);
    exoshapeStateFile<<"ExoshapeAccel Version 1.0\n";

    // Write path to xml driver
    exoshapeStateFile<<"EXPERIMENT FILE\n";
    exoshapeStateFile<<_driverXML<<"\n";

    // Write current iteration progress
    exoshapeStateFile<<"COMPLETED ITERATIONS\n";
    exoshapeStateFile<<_iteration<<"\n";

    // Write data matching, regularity, and stepsize valuesls
    exoshapeStateFile<<"MATCHING REGULARITY STEPSIZE\n";
    for (int i=0; i<_iteration; i++)
    {
        exoshapeStateFile<<_dataMatching[i]<<" "<<_regularity[i]<<" "<<_stepsize[i]<<"\n";
        char buffer[100];
        sprintf (buffer, "%10.10lf %10.10lf %10.10lf", _dataMatching[i], _regularity[i], _stepsize[i]);
    }

    // Write init velocity
    exoshapeStateFile<<"INIT VELOCITY\n";
    for (unsigned int v=0; v<_vectors.size(); v++)
    {
        if (strcmp(_vectorNames[v],"velocity") == 0)
        {
            Array3D<double> curVec = _vectors[v];
            exoshapeStateFile<<curVec.GetLength()<<" "<<curVec.GetWidth()<<"\n";

            for (int j=0; j<_X.GetWidth(); j++)
            {
                for (int k=0; k<curVec.GetLength(); k++)
                {
                    if (k == curVec.GetLength()-1)
                    {
                        exoshapeStateFile<<curVec(k,j,0)<<"\n";
                    }
                    else
                    {
                        exoshapeStateFile<<curVec(k,j,0)<<" ";
                    }
                }
            }
        }
    }

    // Write impulse
    exoshapeStateFile<<"IMPULSE\n";
    for (unsigned int v=0; v<_vectors.size(); v++)
    {
        if (strcmp(_vectorNames[v],"impulse") == 0)
        {
            Array3D<double> curVec = _vectors[v];
            exoshapeStateFile<<curVec.GetLength()<<" "<<curVec.GetWidth()<<" "<<curVec.GetHeight()<<"\n";

            for (int i=0; i<curVec.GetHeight(); i++)
            {
                for (int j=0; j<_X.GetWidth(); j++)
                {
                    for (int k=0; k<curVec.GetLength(); k++)
                    {
                        if (k == curVec.GetLength()-1)
                        {
                            exoshapeStateFile<<curVec(k,j,i)<<"\n";
                        }
                        else
                        {
                            exoshapeStateFile<<curVec(k,j,i)<<" ";
                        }
                    }
                }
            }
        }
    }

    // We are done, close the stream
    exoshapeStateFile.close();

    return true;
}

//----------------------------------------------------------------
// SaveStateVelocity
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//   return - true if successful, false otherwise
//----------------------------------------------------------------
// Saves the state of the regression (velocity)
//----------------------------------------------------------------
bool Shape4DState::SaveStateVelocity()
{	
    // Set the filename
    char filename[500];
    sprintf(filename, "%s%s.exo", _outputDir, _expName);

    // Open the file and set precision
    ofstream exoshapeStateFile(filename);
    exoshapeStateFile.precision(20);
    exoshapeStateFile<<"ExoshapeAccel Version 1.0\n";

    // Write path to xml driver
    exoshapeStateFile<<"EXPERIMENT FILE\n";
    exoshapeStateFile<<_driverXML<<"\n";

    // Write current iteration progress
    exoshapeStateFile<<"COMPLETED ITERATIONS\n";
    exoshapeStateFile<<_iteration<<"\n";

    // Write data matching, regularity, and stepsize valuesls
    exoshapeStateFile<<"MATCHING REGULARITY STEPSIZE\n";
    for (int i=0; i<_iteration; i++)
    {
        exoshapeStateFile<<_dataMatching[i]<<" "<<_regularity[i]<<" "<<_stepsize[i]<<"\n";
        char buffer[100];
        sprintf (buffer, "%10.10lf %10.10lf %10.10lf", _dataMatching[i], _regularity[i], _stepsize[i]);
    }

    // Write momenta
    exoshapeStateFile<<"MOMENTA\n";
    for (unsigned int v=0; v<_vectors.size(); v++)
    {
        if (strcmp(_vectorNames[v],"momenta") == 0)
        {
            Array3D<double> curVec = _vectors[v];
            exoshapeStateFile<<curVec.GetLength()<<" "<<curVec.GetWidth()<<" "<<curVec.GetHeight()<<"\n";

            for (int i=0; i<curVec.GetHeight(); i++)
            {
                for (int j=0; j<_X.GetWidth(); j++)
                {
                    for (int k=0; k<curVec.GetLength(); k++)
                    {
                        if (k == curVec.GetLength()-1)
                        {
                            exoshapeStateFile<<curVec(k,j,i)<<"\n";
                        }
                        else
                        {
                            exoshapeStateFile<<curVec(k,j,i)<<" ";
                        }
                    }
                }
            }
        }
    }

    // We are done, close the stream
    exoshapeStateFile.close();

    return true;
}

//----------------------------------------------------------------
// SaveShapesAndVectors
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//   return - true if successful, false otherwise 
//----------------------------------------------------------------
// Saves regression shapes and vectors
//----------------------------------------------------------------
bool Shape4DState::SaveShapesAndVectors()
{
    // Did writing go successfully?
    bool didWrite = false;

    // Get number of time points
    int T = _X.GetHeight();

    // Point offset for writing vectors
    int pointOffset = 0;

    // Loop over the source shapes
    for (int s=0; s<_numSourceShapes; s++)
    {
        // Loop over the number of time points
        for (int t=0; t<T; t++)
        {
            char filename[500];

            // Is this the final save?
            if (_hasConverged)
            {
                sprintf(filename, "%s%sfinal_time_%.3d.vtk", _outputDir, _outputNames[s], t);
            }
            // Else save the temporary files
            else
            {
                sprintf(filename, "%s%siter_%.3d_time_%.3d.vtk", _outputDir, _outputNames[s], Shape4DState::GetIteration()+1, t);
            }
            VTKPolyDataWriter writer(filename);

            int numSourcePts = _numSourcePtsArray[s];
            Array2D<double> curPts(3,numSourcePts);

            // Select only the points belonging to the current shape and time point
            for (int i=0; i<3; i++)
            {
                for (int j=0; j<numSourcePts; j++)
                {
                    curPts(i,j) = _X(i,pointOffset+j,t);
                }
            }

            // Loop over all the vector fields
            for (unsigned int v=0; v<_vectors.size(); v++)
            {
                // The name of this vector field
                char* vectorName = _vectorNames[v];

                // Get all the current vectors at this time point
                Array3D<double> allVectors(3,numSourcePts, T);
                allVectors = _vectors[v];

                // Hold only the vectors that belong to the current shape and time point
                Array2D<double> curVector(3,numSourcePts);

                // Loop over all the points
                for (int i=0; i<3; i++)
                {
                    for (int j=0; j<numSourcePts; j++)
                    {

                        // Select only the vectors that belong to the current shape and time point
                        curVector(i,j) = allVectors(i,pointOffset+j,t);
                    }

                }

                // Add the vector field to the writer
                writer.AddField(vectorName, curVector);

            }

            // Get the triangles for the source
            Array2D<int> curTris = _sourceTris[s];

            didWrite = writer.WritePointsAndTris(curPts, curTris);
        }

        pointOffset += _numSourcePtsArray[s];
    }

    return didWrite;
}
