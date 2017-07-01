#include <stdio.h>

#include "saveshapesandvectors.h"
#include "regressionacceleration.h"
#include "polydatawriter.h"

int SaveShapesAndVectors::_numSourceShapes = 0;
int* SaveShapesAndVectors::_numSourcePtsArray = 0;
int* SaveShapesAndVectors::_numSourceTrisArray = 0;
vector< Array2D<int> > SaveShapesAndVectors::_sourceTris;

void SaveShapesAndVectors::SetNumSourceShapes(int numSourceShapes)
{
    _numSourceShapes = numSourceShapes;
    _numSourcePtsArray = new int[numSourceShapes];
    _numSourceTrisArray = new int[numSourceShapes];
}

void SaveShapesAndVectors::SetNumSourcePtsArray(int* numSourcePtsArray)
{
    for (int i=0; i<_numSourceShapes; i++)
    {
        _numSourcePtsArray[i] = numSourcePtsArray[i];
    }
}

void SaveShapesAndVectors::SetNumSourceTrisArray(int* numSourceTrisArray)
{
    for (int i=0; i<_numSourceShapes; i++)
    {
        _numSourceTrisArray[i] = numSourceTrisArray[i];
    }
}

void SaveShapesAndVectors::SetSourceTris(vector< Array2D<int> > sourceTris)
{
    _sourceTris = sourceTris;
}

//----------------------------------------------------------------
// SaveShapesAndVectors
//----------------------------------------------------------------
// Inputs: 
//  
//   
// Outputs: 
//   return - 
//----------------------------------------------------------------
// Saves regression shapes and vectors
//----------------------------------------------------------------
bool SaveShapesAndVectors::SaveRegression(Regression* regression)
{	
    Array3D<double> X = ((RegressionAcceleration*) regression)->GetX();
    Array3D<double> dX = ((RegressionAcceleration*) regression)->GetVelocity();
    Array3D<double> accel = ((RegressionAcceleration*) regression)->GetAcceleration();

    // Get total number of shape points
    int nx = (((RegressionAcceleration*) regression)->GetX()).GetWidth();
    // Get number of time points
    int T = (((RegressionAcceleration*) regression)->GetX()).GetHeight();

    // Point offset for writing vectors
    int pointOffset = 0;

    char* outputDir = (char*)"/home/sci/jfishbau/projects/neuro_exp/fishbaugh/exoshape_test/regression/iters/";

    // Loop over the source shapes
    for (int s=0; s<_numSourceShapes; s++)
    {
        // Loop over the number of time points
        for (int t=0; t<T; t++)
        {
            char filename[500];
            sprintf(filename, "%stemp_iteration_time_%.3d.vtk", outputDir, t);
            VTKPolyDataWriter writer(filename);

            int numSourcePts = _numSourcePtsArray[s];
            Array2D<double> curPts(3,numSourcePts);
            Array2D<double> curDX(3,numSourcePts);
            Array2D<double> curAccel(3,numSourcePts);

            for (int i=0; i<3; i++)
            {
                for (int j=0; j<numSourcePts; j++)
                {
                    curPts(i,j) = X(i,pointOffset+j,t);
                    curDX(i,j) = dX(i,pointOffset+j,t);
                    curAccel(i,j) = accel(i,pointOffset+j,t);
                }
            }

            writer.AddField((char*)"velocity", curDX);
            writer.AddField((char*)"acceleration", curAccel);

            Array2D<int> curTris = _sourceTris[s];

            bool didWrite = writer.WritePointsAndTris(curPts, curTris);
        }
    }

    return true;
}

