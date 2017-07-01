//--------------------------------------------------------------------------
// VTKPolyDataReader: A simple and lightweight VTKPolyData reader
// which reads legacy ascii vtk files
//--------------------------------------------------------------------------

#ifndef VTKPOLYDATAREADER_H
#define VTKPOLYDATAREADER_H

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>

#include "array2d.h"

using namespace std;

class VTKPolyDataReader
{

    private:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------
#ifdef USE_VTK
        string m_FileName;
#else
        ifstream* _vtkFile;			// For reading from a file


        //--------------------------------------------------------------------------
        // Helper functions
        //--------------------------------------------------------------------------

        int GotoKeyword(string keyword);
#endif
    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        VTKPolyDataReader(char* filename);
        ~VTKPolyDataReader();

        //--------------------------------------------------------------------------
        // Read VTK file
        //--------------------------------------------------------------------------

        bool ReadPointsAndTris(Array2D<double> &pts, Array2D<int> &tris);

};

#endif // VTKPOLYDATAREADER_H
