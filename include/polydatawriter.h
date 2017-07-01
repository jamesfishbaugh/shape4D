//--------------------------------------------------------------------------
// VTKPolyDataWriter: A simple and lightweight VTKPolyData writer
// which writes legacy ascii vtk files
//--------------------------------------------------------------------------

#ifndef VTKPOLYDATAWRITER_H
#define VTKPOLYDATAWRITER_H

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <string>

#include "array2d.h"
#include "array3d.h"

using namespace std;

class VTKPolyDataWriter
{

    private:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------
#ifdef USE_VTK
        string m_FileName;
#else
        ofstream* _vtkFile;		// For writing to a file
#endif

        vector<char*> _fieldNames;
        vector< Array2D<double> > _fields;

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        VTKPolyDataWriter(char* filename);
        ~VTKPolyDataWriter();

        //--------------------------------------------------------------------------
        // Add fields
        //--------------------------------------------------------------------------
        void AddField(char* fieldName, const Array2D<double>& field);

        //--------------------------------------------------------------------------
        // Write VTK file
        //--------------------------------------------------------------------------

        bool WritePointsAndTris(const Array2D<double>& pts, const Array2D<int>& tris);

};

#endif // VTKPOLYDATAWRITER_H
