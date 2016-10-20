//--------------------------------------------------------------------------
// EXVTKPolyDataWriter: A simple and lightweight VTKPolyData writer 
// which writes legacy ascii vtk files
//--------------------------------------------------------------------------

#include "vtkpolydatawriter.h"

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EXVTKPolyDataWriter
//----------------------------------------------------------------
// Inputs: 
//   filename - full path of VTK file to write
//   
// Outputs: 
//----------------------------------------------------------------
// Init constructor - Recommended constructor
//----------------------------------------------------------------
VTKPolyDataWriter::VTKPolyDataWriter(char* filename)
{
    this->_vtkFile = new ofstream(filename);
}

//----------------------------------------------------------------
// ~EXVTKPolyDataWriter
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
VTKPolyDataWriter::~VTKPolyDataWriter()
{
    delete this->_vtkFile;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// ADD FIELD
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// AddField
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//----------------------------------------------------------------
// Writes a legacy ascii vtk file
//----------------------------------------------------------------
void VTKPolyDataWriter::AddField(char* fieldName, const Array2D<double>& field)
{
    _fieldNames.push_back(fieldName);
    _fields.push_back(field);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// WRITE VTK FILE
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// WritePointsAndTris
//----------------------------------------------------------------
// Inputs: 
//   pts - 2D array (dim, num_pts) of surface vertices
//   tris - 2D array (dim, num_tris) of triangles
//
// Outputs: 
//   return - true if successful, false otherwise
//----------------------------------------------------------------
// Writes a legacy ascii vtk file
//----------------------------------------------------------------
bool VTKPolyDataWriter::WritePointsAndTris(const Array2D<double>& pts, const Array2D<int>& tris)
{
    int numPts = pts.GetWidth();
    int numTris = tris.GetWidth();
    char buffer[100];

    if (!this->_vtkFile->is_open())
    {
        this->_vtkFile->close();
        return false;
    }

    // Write the header information
    *(this->_vtkFile)<<"# vtk DataFile Version 3.0\n";
    *(this->_vtkFile)<<"vtk output\n";
    *(this->_vtkFile)<<"ASCII\n";
    *(this->_vtkFile)<<"DATASET POLYDATA\n";
    sprintf (buffer, "POINTS %d double\n", numPts);
    *(this->_vtkFile)<<buffer;

    // Write the points
    for (int i=0; i<numPts; i++)
    {
        *(this->_vtkFile)<<pts(0,i)<<" "<<pts(1,i)<<" "<<pts(2,i)<<"\n";
    }

    memset(buffer, 0, 100*sizeof(char));
    sprintf(buffer, "POLYGONS %d %d\n", numTris, numTris*4);
    *(this->_vtkFile)<<buffer;

    // Write the triangles
    for (int i=0; i<numTris; i++)
    {
        *(this->_vtkFile)<<"3 "<<tris(0,i)<<" "<<tris(1,i)<<" "<<tris(2,i)<<"\n";
    }

    //memset(buffer, 0, 100*sizeof(char));
    //sprintf(buffer, "CELL_DATA %d\n", numTris);
    //*(this->_vtkFile)<<buffer;
    //memset(buffer, 0, 100*sizeof(char));
    //sprintf(buffer, "POINT_DATA %d\n", numPts);
    //*(this->_vtkFile)<<buffer;

    memset(buffer, 0, 100*sizeof(char));
    sprintf(buffer, "POINT_DATA %d\n", numPts);
    *(this->_vtkFile)<<buffer;

    memset(buffer, 0, 100*sizeof(char));
    sprintf(buffer, "FIELD FieldData %d\n", (int)_fieldNames.size());
    *(this->_vtkFile)<<buffer;

    // Write the fields
    for (unsigned int i=0; i<_fields.size(); i++)
    {
        memset(buffer, 0, 100*sizeof(char));
        sprintf(buffer, "%s %d %d double\n", _fieldNames[i], 3, numPts);
        *(this->_vtkFile)<<buffer;

        // Get the current field
        Array2D<double> curField = this->_fields[i];

        // Write the vectors
        for (int i=0; i<numPts; i++)
        {
            *(this->_vtkFile)<<curField(0,i)<<" "<<curField(1,i)<<" "<<curField(2,i)<<"\n";
        }
    }

    this->_vtkFile->close();

    return true;
}

