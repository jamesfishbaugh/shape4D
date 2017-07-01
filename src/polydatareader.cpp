//--------------------------------------------------------------------------
// EXVTKPolyDataReader: A simple and lightweight VTKPolyData reader 
// which reads legacy ascii vtk files
//--------------------------------------------------------------------------

#include "polydatareader.h"
#include "helper.h"

#ifdef USE_VTK
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EXVTKPolyDataReader
//----------------------------------------------------------------
// Inputs: 
//   filename - full path of VTK file to read
//   
// Outputs: 
//----------------------------------------------------------------
// Init constructor - Recommended constructor
//----------------------------------------------------------------
VTKPolyDataReader::VTKPolyDataReader(char* filename)
{
#ifdef USE_VTK
m_FileName = filename ;
#else
    this->_vtkFile = new ifstream(filename);
#endif
}

//----------------------------------------------------------------
// ~EXVTKPolyDataReader
//----------------------------------------------------------------
// Inputs: 
//
// Outputs: 
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
VTKPolyDataReader::~VTKPolyDataReader()
{
#ifndef USE_VTK
    delete this->_vtkFile;
#endif
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// HELPER FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// GotoKeyword
//----------------------------------------------------------------
// Inputs: 
//   keyword - string to find in the file
//
// Outputs: 
//   return - absolute position of the keyword in the file or
//            -1 if the keyword is not found
//----------------------------------------------------------------
// Returns the absolute position of the keyword in the file or
// -1 if the keyword is not found
//----------------------------------------------------------------
#ifndef USE_VTK
int VTKPolyDataReader::GotoKeyword(string keyword)
{
    int kewordLocation = -1;
    string line;
    // Seek to the beginning so we search the entire file
    this->_vtkFile->seekg(0, ios::beg);

    // Look for the keyword
    while (!this->_vtkFile->eof())
    {
        int tempLocation = this->_vtkFile->tellg();
        getline(*(_vtkFile), line);

        int found = (int) line.find(keyword);

        if (found !=-1)
        {
            kewordLocation = tempLocation;
            break;
        }

    }

    return kewordLocation;
}
#endif
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// READ VTK FILE
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// ReadPointsAndTris
//----------------------------------------------------------------
// Inputs: 
//   pts - 2D array (dim, num_pts) of surface vertices for
//         return by reference
//   tris - 2D array (dim, num_tris) of triangles for return by
//          reference
//
// Outputs: 
//   return - true if successful, false otherwise
//----------------------------------------------------------------
// Read a legacy ascii vtk file
//----------------------------------------------------------------
bool VTKPolyDataReader::ReadPointsAndTris(Array2D<double> &pts, Array2D<int> &tris)
{
#ifdef USE_VTK
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(m_FileName.c_str());
    reader->Update();
    vtkPolyData* polydata = reader->GetOutput();
    if(!polydata)
    {
      return false;
    }
    vtkPoints* readerPoints = polydata->GetPoints();
    Array2D<double> tmpPts(3, polydata->GetNumberOfPoints() );
    double lpts[3];
    for(unsigned int j=0; j < polydata->GetNumberOfPoints() ; j++)
    {
        readerPoints->GetPoint(j, lpts);
        pts(0,j) = lpts[0];
        pts(1,j) = lpts[1];
        pts(2,j) = lpts[2];
    }
    Array2D<int> tmpTris(3,polydata->GetNumberOfCells());
    vtkCell* cell;
    vtkIdList* ptsId;
    for(unsigned int j=0; j < polydata->GetNumberOfCells() ; j++)
    {
        cell = polydata->GetCell(j);
        if(cell->GetNumberOfPoints() != 3)
        {
            printf("***ERROR*** All cells must be triangles.");
            exit(1);
        }
        ptsId = cell->GetPointIds();
        tris(0,j) = ptsId->GetId(0);
        tris(1,j) = ptsId->GetId(1);
        tris(2,j) = ptsId->GetId(2);
    }
    pts = tmpPts;
    tris = tmpTris;
#else
    string line;
    char* charLine;

    if (!this->_vtkFile->is_open())
    {
        this->_vtkFile->close();
        return false;
    }

    // Search for vtk point data
    int pointsLocation = GotoKeyword("POINTS");
    if (pointsLocation == -1)
    {
        this->_vtkFile->close();
        return false;
    }

    // Read in point data
    this->_vtkFile->seekg(pointsLocation, ios::beg);
    getline(*(this->_vtkFile), line);
    charLine = new char[line.size()+1];
    charLine[line.size()] = 0;
    memcpy(charLine, line.c_str(), line.size());

    // The number of points
    int numPts;
    sscanf(charLine,"%*s %d %*s", &numPts);
    Array2D<double> tempPts(3,numPts);

    for (int i=0; i<numPts; i++)
    {
        double curValueX;
        *(this->_vtkFile) >> curValueX;
        double curValueY;
        *(this->_vtkFile) >> curValueY;
        double curValueZ;
        *(this->_vtkFile) >> curValueZ;

        tempPts(0,i) = curValueX;
        tempPts(1,i) = curValueY;
        tempPts(2,i) = curValueZ;
    }

    // Search for vtk polygon data
    int trisLocation = GotoKeyword("POLYGONS");
    if (trisLocation == -1)
    {
        this->_vtkFile->close();
        pts = tempPts;
        //EX2DArray<int> tempTris(3,numPts);
        //tris = tempTris;
        return true;
    }

    // Read in triangle data
    this->_vtkFile->seekg(trisLocation, ios::beg);
    getline(*(this->_vtkFile), line);
    delete [] charLine;
    charLine = new char[line.size()+1];
    charLine[line.size()] = 0;
    memcpy(charLine, line.c_str(), line.size());

    // The number of triangles
    int numTris;
    sscanf(charLine,"%*s %d %*d", &numTris);
    Array2D<int> tempTris(3,numTris);

    for (int i=0; i<numTris; i++)
    {
        int dim;
        *(this->_vtkFile) >> dim;
        int triX;
        *(this->_vtkFile) >> triX;
        double triY;
        *(this->_vtkFile) >> triY;
        double triZ;
        *(this->_vtkFile) >> triZ;

        tempTris(0,i) = triX;
        tempTris(1,i) = triY;
        tempTris(2,i) = triZ;
    }

    this->_vtkFile->close();

    // Make copies for return
    pts = tempPts;
    tris = tempTris;

    // Clear up memory
    delete [] charLine;
#endif
    return true;
}
