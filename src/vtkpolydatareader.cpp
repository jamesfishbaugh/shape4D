//--------------------------------------------------------------------------
// EXVTKPolyDataReader: A simple and lightweight VTKPolyData reader 
// which reads legacy ascii vtk files
//--------------------------------------------------------------------------

#include "vtkpolydatareader.h"
#include "helper.h"

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
	this->_vtkFile = new ifstream(filename);
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
	delete this->_vtkFile;
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
	
	return true;
}
