//--------------------------------------------------------------------------
// EXGrid: A grid used for fast kernel computations
//--------------------------------------------------------------------------

#include "grid.h"

#include <string.h>			// for memcpy	
#include <math.h>			  // For ceil
#include <stdlib.h>			// For abs
#include <fftw3.h>			// For FFT

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EXGrid
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
Grid::Grid()
{
	this->_ratio = 0.2f;
	this->_pas = 0.0f;
	memset(this->_origin, 0.0f, 3*sizeof(double));
	memset(this->_long, 0, 3*sizeof(int));
	this->_fft3kD = 0;
}

//----------------------------------------------------------------
// EXGrid
//----------------------------------------------------------------
// Inputs: 
//   grid - EXGrid object to copy
//
// Outputs:
//----------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------
Grid::Grid(const Grid& grid)
{
	this->_ratio = grid._ratio;
	this->_pas = grid._pas;
	memcpy(this->_origin, grid._origin, 3*sizeof(double));
	memcpy(this->_long, grid._long, 3*sizeof(double));
  if (this->_fft3kD != 0)
	{
		delete [] this->_fft3kD;
	}
	this->_fft3kD = new double[((this->_long[0]/2)+1)*this->_long[1]*this->_long[2]];
	memcpy(this->_fft3kD, grid._fft3kD, (((this->_long[0]/2)+1)*this->_long[1]*this->_long[2])*sizeof(double));
}

//----------------------------------------------------------------
// EXGrid
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
Grid::~Grid()
{
	delete [] this->_fft3kD;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// GETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// GetPas
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//   return - value for grid accuracy
//----------------------------------------------------------------
// Returns the value for grid accuracy
//----------------------------------------------------------------
double Grid::GetPas()
{
	return this->_pas;
}

//----------------------------------------------------------------
// GetOrigin
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//   return - array[3] of the grid origin
//----------------------------------------------------------------
// Returns the grid origin
//----------------------------------------------------------------
double* Grid::GetOrigin()
{
	return this->_origin;
}

//----------------------------------------------------------------
// GetLong
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//   return - array[3] of the grid dimensions
//----------------------------------------------------------------
// Returns the grid dimensions
//----------------------------------------------------------------
int* Grid::GetLong()
{
	return this->_long;
}

//----------------------------------------------------------------
// GetFFT3KD
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//   return - array of FFT of the grid positions
//----------------------------------------------------------------
// Returns the FFT of the grid positions
//----------------------------------------------------------------
double* Grid::GetFFT3KD()
{
	return this->_fft3kD;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// HELPER FUNCTIONS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// ComputeFFT3kD
//----------------------------------------------------------------
// Inputs: 
//   sigmaV - standard deviation of Gaussian kernel
//
// Outputs:
//   return - array of FFT of the grid positions
//----------------------------------------------------------------
// Computes the the FFT of the grid positions
//----------------------------------------------------------------
void Grid::ComputeFFT3kD(double sigmaV)
{
	fftw_complex *k3D = (fftw_complex *) fftw_malloc(this->_long[0]*this->_long[1]*this->_long[2]*sizeof(fftw_complex));
	
	// Initialize with zeros
	for (int i=0; i<this->_long[0]; i++)
	{
		for (int j=0; j<this->_long[1]; j++)
		{
			for (int k=0; k<this->_long[2]; k++)
			{
				int index = i*this->_long[1]*this->_long[2] + j*this->_long[2] + k;
				k3D[index][0] = 0.0;
				k3D[index][1] = 0.0;
			}
		}
	}
	
	int demilongX = this->_long[0]/2;
	int demilongY = this->_long[1]/2;
	int demilongZ = this->_long[2]/2;
	double sigmaV2 = sigmaV*sigmaV;
	double pas2 = this->_pas*this->_pas;
	
	for (int i=0; i<demilongX+1; i++)
	{
		for (int j=0; j<demilongY+1; j++)
		{
			for (int k=0; k<demilongZ+1; k++)
			{
				int index = i*this->_long[1]*this->_long[2] + j*this->_long[2] + k;
				k3D[index][0] = exp(-pas2*(i*i+j*j+k*k)/sigmaV2);
			}
		}
	}

	int vx0[2] = {0, demilongX};
	int vxm[2] = {this->_long[0]-1, demilongX+1};
	int vx2[2] = {1, demilongX-1};
	
	int vy0[2] = {0, demilongY};
	int vym[2] = {this->_long[1]-1, demilongY+1};
	int vy2[2] = {1, demilongY-1};
	
	int vz0[2] = {0, demilongZ};
	int vzm[2] = {this->_long[2]-1, demilongZ+1};
	int vz2[2] = {1, demilongZ-1};
	
	// k(vxm,vy0,vz0) = k(vx2,vy0,vz0);/*
	for (int i=0; i<=abs(vxm[1]-vxm[0]); i++)
	{
    //printf("  ComputeFFT3kD: %d\n", i);
		for (int j=0; j<=abs(vy0[1]-vy0[0]); j++)
		{
      //printf("  ComputeFFT3kD: %d %d\n", i, j);
			for (int k=0; k<=abs(vz0[1]-vz0[0]); k++)
			{
        //printf("  ComputeFFT3kD: %d %d %d\n", i, j, kc);
				int l0 = vxm[0]-i;  int l1 = vy0[0]+j;  int l2 = vz0[0]+k;
				int r0 = vx2[0]+i;  int r1 = vy0[0]+j;  int r2 = vz0[0]+k;
				
				int indexl = l2 + this->_long[2] * (l1 + this->_long[1] * l0);
				int indexr = r2 + this->_long[2] * (r1 + this->_long[1] * r0);
				k3D[indexl][0] = k3D[indexr][0];
			}
		}
	}

	// k(vx0,vym,vz0) = k(vx0,vy2,vz0)
	for (int i=0; i<=abs(vx0[1]-vx0[0]); i++)
	{
		for (int j=0; j<=abs(vym[1]-vym[0]); j++)
		{
			for (int k=0; k<=abs(vz0[1]-vz0[0]); k++)
			{
				int l0 = vx0[0]+i;  int l1 = vym[0]-j;  int l2 = vz0[0]+k;
				int r0 = vx0[0]+i;  int r1 = vy2[0]+j;  int r2 = vz0[0]+k;
				
				int indexl = l2 + this->_long[2] * (l1 + this->_long[1] * l0);
				int indexr = r2 + this->_long[2] * (r1 + this->_long[1] * r0);
				k3D[indexl][0] = k3D[indexr][0];
			}
		}
	}

	// k(vx0,vy0,vzm) = k(vx0,vy0,vz2);
	for (int i=0; i<=abs(vx0[1]-vx0[0]); i++)
	{
		for (int j=0; j<=abs(vy0[1]-vy0[0]); j++)
		{
			for (int k=0; k<=abs(vzm[1]-vzm[0]); k++)
			{
				int l0 = vx0[0]+i;  int l1 = vy0[0]+j;  int l2 = vzm[0]-k;
				int r0 = vx0[0]+i;  int r1 = vy0[0]+j;  int r2 = vz2[0]+k;
				
				int indexl = l2 + this->_long[2] * (l1 + this->_long[1] * l0);
				int indexr = r2 + this->_long[2] * (r1 + this->_long[1] * r0);
				k3D[indexl][0] = k3D[indexr][0];
			}
		}
	}
	
	// k(vxm,vym,vz0) = k(vx2,vy2,vz0);
	for (int i=0; i<=abs(vxm[1]-vxm[0]); i++)
	{
		for (int j=0; j<=abs(vym[1]-vym[0]); j++)
		{
			for (int k=0; k<=abs(vz0[1]-vz0[0]); k++)
			{
				int l0 = vxm[0]-i;  int l1 = vym[0]-j;  int l2 = vz0[0]+k;
				int r0 = vx2[0]+i;  int r1 = vy2[0]+j;  int r2 = vz0[0]+k;
				
				int indexl = l2 + this->_long[2] * (l1 + this->_long[1] * l0);
				int indexr = r2 + this->_long[2] * (r1 + this->_long[1] * r0);
				k3D[indexl][0] = k3D[indexr][0];
			}
		}
	} 

	// k(vxm,vy0,vzm) = k(vx2,vy0,vz2);
	for (int i=0; i<=abs(vxm[1]-vxm[0]); i++)
	{
		for (int j=0; j<=abs(vy0[1]-vy0[0]); j++)
		{
			for (int k=0; k<=abs(vzm[1]-vzm[0]); k++)
			{
				int l0 = vxm[0]-i;  int l1 = vy0[0]+j;  int l2 = vzm[0]-k;
				int r0 = vx2[0]+i;  int r1 = vy0[0]+j;  int r2 = vz2[0]+k;
				
				int indexl = l2 + this->_long[2] * (l1 + this->_long[1] * l0);
				int indexr = r2 + this->_long[2] * (r1 + this->_long[1] * r0);
				k3D[indexl][0] = k3D[indexr][0];
			}
		}
	} 

	// k(vx0,vym,vzm) = k(vx0,vy2,vz2);
	for (int i=0; i<=abs(vx0[1]-vx0[0]); i++)
	{
		for (int j=0; j<=abs(vym[1]-vym[0]); j++)
		{
			for (int k=0; k<=abs(vzm[1]-vzm[0]); k++)
			{
				int l0 = vx0[0]+i;  int l1 = vym[0]-j;  int l2 = vzm[0]-k;
				int r0 = vx0[0]+i;  int r1 = vy2[0]+j;  int r2 = vz2[0]+k;
				
				int indexl = l2 + this->_long[2] * (l1 + this->_long[1] * l0);
				int indexr = r2 + this->_long[2] * (r1 + this->_long[1] * r0);
				k3D[indexl][0] = k3D[indexr][0];
			}
		}
	} 

	// k(vxm,vym,vzm) = k(vx2,vy2,vz2);
	for (int i=0; i<=abs(vxm[1]-vxm[0]); i++)
	{
		for (int j=0; j<=abs(vym[1]-vym[0]); j++)
		{
			for (int k=0; k<=abs(vzm[1]-vzm[0]); k++)
			{
				int l0 = vxm[0]-i;  int l1 = vym[0]-j;  int l2 = vzm[0]-k;
				int r0 = vx2[0]+i;  int r1 = vy2[0]+j;  int r2 = vz2[0]+k;
				
				int indexl = l2 + this->_long[2] * (l1 + this->_long[1] * l0);
				int indexr = r2 + this->_long[2] * (r1 + this->_long[1] * r0);
				k3D[indexl][0] = k3D[indexr][0];
			}
		}
	} 
	
	fftw_complex *fft3k = (fftw_complex *) fftw_malloc(this->_long[0]*this->_long[1]*this->_long[2]*sizeof(fftw_complex));
	fftw_plan foward = fftw_plan_dft_3d(this->_long[0], this->_long[1], this->_long[2], k3D, fft3k, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(foward);
	
  //if (this->_fft3kD != 0)
	//{
	//	delete [] this->_fft3kD;
	//}
	//this->_fft3kD = new double[(demilongX+1)*this->_long[1]*this->_long[2]];
	
	for (int i=0; i<demilongX+1; i++)
	{
		for (int j=0; j<this->_long[1]; j++)
		{
			for (int k=0; k<this->_long[2]; k++)
			{
				int index = i*this->_long[1]*this->_long[2] + j*this->_long[2] + k;
				this->_fft3kD[index] = fft3k[index][0];
			}
		}
	}
	
	// Clean up memory
	fftw_destroy_plan(foward);
  fftw_free(k3D);
	fftw_free(fft3k);
	fftw_cleanup();
	
	//return fft3kD;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// GRID METHODS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// SetGrid 
//----------------------------------------------------------------
// Inputs: 
//   mini - array[3] of minimum coordinates
//   maxi - array[3] of maximum coordinates
//   sigmaV - standard deviation of Gaussian kernel
//
// Outputs: 
//----------------------------------------------------------------
// Computes all properties of the grid
//----------------------------------------------------------------
void Grid::SetGrid(double mini[3], double maxi[3], double sigmaV)
{	
	this->_pas = this->_ratio * sigmaV;
	
	double newlong0, newlong1, newlong2;
	int newlong[3];
	
	newlong0 = (maxi[0]-mini[0])/this->_pas + 3.0f/this->_ratio;
	newlong[0] = (newlong0<=16)*16 + (newlong0>16)*((newlong0<=32)*32 + (newlong0>32)*((newlong0<=64)*64 + (newlong0>64)*2*ceil(newlong0/2.0f)));
	newlong1 = (maxi[1]-mini[1])/this->_pas + 3.0f/this->_ratio;
	newlong[1] = (newlong1<=16)*16 + (newlong1>16)*((newlong1<=32)*32 + (newlong1>32)*((newlong1<=64)*64 + (newlong1>64)*2*ceil(newlong1/2.0f)));
	newlong2 = (maxi[2]-mini[2])/this->_pas + 3.0f/this->_ratio;
	newlong[2] = (newlong2<=16)*16 + (newlong2>16)*((newlong2<=32)*32 + (newlong2>32)*((newlong2<=64)*64 + (newlong2>64)*2*ceil(newlong2/2.0f)));
		
	this->_origin[0] = mini[0] - (newlong[0]*this->_pas - maxi[0]+mini[0]) / 2.0f;
	this->_origin[1] = mini[1] - (newlong[1]*this->_pas - maxi[1]+mini[1]) / 2.0f;
	this->_origin[2] = mini[2] - (newlong[2]*this->_pas - maxi[2]+mini[2]) / 2.0f;

	//printf("originx %0.4f  originy %0.4f  originz %0.4f\n", this->_origin[0], this->_origin[1], this->_origin[2]);
	
	if ((newlong[0] != this->_long[0]) || (newlong[1] != this->_long[1]) || (newlong[2] != this->_long[2]))
	{
		this->_long[0] = newlong[0];
		this->_long[1] = newlong[1];
		this->_long[2] = newlong[2];
		//this->_fft3kD = this->ComputeFFT3kD(sigmaV);
    if (this->_fft3kD != 0)
		{
			delete [] this->_fft3kD;
		}
		this->_fft3kD = new double[((this->_long[0]/2)+1)*this->_long[1]*this->_long[2]];
    this->ComputeFFT3kD(sigmaV);
	}
}

//----------------------------------------------------------------
// ChangeGrid 
//----------------------------------------------------------------
// Inputs: 
//   mini - array[3] of minimum coordinates
//   maxi - array[3] of maximum coordinates
//   sigmaV - standard deviation of Gaussian kernel
//
// Outputs: 
//   return - a boolean specifying if the grid should change
//----------------------------------------------------------------
// Determines if the grid needs to be updated
//----------------------------------------------------------------
bool Grid::ChangeGrid(double mini[3], double maxi[3], double sigmaV)
{
	//printf("In change grid\n");

	//printf("minix %0.4f  miniy %0.4f  miniz %0.4f\n", mini[0], mini[1], mini[2]);
	//printf("maxix %0.4f  maxiy %0.4f  maxiz %0.4f\n", maxi[0], maxi[1], maxi[2]);
	//printf("SigmaV %0.4f\n", sigmaV);

	//printf("originx %0.4f  originy %0.4f  originz %0.4f\n", this->_origin[0], this->_origin[1], this->_origin[2]);
	
	//printf("Origin not defined??\n");

	double sum = 0.0f;
	for (int i=0; i<3; i++)
	{
		sum += (maxi[i] > (this->_origin[i] + this->_pas*this->_long[i] - sigmaV)) | (mini[i] < (this->_origin[i] + sigmaV));
	}
	bool p = (sum != 0);

	//printf("Done with sum\n");
	
	double newLong0;
	double newLong1;
	double newLong2;
	
	if (!p)
	{
		//printf("Inside not p\n");

		newLong0 = (maxi[0]-mini[0])/this->_pas + 3.0f/this->_ratio;
		newLong0 = (newLong0<=16)*16 + (newLong0>16)*((newLong0<=32)*32 + (newLong0>32)*((newLong0<=64)*64 + (newLong0>64)*2*ceil(newLong0/2)));
  	newLong1 = (maxi[1]-mini[1])/this->_pas + 3.0f/this->_ratio;
		newLong1 = (newLong1<=16)*16 + (newLong1>16)*((newLong1<=32)*32 + (newLong1>32)*((newLong1<=64)*64 + (newLong1>64)*2*ceil(newLong1/2)));
		newLong2 = (maxi[2]-mini[2])/this->_pas + 3.0f/this->_ratio;
		newLong2 = (newLong2<=16)*16 + (newLong2>16)*((newLong2<=32)*32 + (newLong2>32)*((newLong2<=64)*64 + (newLong2>64)*2*ceil(newLong2/2)));
		
		p = ((newLong0 < this->_long[0]) + (newLong1 < this->_long[1]) + (newLong2 < this->_long[2])) != 0;
	}

	//printf("Returning p\n");
	
	return p;
}

//----------------------------------------------------------------
// UpdateGrid 
//----------------------------------------------------------------
// Inputs: 
//   mini - array[3] of minimum coordinates
//   maxi - array[3] of maximum coordinates
//   sigmaV - standard deviation of Gaussian kernel
//
// Outputs: 
//----------------------------------------------------------------
// Updates the grid based on new min and max positions
//----------------------------------------------------------------
void Grid::UpdateGrid(double mini[3],double maxi[3], double sigmaV)
{
	//printf("In update grid\n");

	if (this->ChangeGrid(mini, maxi, sigmaV))
	{
		//printf("Need to change grid\n");
		this->SetGrid(mini, maxi, sigmaV);
	}
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// OVERLOADED OPERATORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// operator = 
//----------------------------------------------------------------
// Inputs: 
//   shape - EXGrid object to copy
//
// Outputs:
//   return - copied EXGrid object 
//----------------------------------------------------------------
// Overloaded assignment operator
//----------------------------------------------------------------
Grid& Grid::operator = (const Grid& grid)
{
	if (this != &grid)
	{
		this->_ratio = grid._ratio;
		this->_pas = grid._pas;
		memcpy(this->_origin, grid._origin, 3*sizeof(double));
		memcpy(this->_long, grid._long, 3*sizeof(double));
		if (this->_fft3kD != 0)
		{
			delete [] this->_fft3kD;
		}
		this->_fft3kD = new double[((this->_long[0]/2)+1)*this->_long[1]*this->_long[2]];
		memcpy(this->_fft3kD, grid._fft3kD, (((this->_long[0]/2)+1)*this->_long[1]*this->_long[2])*sizeof(double));
	}
	
	return *this;
}
