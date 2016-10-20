//--------------------------------------------------------------------------
// EXHelper: Static class of helper functions for exoshape
//--------------------------------------------------------------------------

#include "helper.h"

float** Helper::Allocate2DFloatArray(int rows, int cols)
{
    float** array = new float*[rows];
    for (int i=0; i<rows; i++)
    {
        array[i] = new float[cols];
    }
    return array;
}

float*** Helper::Allocate3DFloatArray(int rows, int cols, int time)
{
    float*** array = new float**[rows];
    for (int i=0; i<rows; i++)
    {
        array[i] = new float*[cols];
        for (int j=0; j<cols; j++)
        {
            array[i][j] = new float[time];
        }
    }
    return array;
}

int** Helper::Allocate2DIntArray(int rows, int cols)
{
    int** array = new int*[rows];
    for (int i=0; i<rows; i++)
    {
        array[i] = new int[cols];
    }
    return array;
}

int*** Helper::Allocate3DIntArray(int rows, int cols, int time)
{
    int*** array = new int**[rows];
    for (int i=0; i<rows; i++)
    {
        array[i] = new int*[cols];
        for (int j=0; j<cols; j++)
        {
            array[i][j] = new int[time];
        }
    }
    return array;
}

void Helper::Fill2DFloatArray(float value, int rows, int cols, float** &array)
{
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            array[i][j] = value;
        }
    }
}

void Helper::Fill3DFloatArray(float value, int rows, int cols, int time, float*** &array)
{
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            for (int k=0; k<time; k++)
            {
                array[i][j][k] = value;
            }
        }
    }
}

void Helper::Fill2DIntArray(int value, int rows, int cols, int** &array)
{
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            array[i][j] = value;
        }
    }
}

void Helper::Fill3DIntArray(int value, int rows, int cols, int time, int*** &array)
{
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            for (int k=0; k<time; k++)
            {
                array[i][j][k] = value;
            }
        }
    }
}


void Helper::Copy2DFloatArray(int rows, int cols, float** arrayFrom, float** &arrayTo)
{
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            arrayTo[i][j] = arrayFrom[i][j];
        }
    }
}

void Helper::Copy3DFloatArray(int rows, int cols, int time, float*** arrayFrom, float*** &arrayTo)
{
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            for (int k=0; k<time; k++)
            {
                arrayTo[i][j][k] = arrayFrom[i][j][k];
            }
        }
    }
}

void Helper::Copy2DIntArray(int rows, int cols, int** arrayFrom, int** &arrayTo)
{
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            arrayTo[i][j] = arrayFrom[i][j];
        }
    }
}

void Helper::Copy3DIntArray(int rows, int cols, int time, int*** arrayFrom, int*** &arrayTo)
{
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            for (int k=0; k<time; k++)
            {
                arrayTo[i][j][k] = arrayFrom[i][j][k];
            }
        }
    }
}

float Helper::Sum2DFloatArray(int rows, int cols, float** array)
{
    float sum = 0.0f;

    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            sum += array[i][j];
        }
    }

    return sum;
}

float Helper::Sum3DFloatArray(int rows, int cols, int time, float*** array)
{
    float sum = 0.0f;

    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            for (int k=0; k<time; k++)
            {
                sum += array[i][j][k];
            }
        }
    }

    return sum;
}

int Helper::Sum2DIntArray(int rows, int cols, int** array)
{
    int sum = 0.0f;

    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            sum += array[i][j];
        }
    }

    return sum;
}

int Helper::Sum3DIntArray(int rows, int cols, int time, int*** array)
{
    int sum = 0.0f;

    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            for (int k=0; k<time; k++)
            {
                sum += array[i][j][k];
            }
        }
    }

    return sum;
}

void Helper::Deallocate2DFloatArray(float** array)
{
    delete [] *array;
    delete [] array;
}

void Helper::Deallocate3DFloatArray(float*** array)
{
    delete [] **array;
    delete [] *array;
    delete [] array;
}

void Helper::Deallocate2DIntArray(int** array)
{
    delete [] *array;
    delete [] array;
}

void Helper::Deallocate3DIntArray(int*** array)
{
    delete [] **array;
    delete [] *array;
    delete [] array;
}

