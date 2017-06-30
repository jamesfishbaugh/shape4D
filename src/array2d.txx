//--------------------------------------------------------------------------
// EX2DArray:  Templated 2D array class stored as a 1D array (row-major)
//--------------------------------------------------------------------------

#include "array2d.h"
#include <string.h>			// For memcpy

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EX2DArray
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
template <class datatype>
Array2D<datatype>::Array2D()
{
    this->_length = 0;
    this->_width = 0;
    this->_lengthTimesWidth = 0;

    this->_array =  0;
}

//----------------------------------------------------------------
// EX2DArray
//----------------------------------------------------------------
// Inputs: 
//   length - the length of the array
//   width - the width of the array
//
// Outputs: 
//----------------------------------------------------------------
// Init constructor
//----------------------------------------------------------------
template <class datatype>
Array2D<datatype>::Array2D(int length, int width)
{
    this->_length = length;
    this->_width = width;
    this->_lengthTimesWidth = length*width;

    this->_array = new datatype[this->_lengthTimesWidth];
}

//----------------------------------------------------------------
// EX2DArray
//----------------------------------------------------------------
// Inputs: 
//   array - EX2DArray object to copy
//
// Outputs:
//----------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------
template <class datatype>
Array2D<datatype>::Array2D(const Array2D& array)
{
    this->_length = array._length;
    this->_width = array._width;
    this->_lengthTimesWidth = array._lengthTimesWidth;

    this->_array = new datatype[this->_lengthTimesWidth];
    memcpy(this->_array, array._array, this->_lengthTimesWidth*sizeof(datatype));
}

//----------------------------------------------------------------
// ~EX2DArray
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
template <class datatype>
Array2D<datatype>::~Array2D()
{
    delete [] this->_array;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// GETTERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// GetLength
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//   return - length of the array
//----------------------------------------------------------------
// Returns the length of the array
//----------------------------------------------------------------
template <class datatype>
const int Array2D<datatype>::GetLength() const
{
    return this->_length;
}

//----------------------------------------------------------------
// GetWidth
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//   return - width of the array
//----------------------------------------------------------------
// Returns the width of the array
//----------------------------------------------------------------
template <class datatype>
const int Array2D<datatype>::GetWidth() const
{
    return this->_width;
}

//----------------------------------------------------------------
// GetArray
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//   return - the array
//----------------------------------------------------------------
// Returns the array
//----------------------------------------------------------------
template <class datatype>
datatype* Array2D<datatype>::GetArray()
{
    return this->_array;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// ACCESSORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// GetAt
//----------------------------------------------------------------
// Inputs: 
//   l - length index of desired value
//   w - width index of desired value
//
// Outputs:
//   return - value of the array at the given position
//----------------------------------------------------------------
// Returns the value of the array at the given position
//----------------------------------------------------------------
template <class datatype>
datatype Array2D<datatype>::GetAt(int l, int w)
{
    int index = l*this->_width + w;
    return this->_array[index];
}

//----------------------------------------------------------------
// SetAt
//----------------------------------------------------------------
// Inputs: 
//   l - length index of desired value
//   w - width index of desired value
//   value - the value to set
//
// Outputs:
//----------------------------------------------------------------
// Sets the value of the array at the given position
//----------------------------------------------------------------
template <class datatype>
void Array2D<datatype>::SetAt(int l, int w, datatype value)
{
    int index = l*this->_width + w;
    this->_array[index] = value;
}	

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// FILL METHODS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// FillArray
//----------------------------------------------------------------
// Inputs: 
//   value - value to fill the array with
//
// Outputs:
//----------------------------------------------------------------
// Fills the array at with the given value
//----------------------------------------------------------------
template <class datatype>
void Array2D<datatype>::FillArray(datatype value)
{
    memset(this->_array, value, this->_lengthTimesWidth*sizeof(datatype));
}

//----------------------------------------------------------------
// FillArray
//----------------------------------------------------------------
// Inputs: 
//   array - array to fill the array with
//
// Outputs:
//----------------------------------------------------------------
// Fills the array at with the given array
//----------------------------------------------------------------
template <class datatype>
void Array2D<datatype>::FillArray(datatype* array)
{
    memcpy(this->_array, array, this->_lengthTimesWidth*sizeof(datatype));
}

//----------------------------------------------------------------
// FillArray
//----------------------------------------------------------------
// Inputs: 
//   array - array to fill the array with
//
// Outputs:
//----------------------------------------------------------------
// Fills the array at with the given 2D array
//----------------------------------------------------------------
template <class datatype>
void Array2D<datatype>::FillArray(datatype** array)
{
    for (int i=0; i<this->_length; i++)
    {
        for (int j=0; j<this->_width; j++)
        {
            int index = i*this->_width + j;
            this->_array[index] = array[i][j];
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// ARRAY METHODS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// Sum
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//   return - sum of array elements
//----------------------------------------------------------------
// Returns the sum of the array elements
//----------------------------------------------------------------
template <class datatype>
datatype Array2D<datatype>::Sum()
{
    datatype sum = 0;

    for (int i=0; i<this->_lengthTimesWidth; i++)
    {
        sum += this->_array[i];
    }

    return sum;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// OVERLOADED OPERATORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// operator =
//----------------------------------------------------------------
// Inputs: 
//   array - EX2DArray object to copy
//
// Outputs:
//   return - copied EX2DArray object
//----------------------------------------------------------------
// Overloaded assignment operator
//----------------------------------------------------------------
template <class datatype>
Array2D<datatype>&  Array2D<datatype>::operator = (const Array2D& array)
{
    if (this != &array)
    {
        this->_length = array._length;
        this->_width = array._width;
        this->_lengthTimesWidth = array._lengthTimesWidth;

        if (this->_array != 0)
        {
            delete [] this->_array;
        }
        this->_array = new datatype[this->_lengthTimesWidth];
        memcpy(this->_array, array._array, this->_lengthTimesWidth*sizeof(datatype));
    }

    return *this;
}

//----------------------------------------------------------------
// operator ==
//----------------------------------------------------------------
// Inputs: 
//   array - array to compare with
//
// Outputs:
//   return - true if equal, false otherwise
//----------------------------------------------------------------
// Overloaded comparison operator
//----------------------------------------------------------------
template <class datatype>
bool Array2D<datatype>::operator == (const Array2D& array)
{
    return (this == &array);
}

//----------------------------------------------------------------
// operator !=
//----------------------------------------------------------------
// Inputs: 
//   array - array to compare with
//
// Outputs:
//   return - true if not equal, false otherwise
//----------------------------------------------------------------
// Overloaded comparison operator
//----------------------------------------------------------------
template <class datatype>
bool Array2D<datatype>::operator != (const Array2D& array)
{
    return (this != &array);
}


//----------------------------------------------------------------
// operator ()
//----------------------------------------------------------------
// Inputs: 
//   l - length index of desired value
//   w - width index of desired value
//
// Outputs:
//   return - reference to value at given index
//----------------------------------------------------------------
// Overloaded () operator which allows access and modifcation
// of the array
//----------------------------------------------------------------
template <class datatype>
datatype& Array2D<datatype>::operator()(int l, int w)
{
    int index = l*this->_width + w;
    return this->_array[index];
}

//----------------------------------------------------------------
// operator ()
//----------------------------------------------------------------
// Inputs: 
//   l - length index of desired value
//   w - width index of desired value
//
// Outputs:
//   return - value at given index
//----------------------------------------------------------------
// Overloaded () operator which only allows access to the array
//----------------------------------------------------------------
template <class datatype>
const datatype& Array2D<datatype>::operator()(int l, int w) const
{
    int index = l*this->_width + w;
    return this->_array[index];
}

//----------------------------------------------------------------
// operator +
//----------------------------------------------------------------
// Inputs: 
//   value - value to add to the array
//
// Outputs:
//   return - array with elements added by value
//----------------------------------------------------------------
// Overloaded + operator with a scalar
//----------------------------------------------------------------
template <class datatype>
Array2D<datatype>  Array2D<datatype>::operator + (datatype value)
{
    Array2D<datatype> output(this->_length, this->_width);

    for (int i=0; i<this->_lengthTimesWidth; i++)
    {
        output._array[i] = this->_array[i] + value;
    }

    return output;
}

//----------------------------------------------------------------
// operator +
//----------------------------------------------------------------
// Inputs: 
//   array - array to add to the array
//
// Outputs:
//   return - array with elements added by elements of array
//----------------------------------------------------------------
// Overloaded + operator which implements elementwise addition
// between two arrays
//----------------------------------------------------------------
template <class datatype>
Array2D<datatype> Array2D<datatype>::operator + (const Array2D& array)
{
    Array2D<datatype> output(this->_length, this->_width);

    for (int i=0; i<this->_lengthTimesWidth; i++)
    {
        output._array[i] = this->_array[i] + array._array[i];
    }

    return output;
}

//----------------------------------------------------------------
// operator -
//----------------------------------------------------------------
// Inputs: 
//   value - value to subtract from the array
//
// Outputs:
//   return - array with elements subtracted by value
//----------------------------------------------------------------
// Overloaded - operator with a scalar
//----------------------------------------------------------------
template <class datatype>
Array2D<datatype>  Array2D<datatype>::operator - (datatype value)
{
    Array2D<datatype> output(this->_length, this->_width);

    for (int i=0; i<this->_lengthTimesWidth; i++)
    {
        output._array[i] = this->_array[i] - value;
    }

    return output;
}

//----------------------------------------------------------------
// operator -
//----------------------------------------------------------------
// Inputs: 
//   array - array to subtract from the array
//
// Outputs:
//   return - array with elements subtracted by elements of 
//            array
//----------------------------------------------------------------
// Overloaded - operator which implements elementwise 
// subtraction between two arrays
//----------------------------------------------------------------
template <class datatype>
Array2D<datatype> Array2D<datatype>::operator - (const Array2D& array) const
{
    Array2D<datatype> output(this->_length, this->_width);

    for (int i=0; i<this->_lengthTimesWidth; i++)
    {
        output._array[i] = this->_array[i] - array._array[i];
    }

    return output;
}

//----------------------------------------------------------------
// operator *
//----------------------------------------------------------------
// Inputs: 
//   value - value to multiply to the array
//
// Outputs:
//   return - array with elements multiplied by value
//----------------------------------------------------------------
// Overloaded * operator with a scalar
//----------------------------------------------------------------
template <class datatype>
Array2D<datatype>  Array2D<datatype>::operator * (datatype value) const
{
    Array2D<datatype> output(this->_length, this->_width);

    for (int i=0; i<this->_lengthTimesWidth; i++)
    {
        output._array[i] = this->_array[i] * value;
    }

    return output;
}

//----------------------------------------------------------------
// operator *
//----------------------------------------------------------------
// Inputs: 
//   array - array to multiply the array with
//
// Outputs:
//   return - array with elements multiplied by elements of 
//            array
//----------------------------------------------------------------
// Overloaded * operator which implements elementwise 
// multiplication between two arrays
//---------------------------------------------------------------
template <class datatype>
Array2D<datatype> Array2D<datatype>::operator * (const Array2D& array)
{
    Array2D<datatype> output(this->_length, this->_width);

    for (int i=0; i<this->_lengthTimesWidth; i++)
    {
        output._array[i] = this->_array[i] * array._array[i];
    }

    return output;
}

//----------------------------------------------------------------
// operator /
//----------------------------------------------------------------
// Inputs: 
//   value - value to divide the array by
//
// Outputs:
//   return - array with elements divided by value
//----------------------------------------------------------------
// Overloaded / operator with a scalar
//----------------------------------------------------------------
template <class datatype>
Array2D<datatype>  Array2D<datatype>::operator / (datatype value)
{
    Array2D<datatype> output(this->_length, this->_width);

    for (int i=0; i<this->_lengthTimesWidth; i++)
    {
        output._array[i] = this->_array[i] / value;
    }

    return output;
}

//----------------------------------------------------------------
// operator /
//----------------------------------------------------------------
// Inputs: 
//   array - array to divide the array by
//
// Outputs:
//   return - array with elements divided by elements of array
//----------------------------------------------------------------
// Overloaded * operator which implements elementwise division
// between two arrays
//---------------------------------------------------------------
template <class datatype>
Array2D<datatype> Array2D<datatype>::operator / (const Array2D& array)
{
    Array2D<datatype> output(this->_length, this->_width);

    for (int i=0; i<this->_lengthTimesWidth; i++)
    {
        output._array[i] = this->_array[i] / array._array[i];
    }

    return output;
}
