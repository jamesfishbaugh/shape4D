//--------------------------------------------------------------------------
// EX3DArray:  Templated 3D array class stored as a 1D array (row-major)
//--------------------------------------------------------------------------

#include "array3d.h"
#include <string.h>			// For memcpy

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EX3DArray
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
template <class datatype>
Array3D<datatype>::Array3D()
{
    this->_length = 0;
    this->_width = 0;
    this->_height = 0;
    this->_lengthTimesWidth = 0;
    this->_widthTimesHeight = 0;
    this->_lengthTimesWidthTimesHeight = 0;

    this->_array = 0;
}

//----------------------------------------------------------------
// EX3DArray
//----------------------------------------------------------------
// Inputs: 
//   length - the length of the array
//   width - the width of the array
//   height - the height of the array
//
// Outputs: 
//----------------------------------------------------------------
// Init constructor
//----------------------------------------------------------------
template <class datatype>
Array3D<datatype>::Array3D(int length, int width, int height)
{
    this->_length = length;
    this->_width = width;
    this->_height = height;
    this->_lengthTimesWidth = length*width;
    this->_widthTimesHeight = width*height;
    this->_lengthTimesWidthTimesHeight = length*width*height;

    this->_array = new datatype[length*width*height];
}

//----------------------------------------------------------------
// EX3DArray
//----------------------------------------------------------------
// Inputs: 
//   array - EX3DArray object to copy
//
// Outputs:
//----------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------
template <class datatype>
Array3D<datatype>::Array3D(const Array3D& array)
{
    this->_length = array._length;
    this->_width = array._width;
    this->_height = array._height;
    this->_lengthTimesWidth = array._lengthTimesWidth;
    this->_widthTimesHeight = array._widthTimesHeight;
    this->_lengthTimesWidthTimesHeight = array._lengthTimesWidthTimesHeight;

    this->_array = new datatype[this->_lengthTimesWidthTimesHeight];
    memcpy(this->_array, array._array, this->_lengthTimesWidthTimesHeight*sizeof(datatype));
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
Array3D<datatype>::~Array3D()
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
int Array3D<datatype>::GetLength() const
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
int Array3D<datatype>::GetWidth() const
{
    return this->_width;
}

//----------------------------------------------------------------
// GetHeight
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//   return - height of the array
//----------------------------------------------------------------
// Returns the height of the array
//----------------------------------------------------------------
template <class datatype>
int Array3D<datatype>::GetHeight() const
{
    return this->_height;
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
datatype* Array3D<datatype>::GetArray()
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
//   h - height index of desired value
//
// Outputs:
//   return - value of the array at the given position
//----------------------------------------------------------------
// Returns the value of the array at the given position
//----------------------------------------------------------------
template <class datatype>
datatype Array3D<datatype>::GetAt(int l, int w, int h) const
{
    int index = l*this->_widthTimesHeight + w*this->_height + h;
    return this->_array[index];
}

//----------------------------------------------------------------
// SetAt
//----------------------------------------------------------------
// Inputs: 
//   l - length index of desired value
//   w - width index of desired value
//   h - height index of desired value
//   value - the value to set
//
// Outputs:
//----------------------------------------------------------------
// Sets the value of the array at the given position
//----------------------------------------------------------------
template <class datatype>
void Array3D<datatype>::SetAt(int l, int w, int h, datatype value)
{
    int index = l*this->_widthTimesHeight + w*this->_height + h;
    this->_array[index] = value;
}	

//----------------------------------------------------------------
// Get2DSliceAtHeight
//----------------------------------------------------------------
// Inputs: 
//   h - height index to slice at
//
// Outputs:
//   return - 2D array (length, width)
//----------------------------------------------------------------
// Returns a 2D array at the given height index
//----------------------------------------------------------------
template <class datatype>
Array2D<datatype> Array3D<datatype>::Get2DSliceAtHeight(int h) const
{
    Array2D<datatype> output(this->_length, this->_width);

    for (int i=0; i<this->_length; i++)
    {
        for (int j=0; j<this->_width; j++)
        {
            output(i,j) = this->GetAt(i,j,h);
        }
    }

    return output;
}

//----------------------------------------------------------------
// Set2DSliceAtHeight
//----------------------------------------------------------------
// Inputs: 
//   array - 2D array (length, width) to set
//   h - height index to slice at
//
// Outputs:
//----------------------------------------------------------------
// Sets a 2D array at the given height index
//----------------------------------------------------------------
template <class datatype>
void Array3D<datatype>::Set2DSliceAtHeight(const Array2D<datatype>& array, int h)
{
    for (int i=0; i<this->_length; i++)
    {
        for (int j=0; j<this->_width; j++)
        {
            this->SetAt(i,j,h,array(i,j));
        }
    }
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
void Array3D<datatype>::FillArray(datatype value)
{
    memset(this->_array, value, this->_lengthTimesWidthTimesHeight*sizeof(datatype));
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
void Array3D<datatype>::FillArray(datatype* array)
{

    memcpy(this->_array, array, this->_lengthTimesWidthTimesHeight*sizeof(datatype));
}

//----------------------------------------------------------------
// FillArray
//----------------------------------------------------------------
// Inputs: 
//   array - array to fill the array with
//
// Outputs:
//----------------------------------------------------------------
// Fills the array at with the given 3D array
//----------------------------------------------------------------
template <class datatype>
void Array3D<datatype>::FillArray(datatype*** array)
{
    for (int l=0; l<this->_length; l++)
    {
        for (int w=0; w<this->_width; w++)
        {
            for (int h=0; h<this->_height; h++)
            {
                int index = l*this->_widthTimesHeight + w*this->_height + h;
                this->_array[index] = array[l][w][h];
            }
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
datatype Array3D<datatype>::Sum()
{
    datatype sum = 0;

    for (int i=0; i<this->_lengthTimesWidthTimesHeight; i++)
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
//   array - EX3DArray object to copy
//
// Outputs:
//   return - copied EX3DArray object
//----------------------------------------------------------------
// Overloaded assignment operator
//----------------------------------------------------------------
template <class datatype>
Array3D<datatype>&  Array3D<datatype>::operator = (const Array3D& array)
{
    if (this != &array)
    {
        this->_length = array._length;
        this->_width = array._width;
        this->_height = array._height;
        this->_lengthTimesWidth = array._lengthTimesWidth;
        this->_widthTimesHeight = array._widthTimesHeight;
        this->_lengthTimesWidthTimesHeight = array._lengthTimesWidthTimesHeight;

        if (this->_array != 0)
        {
            delete [] this->_array;
        }
        this->_array = new datatype[this->_lengthTimesWidthTimesHeight];
        memcpy(this->_array, array._array, this->_lengthTimesWidthTimesHeight*sizeof(datatype));
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
bool Array3D<datatype>::operator == (const Array3D& array)
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
bool Array3D<datatype>::operator != (const Array3D& array)
{
    return (this != &array);
}

//----------------------------------------------------------------
// operator ()
//----------------------------------------------------------------
// Inputs: 
//   l - length index of desired value
//   w - width index of desired value
//   h - height index of desired value
//
// Outputs:
//   return - reference to value at given index
//----------------------------------------------------------------
// Overloaded () operator which allows access and modifcation
// of the array
//----------------------------------------------------------------
template <class datatype>
datatype& Array3D<datatype>::operator()(int l, int w, int h)
{
    int index = l*this->_widthTimesHeight + w*this->_height + h;
    return this->_array[index];
}

//----------------------------------------------------------------
// operator ()
//----------------------------------------------------------------
// Inputs: 
//   l - length index of desired value
//   w - width index of desired value
//   h - height index of desired value
//
// Outputs:
//   return - value at given index
//----------------------------------------------------------------
// Overloaded () operator which only allows access to the array
//----------------------------------------------------------------
template <class datatype>
const datatype& Array3D<datatype>::operator()(int l, int w, int h) const
{
    int index = l*this->_widthTimesHeight + w*this->_height + h;
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
Array3D<datatype> Array3D<datatype>::operator + (datatype value)
{
    Array3D<datatype> output(this->_length, this->_width, this->_height);

    for (int i=0; i<this->_lengthTimesWidthTimesHeight; i++)
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
Array3D<datatype> Array3D<datatype>::operator + (const Array3D& array)
{
    Array3D<datatype> output(this->_length, this->_width, this->_height);

    for (int i=0; i<this->_lengthTimesWidthTimesHeight; i++)
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
Array3D<datatype> Array3D<datatype>::operator - (datatype value)
{
    Array3D<datatype> output(this->_length, this->_width, this->_height);

    for (int i=0; i<this->_lengthTimesWidthTimesHeight; i++)
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
Array3D<datatype> Array3D<datatype>::operator - (const Array3D& array) const
{
    Array3D<datatype> output(this->_length, this->_width, this->_height);

    for (int i=0; i<this->_lengthTimesWidthTimesHeight; i++)
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
Array3D<datatype> Array3D<datatype>::operator * (datatype value) const
{
    Array3D<datatype> output(this->_length, this->_width, this->_height);

    for (int i=0; i<this->_lengthTimesWidthTimesHeight; i++)
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
Array3D<datatype> Array3D<datatype>::operator * (const Array3D& array)
{
    Array3D<datatype> output(this->_length, this->_width, this->_height);

    for (int i=0; i<this->_lengthTimesWidthTimesHeight; i++)
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
Array3D<datatype> Array3D<datatype>::operator / (datatype value)
{
    Array3D<datatype> output(this->_length, this->_width, this->_height);

    for (int i=0; i<this->_lengthTimesWidthTimesHeight; i++)
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
Array3D<datatype> Array3D<datatype>::operator / (const Array3D& array)
{
    Array3D<datatype> output(this->_length, this->_width, this->_height);

    for (int i=0; i<this->_lengthTimesWidthTimesHeight; i++)
    {
        output._array[i] = this->_array[i] / array._array[i];
    }

    return output;
}
