//--------------------------------------------------------------------------
// EX1DArray:  Templated 1D array class
//--------------------------------------------------------------------------

#include "array1d.h"
#include <string.h>			// For memcpy

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EX1DArray
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
template <class datatype>
Array1D<datatype>::Array1D()
{
    this->_length = 0;

    this->_array =  0;
}

//----------------------------------------------------------------
// EX1DArray
//----------------------------------------------------------------
// Inputs: 
//   length - the size of the array
//
// Outputs: 
//----------------------------------------------------------------
// Init constructor
//----------------------------------------------------------------
template <class datatype>
Array1D<datatype>::Array1D(int length)
{
    this->_length = length;

    this->_array = new datatype[length];
}

//----------------------------------------------------------------
// EX1DArray
//----------------------------------------------------------------
// Inputs: 
//   array - EX1DArray object to copy
//
// Outputs:
//----------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------
template <class datatype>
Array1D<datatype>::Array1D(const Array1D& array)
{
    this->_length = array._length;

    this->_array = new datatype[this->_length];
    memcpy(this->_array, array._array, this->_length*sizeof(datatype));
}

//----------------------------------------------------------------
// ~EX1DArray
//----------------------------------------------------------------
// Inputs: 
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
template <class datatype>
Array1D<datatype>::~Array1D()
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
//   return - size of the array
//----------------------------------------------------------------
// Returns the size of the array
//----------------------------------------------------------------
template <class datatype>
const int Array1D<datatype>::GetLength() const
{
    return this->_length;
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
datatype* Array1D<datatype>::GetArray()
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
//   l - index of desired value
//
// Outputs:
//   return - value of the array at the given position
//----------------------------------------------------------------
// Returns the value of the array at the given position
//----------------------------------------------------------------
template <class datatype>
datatype Array1D<datatype>::GetAt(int l)
{
    return this->_array[l];
}

//----------------------------------------------------------------
// SetAt
//----------------------------------------------------------------
// Inputs: 
//   l - index of desired value
//   value - the value to set
//
// Outputs:
//----------------------------------------------------------------
// Sets the value of the array at the given position
//----------------------------------------------------------------
template <class datatype>
void Array1D<datatype>::SetAt(int l, datatype value)
{
    this->_array[l] = value;
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
void Array1D<datatype>::FillArray(datatype value)
{
    memset(this->_array, value, this->_length*sizeof(datatype));
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
void Array1D<datatype>::FillArray(datatype* array)
{
    memcpy(this->_array, array, this->_length*sizeof(datatype));
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
datatype Array1D<datatype>::Sum()
{
    datatype sum = 0;

    for (int i=0; i<this->_length; i++)
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
//   array - EX1DArray object to copy
//
// Outputs:
//   return - copied EX1DArray object
//----------------------------------------------------------------
// Overloaded assignment operator
//----------------------------------------------------------------
template <class datatype>
Array1D<datatype>&  Array1D<datatype>::operator = (const Array1D& array)
{
    if (this != &array)
    {
        this->_length = array._length;

        if (this->_array != 0)
        {
            delete [] this->_array;
        }
        this->_array = new datatype[this->_length];
        memcpy(this->_array, array._array, this->_length*sizeof(datatype));
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
bool Array1D<datatype>::operator == (const Array1D& array)
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
bool Array1D<datatype>::operator != (const Array1D& array)
{
    return (this != &array);
}

//----------------------------------------------------------------
// operator ()
//----------------------------------------------------------------
// Inputs: 
//   l - index to access
//
// Outputs:
//   return - reference to value at given index
//----------------------------------------------------------------
// Overloaded () operator which allows access and modifcation
// of the array
//----------------------------------------------------------------
template <class datatype>
datatype& Array1D<datatype>::operator()(int l)
{
    return this->_array[l];
}

//----------------------------------------------------------------
// operator ()
//----------------------------------------------------------------
// Inputs: 
//   l - index to access
//
// Outputs:
//   return - value at given index
//----------------------------------------------------------------
// Overloaded () operator which only allows access to the array
//----------------------------------------------------------------
template <class datatype>
const datatype& Array1D<datatype>::operator()(int l) const
{
    return this->_array[l];
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
Array1D<datatype> Array1D<datatype>::operator + (datatype value)
{
    Array1D<datatype> output(this->_length);

    for (int i=0; i<this->_length; i++)
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
Array1D<datatype> Array1D<datatype>::operator + (const Array1D& array)
{
    Array1D<datatype> output(this->_length);

    for (int i=0; i<this->_length; i++)
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
Array1D<datatype> Array1D<datatype>::operator - (datatype value)
{
    Array1D<datatype> output(this->_length);

    for (int i=0; i<this->_length; i++)
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
Array1D<datatype> Array1D<datatype>::operator - (const Array1D& array)
{
    Array1D<datatype> output(this->_length);

    for (int i=0; i<this->_length; i++)
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
Array1D<datatype> Array1D<datatype>::operator * (datatype value)
{
    Array1D<datatype> output(this->_length);

    for (int i=0; i<this->_length; i++)
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
Array1D<datatype> Array1D<datatype>::operator * (const Array1D& array)
{
    Array1D<datatype> output(this->_length);

    for (int i=0; i<this->_length; i++)
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
Array1D<datatype> Array1D<datatype>::operator / (datatype value)
{
    Array1D<datatype> output(this->_length);

    for (int i=0; i<this->_length; i++)
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
Array1D<datatype> Array1D<datatype>::operator / (const Array1D& array)
{
    Array1D<datatype> output(this->_length);

    for (int i=0; i<this->_length; i++)
    {
        output._array[i] = this->_array[i] / array._array[i];
    }

    return output;
}


