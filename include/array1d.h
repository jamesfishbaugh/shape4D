//--------------------------------------------------------------------------
// Array1D:  Templated 1D array class
//--------------------------------------------------------------------------

#ifndef ARRAY1D_H
#define ARRAY1D_H

template<class datatype>
class Array1D
{

    private:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------

        datatype* _array;		// The array which forms the basis of this class
        int _length;			// The length of the array

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        Array1D();
        Array1D(int length);
        Array1D(const Array1D& array);
        ~Array1D();

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        int GetLength() const;
        datatype* GetArray();

        //--------------------------------------------------------------------------
        // Accessors
        //--------------------------------------------------------------------------

        datatype GetAt(int l);
        void SetAt(int l, datatype value);

        //--------------------------------------------------------------------------
        // Fill methods
        //--------------------------------------------------------------------------

        void FillArray(datatype value);
        void FillArray(datatype* array);

        //--------------------------------------------------------------------------
        // Array methods
        //--------------------------------------------------------------------------

        datatype Sum();

        //--------------------------------------------------------------------------
        // Overloaded operators
        //--------------------------------------------------------------------------

        Array1D<datatype>& operator = (const Array1D& array);
        bool operator == (const Array1D& array);
        bool operator != (const Array1D& array);
        datatype& operator()(int l);
        const datatype& operator()(int l) const;
        Array1D<datatype> operator + (datatype value);
        Array1D<datatype> operator + (const Array1D& array);
        Array1D<datatype> operator - (datatype value);
        Array1D<datatype> operator - (const Array1D& array);
        Array1D<datatype> operator * (datatype value);
        Array1D<datatype> operator * (const Array1D& array);
        Array1D<datatype> operator / (datatype value);
        Array1D<datatype> operator / (const Array1D& array);

};

#include "array1d.txx"

#endif // ARRAY1D_H


