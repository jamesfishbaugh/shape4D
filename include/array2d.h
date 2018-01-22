//--------------------------------------------------------------------------
// Array2D:  Templated 2D array class stored as a 1D array (row-major)
//--------------------------------------------------------------------------

#ifndef ARRAY2D_H
#define ARRAY2D_H

template<class datatype>
class Array2D
{

    private:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------

        datatype* _array;			// The array which forms the basis of this class
        int _length;				// The length of the array
        int _width;					// The width of the array
        int _lengthTimesWidth;		// Total size of the array

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        Array2D();
        Array2D(int length, int width);
        Array2D(const Array2D& array);
        ~Array2D();

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        int GetLength() const;
        int GetWidth() const;
        datatype* GetArray();

        //--------------------------------------------------------------------------
        // Accessors
        //--------------------------------------------------------------------------

        datatype GetAt(int l, int w);
        void SetAt(int l, int w, datatype value);

        //--------------------------------------------------------------------------
        // Fill methods
        //--------------------------------------------------------------------------

        void FillArray(datatype value);
        void FillArray(datatype* array);
        void FillArray(datatype** array);

        //--------------------------------------------------------------------------
        // Array methods
        //--------------------------------------------------------------------------

        datatype Sum();

        //--------------------------------------------------------------------------
        // Overloaded operators
        //--------------------------------------------------------------------------

        Array2D<datatype>& operator = (const Array2D& array);
        bool operator == (const Array2D& array);
        bool operator != (const Array2D& array);
        datatype& operator()(int l, int w);
        const datatype& operator()(int l, int w) const;
        Array2D<datatype> operator + (datatype value);
        Array2D<datatype> operator + (const Array2D& array);
        Array2D<datatype> operator - (datatype value);
        Array2D<datatype> operator - (const Array2D& array) const;
        Array2D<datatype> operator * (datatype value) const;
        Array2D<datatype> operator * (const Array2D& array);
        Array2D<datatype> operator / (datatype value);
        Array2D<datatype> operator / (const Array2D& array);

};

#include "array2d.txx"

#endif // ARRAY2D_H


