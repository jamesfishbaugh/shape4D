//--------------------------------------------------------------------------
// Array3D:  Templated 3D array class stored as a 1D array (row-major)
//--------------------------------------------------------------------------

#ifndef ARRAY3D_H
#define ARRAY3D_H

#include "array2d.h"

template<class datatype>
class Array3D
{

    private:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------

        datatype* _array;						// The array which forms the basis of this class
        int _length;							// The length of the array
        int _width;								// The width of the array
        int _height;							// The height of the array
        int _lengthTimesWidth;					// Length times height
        int _widthTimesHeight;					// Width times height
        int _lengthTimesWidthTimesHeight;		// Total size of the array

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        Array3D();
        Array3D(int length, int width, int height);
        Array3D(const Array3D& array);
        ~Array3D();

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        int GetLength() const;
        int GetWidth() const;
        int GetHeight() const;
        datatype* GetArray();

        //--------------------------------------------------------------------------
        // Accessors
        //--------------------------------------------------------------------------

        datatype GetAt(int l, int w, int h) const;
        void SetAt(int l, int w, int h, datatype value);
        Array2D<datatype> Get2DSliceAtHeight(int h) const;
        void Set2DSliceAtHeight(const Array2D<datatype>& array, int h);

        //--------------------------------------------------------------------------
        // Fill methods
        //--------------------------------------------------------------------------

        void FillArray(datatype value);
        void FillArray(datatype* array);
        void FillArray(datatype*** array);

        //--------------------------------------------------------------------------
        // Array methods
        //--------------------------------------------------------------------------

        datatype Sum();

        //--------------------------------------------------------------------------
        // Overloaded operators
        //--------------------------------------------------------------------------

        Array3D<datatype>& operator = (const Array3D& array);
        bool operator == (const Array3D& array);
        bool operator != (const Array3D& array);
        datatype& operator()(int l, int w, int h);
        const datatype& operator()(int l, int w, int h) const;
        Array3D<datatype> operator + (datatype value);
        Array3D<datatype> operator + (const Array3D& array);
        Array3D<datatype> operator - (datatype value);
        Array3D<datatype> operator - (const Array3D& array) const;
        Array3D<datatype> operator * (datatype value) const;
        Array3D<datatype> operator * (const Array3D& array);
        Array3D<datatype> operator / (datatype value);
        Array3D<datatype> operator / (const Array3D& array);

};

#include "array3d.txx"

#endif // ARRAY3D_H


