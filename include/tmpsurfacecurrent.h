//----------------------------------------------------------------
// SurfaceCurrent: Shape data represented as a surface current
//----------------------------------------------------------------

#ifndef TMPSURFACECURRENT_H
#define TMPSURFACECURRENT_H

#include "shapeobject.h"

class tmpSurfaceCurrent : public ShapeObject
{

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        tmpSurfaceCurrent();
        tmpSurfaceCurrent(const tmpSurfaceCurrent& shape);
        ~tmpSurfaceCurrent();

        //--------------------------------------------------------------------------
        // Interface to compute data matching metric and gradient
        //--------------------------------------------------------------------------

        virtual double Matching(const ShapeObject* shape);
        virtual Array2D<double> GradMatching(const ShapeObject* shape);

        //--------------------------------------------------------------------------
        // Overloaded operators
        //--------------------------------------------------------------------------
        tmpSurfaceCurrent& operator = (const tmpSurfaceCurrent& shape);

};

#endif // TMPSURFACECURRENT_H
