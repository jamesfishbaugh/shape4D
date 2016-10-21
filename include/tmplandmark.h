//----------------------------------------------------------------
// tmpLandmarks: Shape data represented as landmarks
//----------------------------------------------------------------

#ifndef TMPLANDMARKS_H
#define TMPLANDMARKS_H

#include "shapeobject.h"

class tmpLandmarks : public ShapeObject
{

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        tmpLandmarks();
        tmpLandmarks(const tmpLandmarks& shape);
        ~tmpLandmarks();

        //--------------------------------------------------------------------------
        // Interface to compute data matching metric and gradient
        //--------------------------------------------------------------------------

        virtual double Matching(const ShapeObject* shape);
        virtual Array2D<double> GradMatching(const ShapeObject* shape);

        //--------------------------------------------------------------------------
        // Overloaded operators
        //--------------------------------------------------------------------------
        tmpLandmarks& operator = (const tmpLandmarks& shape);

};

#endif // TMPLANDMARKS_H
