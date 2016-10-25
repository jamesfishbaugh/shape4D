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

        virtual ShapeObject* Copy() const;

        //--------------------------------------------------------------------------
        // Interface to compute data matching metric and gradient
        //--------------------------------------------------------------------------

        virtual double Matching(const Array2D<double>& shape) const;
        virtual Array2D<double> GradMatching(const Array2D<double>& shape) const;

        //--------------------------------------------------------------------------
        // Overloaded operators
        //--------------------------------------------------------------------------
        tmpLandmarks& operator = (const tmpLandmarks& shape);

};

#endif // TMPLANDMARKS_H
