//--------------------------------------------------------------------------
// MultiObjectComplex:  Collection of multiple ShapeObjects
//--------------------------------------------------------------------------

#ifndef MULTIOBJECTCOMPLEX_H
#define MULTIOBJECTCOMPLEX_H

#include "shapeobject.h"
#include <vector>

using namespace std;

class MultiObjectComplex
{

    private:

        vector<ShapeObject*> _shapes;

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructors
        //--------------------------------------------------------------------------

        MultiObjectComplex();
        MultiObjectComplex(const MultiObjectComplex& multiObject);
        ~MultiObjectComplex();

        //--------------------------------------------------------------------------
        // Accessors
        //--------------------------------------------------------------------------

        void AddShape(ShapeObject* shape);

        void SetShapeAt(int index, ShapeObject* shape);

        int GetNumberOfShapes();
        ShapeObject* GetShapeAt(int index);

};

#endif // MULTIOBJECTCOMPLEX_H

