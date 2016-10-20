//--------------------------------------------------------------------------
// SaveShapesAndVectors:  Save the current shapes and vectors
//--------------------------------------------------------------------------

#ifndef SAVESHAPESANDVECTORS_H
#define SAVESHAPESANDVECTORS_H

#include "regression.h"
#include "array3d.h"
#include <vector>

using namespace std;

class SaveShapesAndVectors
{

    private:

        static int _T;
        static int _numSourceShapes;
        static int* _numSourcePtsArray;
        static int* _numSourceTrisArray;
        static vector< Array2D<int> > _sourceTris;

    public:

        static void SetT(int T);
        static void SetNumSourceShapes(int numSourceShapes);
        static void SetNumSourcePtsArray(int* numSourcePtsArray);
        static void SetNumSourceTrisArray(int* numSourceTrisArray);
        static void SetSourceTris(vector< Array2D<int> > sourceTris);

        static bool SaveRegression(Regression* regression);

};

#endif // SAVESHAPESANDVECTORS_H
