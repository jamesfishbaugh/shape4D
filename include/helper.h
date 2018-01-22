//--------------------------------------------------------------------------
// Helper: Static class of helper functions for exoshape
//--------------------------------------------------------------------------

#ifndef HELPER_H
#define HELPER_H


class Helper
{

    public:

        static float** Allocate2DFloatArray(int rows, int cols);
        static float*** Allocate3DFloatArray(int rows, int cols, int time);

        static int** Allocate2DIntArray(int rows, int cols);
        static int*** Allocate3DIntArray(int rows, int cols, int time);

        static void Fill2DFloatArray(float value, int rows, int cols, float** &array);
        static void Fill3DFloatArray(float value, int rows, int cols, int time, float*** &array);

        static void Fill2DIntArray(int value, int rows, int cols, int** &array);
        static void Fill3DIntArray(int value, int rows, int cols, int time, int*** &array);

        static void Copy2DFloatArray(int rows, int cols, float** arrayFrom, float** &arrayTo);
        static void Copy3DFloatArray(int rows, int cols, int time, float*** arrayFrom, float*** &arrayTo);

        static void Copy2DIntArray(int rows, int cols, int** arrayFrom, int** &arrayTo);
        static void Copy3DIntArray(int rows, int cols, int time, int*** arrayFrom, int*** &arrayTo);

        static float Sum2DFloatArray(int rows, int cols, float** array);
        static float Sum3DFloatArray(int rows, int cols, int time, float*** array);

        static int Sum2DIntArray(int rows, int cols, int** array);
        static int Sum3DIntArray(int rows, int cols, int time, int*** array);

        static void Deallocate2DFloatArray(float** array);
        static void Deallocate3DFloatArray(float*** array);

        static void Deallocate2DIntArray(int** array);
        static void Deallocate3DIntArray(int*** array);
};

#endif // HELPER_H
