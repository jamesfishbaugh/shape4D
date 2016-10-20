//--------------------------------------------------------------------------
// Grid: A grid used for fast kernel computations
//--------------------------------------------------------------------------

#ifndef GRID_H
#define GRID_H

class Grid
{

    private:

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------

        double _ratio;			    // Accuracy of the grid
        double _pas;				// Accuracy of the grid
        double _origin[3];		    // Grid origin
        int _long[3];			    // Grid dimensions
        double* _fft3kD;			// FFT of the grid positions

        //--------------------------------------------------------------------------
        // Grid methods
        //--------------------------------------------------------------------------

        bool ChangeGrid(double mini[3],double maxi[3], double sigmaV);
        void SetGrid(double mini[3], double maxi[3], double sigmaV);

        //--------------------------------------------------------------------------
        // Helper functions
        //--------------------------------------------------------------------------

        void ComputeFFT3kD(double sigmaV);

    public:

        //--------------------------------------------------------------------------
        // Constructors/Destructor
        //--------------------------------------------------------------------------

        Grid();
        Grid(const Grid& grid);
        ~Grid();

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        double GetPas();
        double* GetOrigin();
        int* GetLong();
        double* GetFFT3KD();

        //--------------------------------------------------------------------------
        // Grid methods
        //--------------------------------------------------------------------------

        void UpdateGrid(double mini[3],double maxi[3], double sigmaV);

        //--------------------------------------------------------------------------
        // Overloaded operators
        //--------------------------------------------------------------------------

        Grid& operator = (const Grid& grid);

};

#endif // GRID_H
