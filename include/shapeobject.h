//--------------------------------------------------------------------------
// ShapeObject:  Virtual base class representing a shape object
//--------------------------------------------------------------------------

#ifndef SHAPEOBJECT_H
#define SHAPEOBJECT_H

#include "array1d.h"
#include "array2d.h"
#include "array3d.h"

class ShapeObject
{

    protected:

        //--------------------------------------------------------------------------
        // Protected member variables
        //--------------------------------------------------------------------------

        int _numPoints;					// Number of points
        Array2D<double> _points;        // Shape points
        int _numEdges;					// Number of triangles
        Array2D<int> _edges;			// Edges of the shape
        double _sigmaW;                 // Size of the kernel
        double _timept;                 // Time associated with the target
        int _timeIndex;                 // Time index associated with the target
        double _weight;                 // Relative importance of the target

    public:

        //--------------------------------------------------------------------------
        // Destructor
        //--------------------------------------------------------------------------
        ShapeObject();
        ShapeObject(const ShapeObject& shape);
        virtual ~ShapeObject();

        virtual ShapeObject* Copy() const = 0;

        //--------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------

        void SetPoints(const Array2D<double>& points);
        void SetEdges(const Array2D<int>& edges);
        void SetSigmaW(double sigmaW);
        void SetTimept(double timept);
        void SetTimeIndex(int timeIndex);
        void SetWeight(double weight);

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        int GetNumberOfPoints() const;
        const Array2D<double> GetPoints() const;
        int GetNumberOfEdges() const;
        const Array2D<int> GetEdges() const;
        double GetSigmaW() const;
        double GetTimept() const;
        double GetTimeIndex() const;
        double GetWeight() const;

        //--------------------------------------------------------------------------
        // Interface to compute data matching metric and gradient
        //--------------------------------------------------------------------------

        virtual double Matching(const Array2D<double>& shape) const = 0;
        virtual Array2D<double> GradMatching(const Array2D<double>& shape) const = 0 ;

};

#endif // SHAPEOBJECT_H
