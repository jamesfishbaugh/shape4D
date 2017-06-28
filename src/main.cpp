//--------------------------------------------------------------------------
// main.cpp:  Entry point for the application.
//--------------------------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>

#include "runexperiment.h"

#include "shapeobject.h"
#include "tmplandmark.h"
#include "tmpsurfacecurrent.h"
#include "multiobjectcomplex.h"

#include "shape4DCLP.h"

//----------------------------------------------------------------
// main
//----------------------------------------------------------------
// Inputs: 
//   argc - number of command line arguments
//	 argv - array of command line arguments   
//
// Outputs: 
//   return - value indicates success or failure
//----------------------------------------------------------------
// Entry point of application
//----------------------------------------------------------------
int main(int argc, char *argv[])
{	
    printf("\n");

    // Some testing
    ShapeObject* landmark1 = (ShapeObject*) new tmpLandmarks();
    Array2D<double> points1(3, 200);
    points1.FillArray(0.0);
    landmark1->SetPoints(points1);
    landmark1->SetSigmaW(20.0);
    landmark1->SetTimept(2.0);
    landmark1->SetTimeIndex(0);
    landmark1->SetWeight(1.0);

    ShapeObject* surface1 = (ShapeObject*) new tmpSurfaceCurrent();
    Array2D<double> points2(3, 360);
    points1.FillArray(0.0);
    Array2D<int> edges2(3,600);
    edges2.FillArray(1);
    surface1->SetPoints(points2);
    surface1->SetEdges(edges2);
    surface1->SetSigmaW(10.0);
    surface1->SetTimept(6.0);
    surface1->SetTimeIndex(30);
    surface1->SetWeight(1.0);

    MultiObjectComplex multiObject;
    multiObject.AddShape(landmark1);
    multiObject.AddShape(surface1);

    ShapeObject* newObject = multiObject.GetShapeAt(1);
    ShapeObject* anotherObject = multiObject.GetShapeAt(0);

    multiObject.SetShapeAt(0, newObject);
    multiObject.SetShapeAt(1, anotherObject);

    ShapeObject* firstObject = multiObject.GetShapeAt(0);
    ShapeObject* secondObject = multiObject.GetShapeAt(1);

    double matchValue = firstObject->Matching(secondObject->GetPoints());

    // Check command line arguments
    if (argc < 2)
    {
        printf("Usage:  shape4D driver_file.xml    OR\n        shape4D --continue progress_file.exo\n");
        exit(1);
    }

    RunExperiment experiment;
    if (argc == 2)
    {
        experiment.StartExperiment(argv[1]);
    }
    else
    {
        experiment.ContinueExperiment(argv[2]);
    }

    return 0;
}
