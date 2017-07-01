//--------------------------------------------------------------------------
// EXRunExperiment:  Driver for running shape regression experiments
//--------------------------------------------------------------------------

#include "runexperiment.h"

#include <time.h>			// For time
#include <vector>			// For vector
#include <stdlib.h>			// For file paths
#include <sys/stat.h>       // For checking if a directory exists
#include <string>
#include <fstream>

// Data structures
#include "array1d.h"
#include "array3d.h"

// File IO
#include "tinyxml.h"
#include "polydatareader.h"
#include "polydatawriter.h"
#include "shape4dstate.h"

// Data for regression
#include "array2d.h"
#include "regressionparams.h"
#include "targetdata.h"
#include "surfacecurrent.h"
#include "landmarks.h"

// Algorithms for regression
#include "regression.h"
#include "regressionvelocity.h"
#include "regressionacceleration.h"

#include "shapeobject.h"
#include "tmplandmark.h"
#include "tmpsurfacecurrent.h"
#include "multiobjectcomplex.h"

#ifdef _MSC_VER
#define realpath(N,R) _fullpath((R),(N),260)
#endif
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// EXRunExperiment
//----------------------------------------------------------------
// Inputs: 
//     
// Outputs: 
//----------------------------------------------------------------
// Default constructor 
//----------------------------------------------------------------
RunExperiment::RunExperiment()
{
    this->_continueExp = false;
    this->_skipRegressionVelocity = false;
}

//----------------------------------------------------------------
// ~EXRunExperiment
//----------------------------------------------------------------
// Inputs: 
//   
// Outputs: 
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
RunExperiment::~RunExperiment()
{
    //delete this->_initV0;
    //delete this->_impulse;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// EXPERIMENT DRIVERS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// RunExperiment
//----------------------------------------------------------------
// Inputs: 
//   pathToFile - path to XML driver file
//
// Outputs: 
//----------------------------------------------------------------
// Runs an experiment defined in an XML file
//----------------------------------------------------------------
void RunExperiment::StartExperiment(char* pathToFile)
{
    // Load the xml driver
    TiXmlDocument xmlDriver(pathToFile);
    bool loadOkay = xmlDriver.LoadFile();
    if (!loadOkay)
    {
        printf("***ERROR*** Could not load driver file: '%s'\n\n", xmlDriver.ErrorDesc());
        exit(1);
    }

    //--------------------------------------------------------------------------------
    // <experiment> </experiment>
    //--------------------------------------------------------------------------------

    // Grab the xml node named "experiment"
    TiXmlNode* experimentNode = xmlDriver.FirstChild("experiment");
    TiXmlNode* algNode = 0;

    // For using an initial velocity other than zero
    Array2D<double> initV0;
    // Check for warnings
    bool warnings = false;

    // Check to see if an "experiment" tag exists
    if (experimentNode == 0)
    {
        printf("***ERROR*** No <experiment> tag found\n\n");
        exit(1);
    }

    // Get the name of the experiment
    char* experimentName = (char*)experimentNode->ToElement()->Attribute("name");

    // If the experiment name exists
    if (experimentName != 0)
    {
        // Set the name of the .exo progress file
        Shape4DState::SetExperimentName(experimentName);
    }
    else
    {
        // Else give it a default name
        Shape4DState::SetExperimentName((char*)"progress");
        printf("***WARNING*** No name specified in <experiment name=\"\">.  Saving progress as \"progress.exo\".\n");
        warnings = true;
    }

    //--------------------------------------------------------------------------------
    // <algorithm> </algorithm>
    //--------------------------------------------------------------------------------

    int numAlgorithms = 0;
    for(algNode = experimentNode->FirstChild("algorithm"); algNode; algNode = algNode->NextSibling("algorithm"))
    {
        numAlgorithms++;
    }

    // Loop over the algorithms
    for(algNode = experimentNode->FirstChild("algorithm"); algNode; algNode = algNode->NextSibling("algorithm"))
    {

        //--------------------------------------------------------------------------------
        // <algorithm name="">
        //--------------------------------------------------------------------------------

        // Get the name which specifies which algorithm to run
        char* algorithmType = (char*)algNode->ToElement()->Attribute("name");
        // If none is specified we will run an acceleration controled regression
        if (algorithmType == 0)
        {
            printf("***WARNING*** No name specified in <algorithm name=\"\">.  Using \"RegressionAccel\".\n");
            algorithmType = (char*) "RegressionAccel";
            warnings = true;
        }
        if ((strcmp(algorithmType,"RegressionVelocity") != 0) && (strcmp(algorithmType,"RegressionAccel") != 0))
        {
            printf("***WARNING*** No algorithm named \"%s\".  Using \"RegressionAccel\".\n", algorithmType);
            algorithmType = (char*) "RegressionAccel";
            warnings = true;
        }

        // If we are continuing and the RegressionVelocity portion has already completed, we skip it
        if ((strcmp(algorithmType,"RegressionVelocity") == 0))
        {
            if (this->_skipRegressionVelocity)
            {
                continue;
            }
        }

        //--------------------------------------------------------------------------------
        // <source> </source>
        //--------------------------------------------------------------------------------

        // Read the "source" element
        TiXmlNode* sourceNode = algNode->FirstChild("source");
        if (sourceNode == 0)
        {
            printf("***ERROR*** No <source> tag found.  A source is required.\n\n");
            exit(1);
        }

        //--------------------------------------------------------------------------------
        // <input> </input>
        //--------------------------------------------------------------------------------

        // Read the "input" element
        TiXmlNode* inputNode = sourceNode->FirstChild("input");
        if (inputNode == 0)
        {
            printf("***ERROR*** No <input> tag found.  A source requires input data.\n\n");
            exit(1);
        }

        //--------------------------------------------------------------------------------
        // <shape> </shape>
        //--------------------------------------------------------------------------------

        // Figure out how many shapes make up the source
        TiXmlNode* node = 0;
        int numSourceShapes = 0;
        for(node = inputNode->FirstChild("shape"); node; node = node->NextSibling("shape"))
        {
            numSourceShapes++;
        }
        if (numSourceShapes == 0)
        {
            printf("***ERROR*** No <shape> tag found.  A source must include at least one path to data.\n\n");
            exit(1);
        }

        // Read the source data paths and store them in an array
        char** sourcePathArray = new char*[numSourceShapes];
        node = 0;
        int count = 0;
        for(node = inputNode->FirstChild("shape"); node; node = node->NextSibling("shape"))
        {
            sourcePathArray[count] = (char*)node->ToElement()->GetText();
            count++;
        }

        //--------------------------------------------------------------------------------
        // <output> </output>
        //--------------------------------------------------------------------------------

        // Read the "output" element
        TiXmlNode* outputNode = sourceNode->FirstChild("output");
        if (outputNode == 0)
        {
            printf("***WARNING*** No <output> tag found.  Results will not be saved.\n");
            warnings = true;
            outputNode = (TiXmlNode*) new TiXmlElement("NO_OUTPUT");
            Shape4DState::SetShouldSave(false);
        }

        //--------------------------------------------------------------------------------
        // <saveProgress> </saveProgress>
        //--------------------------------------------------------------------------------

        int saveProgressIters = 0;
        char* saveProgress = (char*)"";
        if (strcmp(outputNode->Value(), "NO_OUTPUT") != 0)
        {
            if (outputNode->FirstChild("saveProgress") == 0)
            {
                printf("***WARNING*** No <saveProgress> tag found.  Progress will not be saved.\n");
                warnings = true;
                Shape4DState::SetShouldSave(false);
            }
            else
            {
                saveProgress = (char*)outputNode->FirstChild("saveProgress")->ToElement()->GetText();
                sscanf(saveProgress, "%d", &saveProgressIters);
                Shape4DState::SetShouldSave(true);
                Shape4DState::SetSaveProgressEveryNIterations(saveProgressIters);
            }
        }


        //--------------------------------------------------------------------------------
        // <dir> </dir>
        //--------------------------------------------------------------------------------

        char* outputDir = (char*)"";
        if (strcmp(outputNode->Value(), "NO_OUTPUT") != 0)
        {
            if (outputNode->FirstChild("dir") == 0)
            {
                printf("***ERROR*** No <dir> tag found.  Please include a directory for output.\n\n");
                exit(1);
            }
            else
            {
                outputDir = (char*)outputNode->FirstChild("dir")->ToElement()->GetText();

                // Make sure output directory exists
                struct stat st;
                if(stat(outputDir,&st) != 0)
                {
                    printf("***ERROR*** Output directory \"%s\" does not exist, please specifiy an existing directory for output.\n\n", outputDir);
                    exit(1);
                }
            }
        }

        //--------------------------------------------------------------------------------
        // <prefix> </prefix>
        //--------------------------------------------------------------------------------

        // Read the output paths and store them in an array
        char** outputPathArray = new char*[numSourceShapes];
        // First initialize the values to "" in case outputs aren't included
        // Later we check if the outputPath is "" to decide if we need to write results
        for (int s=0; s<numSourceShapes; s++)
        {
            outputPathArray[s] = (char*)"";
        }
        node = 0;
        count = 0;
        for(node = outputNode->FirstChild("prefix"); node; node = node->NextSibling("prefix"))
        {
            outputPathArray[count] = (char*)node->ToElement()->GetText();
            count++;
        }
        // Let the user know some or all of the results will not be saved
        if (count < numSourceShapes)
        {
            if (count == 0)
            {
                printf("***WARNING*** No <prefix> tags found.  Shape results will not be saved.\n");
            }
            else
            {
                printf("***WARNING*** Less output <prefix> tags than input <shape> tags.  Some shape results will not be saved.\n");
            }
            warnings = true;
        }

        //--------------------------------------------------------------------------------
        // <sigmaV> </sigmaV>
        //--------------------------------------------------------------------------------

        double sigmaV;
        if (sourceNode->FirstChild("sigmaV") == 0)
        {
            printf("***ERROR*** No <sigmaV> tag found.  A source must include a value for sigmaV.\n\n");
            exit(1);
        }
        char* sigmaVString = (char*)sourceNode->FirstChild("sigmaV")->ToElement()->GetText();
        sscanf(sigmaVString, "%lf", &sigmaV);

        //--------------------------------------------------------------------------------
        // <gammaR> </gammaR>
        //--------------------------------------------------------------------------------

        double gammaR;
        if (sourceNode->FirstChild("gammaR") == 0)
        {
            printf("***ERROR*** No <gammaR> tag found.  A source must include a value for gammaR.\n\n");
            exit(1);
        }
        char* gammaRString = (char*)sourceNode->FirstChild("gammaR")->ToElement()->GetText();
        sscanf(gammaRString, "%lf", &gammaR);

        //--------------------------------------------------------------------------------
        // <t0> </t0>
        //--------------------------------------------------------------------------------

        double t0;
        if (sourceNode->FirstChild("t0") == 0)
        {
            printf("***ERROR*** No <t0> tag found.  A source must include an initial time point.\n\n");
            exit(1);
        }
        char* t0String = (char*)sourceNode->FirstChild("t0")->ToElement()->GetText();
        sscanf(t0String, "%lf", &t0);

        //--------------------------------------------------------------------------------
        // <tn> </tn>
        //--------------------------------------------------------------------------------

        double tn;
        if (sourceNode->FirstChild("tn") == 0)
        {
            printf("***ERROR*** No <tn> tag found.  A source must include a final time point.\n\n");
            exit(1);
        }
        char* tnString = (char*)sourceNode->FirstChild("tn")->ToElement()->GetText();
        sscanf(tnString, "%lf", &tn);

        //--------------------------------------------------------------------------------
        // <T> </T>
        //--------------------------------------------------------------------------------

        int T;
        if (sourceNode->FirstChild("T") == 0)
        {
            printf("***ERROR*** No <T> tag found.  A source must include the number of points for the time discretization.\n\n");
            exit(1);
        }
        char* TString = (char*)sourceNode->FirstChild("T")->ToElement()->GetText();
        sscanf(TString, "%d", &T);

        //--------------------------------------------------------------------------------
        // <useInitV0> </useInitV0>
        //--------------------------------------------------------------------------------

        // Attempt to read the "useInitV0" option
        int usePGV0 = 0;
        if (sourceNode->FirstChild("useInitV0") != 0)
        {
            char* usePGV0String = (char*)sourceNode->FirstChild("useInitV0")->ToElement()->GetText();
            sscanf(usePGV0String, "%d", &usePGV0);

            if ((numAlgorithms == 1) && (usePGV0 == 1))
            {
                printf("***WARNING*** <useInitV0> 1 </useInitV0> specified but no <algorithm name=\"RegressionVelocity\"> found.  Using 0 for initial velocity instead.\n");
                usePGV0 = 0;
                warnings = true;
            }
        }

        //--------------------------------------------------------------------------------
        // <v0weight> </v0weight>
        //--------------------------------------------------------------------------------
        double v0weight = 1.0;
        if (sourceNode->FirstChild("v0weight") != 0)
        {
            char* v0weightString = (char*)sourceNode->FirstChild("v0weight")->ToElement()->GetText();
            sscanf(v0weightString, "%lf", &v0weight);

            //if ((v0weight != 0) && (usePGV0 == 1))
            //{
            //	printf("***WARNING*** Conflict between <useInitV0> 1 </useInitV0> and <v0weight> >0 </v0weight>.  Using 0 for v0weight.\n");
            //	usePGV0 = 0;
            //	warnings = true;
            //}
        }

        //--------------------------------------------------------------------------------
        // <estimateBaseline> </estimateBaseline>
        //--------------------------------------------------------------------------------

        // Attempt to read the "estimateBaseline" option
        int estimateBaseline = 0;
        if (sourceNode->FirstChild("estimateBaseline") != 0)
        {
            char* estimateBaselineString = (char*)sourceNode->FirstChild("estimateBaseline")->ToElement()->GetText();
            sscanf(estimateBaselineString, "%d", &estimateBaseline);
        }

        //--------------------------------------------------------------------------------
        // <kernelType> </kernelType>
        //--------------------------------------------------------------------------------

        // Attempt to read the "kernelType" option
        char* kernelType;
        if (sourceNode->FirstChild("kernelType") != 0)
        {
            kernelType = (char*)sourceNode->FirstChild("kernelType")->ToElement()->GetText();

            if ((strcmp(kernelType,"exact") != 0) && (strcmp(kernelType,"p3m") != 0))
            {
                printf("***WARNING*** Available kernels types are \"exact\" and \"p3m\". Defaulting to \"exact\" kernel\n");
                kernelType = (char*) "exact";
                warnings = true;
            }
        }
        else
        {
            kernelType = (char*) "exact";
        }

        //--------------------------------------------------------------------------------
        // <maxIters> </maxIters>
        //--------------------------------------------------------------------------------
        int maxIters = 500;
        if (sourceNode->FirstChild("maxIters") != 0)
        {
            char* maxItersString = (char*)sourceNode->FirstChild("maxIters")->ToElement()->GetText();
            sscanf(maxItersString, "%d", &maxIters);
        }

        //--------------------------------------------------------------------------------
        // <breakRatio> </breakRatio>
        //--------------------------------------------------------------------------------
        double breakRatio = 1e-6;
        if (sourceNode->FirstChild("breakRatio") != 0)
        {
            char* breakRatioString = (char*)sourceNode->FirstChild("breakRatio")->ToElement()->GetText();
            sscanf(breakRatioString, "%lf", &breakRatio);
        }

        //--------------------------------------------------------------------------------
        // <useFista> </useFista>
        //--------------------------------------------------------------------------------

        // Attempt to read the "useFista" option
        int useFista = 1;
        if (sourceNode->FirstChild("useFista") != 0)
        {
            char* useFistaString = (char*)sourceNode->FirstChild("useFista")->ToElement()->GetText();
            sscanf(useFistaString, "%d", &useFista);
        }

        //--------------------------------------------------------------------------------
        // <baselineSmoothing> </baselineSmoothing>
        //--------------------------------------------------------------------------------

        // Attempt to read the "baselineSmoothing" option
        double baselineSmoothing = 1.0;
        if (sourceNode->FirstChild("baselineSmoothing") != 0)
        {
            char* baselineSmoothingString = (char*)sourceNode->FirstChild("baselineSmoothing")->ToElement()->GetText();
            sscanf(baselineSmoothingString, "%lf", &baselineSmoothing);
        }

        //--------------------------------------------------------------------------------
        // We are done reading source information and can now read source data and
        // build the source object
        //--------------------------------------------------------------------------------

        // Vectors to store all arrays of source pts and tris
        vector< Array2D<double> > allSourcePts(numSourceShapes);
        vector< Array2D<int> > allSourceTris(numSourceShapes);

        // Keep track of the sums and individual numbers of pts and tris
        int totalSourcePts = 0;
        int* numSourcePtsArray = new int[numSourceShapes];
        int* numSourceTrisArray = new int[numSourceShapes];

        MultiObjectComplex theSource;

        // Read all the source data files
        for (int i=0; i<numSourceShapes; i++)
        {
            Array2D<double> sourcePts;
            Array2D<int> sourceTris;
            VTKPolyDataReader sourceReader(sourcePathArray[i]);
            bool didReadSource = sourceReader.ReadPointsAndTris(sourcePts, sourceTris);
            // Did we fail to read the data file
            if (!didReadSource)
            {
                printf("***ERROR*** Could not read source VTK file: \"%s\"\n\n", sourcePathArray[i]);
                exit(1);
            }

            // Keep track of the number of pts and tris
            totalSourcePts += sourcePts.GetWidth();
            numSourcePtsArray[i] = sourcePts.GetWidth();
            numSourceTrisArray[i] = sourceTris.GetWidth();

            // Add current arrays to the vectors
            allSourcePts[i] = sourcePts;
            allSourceTris[i] = sourceTris;

            ShapeObject* curSourceShape = (ShapeObject*) new tmpSurfaceCurrent();
            curSourceShape->SetPoints(sourcePts);
            curSourceShape->SetEdges(sourceTris);
            curSourceShape->SetTimeIndex(0);
            curSourceShape->SetTimept(t0);

            theSource.AddShape(curSourceShape);
        }

        // Now lets put all the source points into one array
        Array2D<double> sourcePts(3, totalSourcePts);
        count = 0;
        for (int i=0; i<numSourceShapes; i++)
        {
            Array2D<double> curPts = allSourcePts[i];

            for (int j=0; j<curPts.GetWidth(); j++)
            {
                sourcePts(0,count) = curPts(0,j);
                sourcePts(1,count) = curPts(1,j);
                sourcePts(2,count) = curPts(2,j);
                count++;
            }
        }

        // Build the time discretization array
        double inc = (tn-t0)/(T-1);
        Array1D<double> xtime(T);
        double curTime = t0;
        for (int i=0; i<T; i++)
        {
            xtime(i) = curTime;
            curTime += inc;
        }

        // Create the source object
        RegressionParams source(sourcePts, sigmaV, gammaR, xtime);
        source.SetKernelType(kernelType);
        source.SetMaxIters(maxIters);
        source.SetBreakRatio(breakRatio);
        source.SetBaselineSmoothing(baselineSmoothing);

        //--------------------------------------------------------------------------------
        // <targets> </targets>
        //--------------------------------------------------------------------------------

        TiXmlNode* targetsNode = algNode->FirstChild("targets");
        if (targetsNode == 0)
        {
            printf("***ERROR*** No <targets> tag found.  At least one target is required.\n\n");
            exit(1);
        }

        //--------------------------------------------------------------------------------
        // <target> </target>
        //--------------------------------------------------------------------------------

        // First count the number of targets
        node = 0;
        int numTargets = 0;
        for(node = targetsNode->FirstChild("target"); node; node = node->NextSibling())
        {
            numTargets++;
        }
        if (numTargets == 0)
        {
            printf("***ERROR*** No <target> tag found.  At least one target is required.\n\n");
            exit(1);
        }

        // Create the target array and read target data
        count = 0;
        TargetData** targetArray = new TargetData*[numTargets];
        for(node = targetsNode->FirstChild("target"); node; node = node->NextSibling())
        {
            //--------------------------------------------------------------------------------
            // <shape> </shape>
            //--------------------------------------------------------------------------------

            if (node->FirstChild("shape") == 0)
            {
                printf("***ERROR*** No <shape> tag found.  A target must include a path to data.\n\n");
                exit(1);
            }
            char* targetDataPath = (char*)node->FirstChild("shape")->ToElement()->GetText();

            //--------------------------------------------------------------------------------
            // <type> </type>
            //--------------------------------------------------------------------------------

            if (node->FirstChild("type") == 0)
            {
                printf("***ERROR*** No <type> tag found.  A target must include a type.  Currently supported types are \"SURFACE\" and \"LANDMARKS\".\n\n");
                exit(1);
            }
            // Read target type which specifies the data representation and metric for the target
            char* targetTypeString = (char*)node->FirstChild("type")->ToElement()->GetText();
            // Convert the targetTypeString to its #DEFINE value
            int targetType;
            if (strcmp(targetTypeString, "SURFACE") == 0)
            {
                targetType = SURFACE;
            }
            else if (strcmp(targetTypeString, "LANDMARKS") == 0)
            {
                targetType = LANDMARKS;
            }
            else
            {
                printf("***ERROR*** Unknown target type: \"%s\".  Currently supported types are \"SURFACE\" and \"LANDMARKS\".\n\n", targetTypeString);
                exit(1);
            }

            //--------------------------------------------------------------------------------
            // <tris> </tris>
            //--------------------------------------------------------------------------------

            int tris = 0;
            // No source connectivity information provided
            if (node->FirstChild("tris") == 0)
            {
                printf("***WARNING*** No <tris> tag found.  Assuming connectivity with first source data.  This could lead to strange behavior.\n");
                warnings = true;
            }
            else
            {
                char* trisString = (char*)node->FirstChild("tris")->ToElement()->GetText();
                sscanf(trisString, "%d", &tris);
            }

            //--------------------------------------------------------------------------------
            // <sigmaW> </sigmaW>
            //--------------------------------------------------------------------------------

            double sigmaW;
            if (node->FirstChild("sigmaW") == 0)
            {
                if (targetType != LANDMARKS)
                {
                    printf("***ERROR*** No <sigmaW> tag found.  A target must include a value for sigmaW.\n\n");
                    exit(1);
                }
            }
            else
            {
                char* sigmaWString = (char*)node->FirstChild("sigmaW")->ToElement()->GetText();
                sscanf(sigmaWString, "%lf", &sigmaW);
            }

            //--------------------------------------------------------------------------------
            // <timept> </timept>
            //--------------------------------------------------------------------------------

            double timept;
            if (node->FirstChild("timept") == 0)
            {
                printf("***ERROR*** No <timept> tag found.  A target must include a time point.\n\n");
                exit(1);
            }
            char* timeptString = (char*)node->FirstChild("timept")->ToElement()->GetText();
            sscanf(timeptString, "%lf", &timept);

            if ((timept > tn) || (timept < t0))
            {
                printf("***ERROR*** <timept> %0.2f </timept> is out of the time range [%0.2f, %0.2f].\n\n", timept, t0, tn);
                exit(1);
            }

            //--------------------------------------------------------------------------------
            // <weight> </weight>
            //--------------------------------------------------------------------------------

            // Check to see if the user included a weight for this target
            double weight = 1.0;
            if (node->FirstChild("weight") != 0)
            {
                char* weightString = (char*)node->FirstChild("weight")->ToElement()->GetText();
                sscanf(weightString, "%lf", &weight);
            }

            //--------------------------------------------------------------------------------
            // We are done reading target information and can now read target data and
            // build the target object
            //--------------------------------------------------------------------------------

            // Read the target data
            Array2D<double> targetPts;
            Array2D<int> targetTris;
            VTKPolyDataReader targetReader(targetDataPath);
            bool didReadTarget = targetReader.ReadPointsAndTris(targetPts, targetTris);
            //printf("Num points = %d\n", targetTris.GetWidth());

            // Did we fail to read the data file
            if (!didReadTarget)
            {
                printf("***ERROR*** Could not read target VTK file: \"%s\"\n\n", targetDataPath);
                exit(1);
            }

            // Build the target object and add it to our array
            Array2D<int> sourceConnectivity = allSourceTris[tris];
            // Add the offset necessary since the source array can be a concatination of multiple shape data
            if (tris != 0)
            {
                int offset = 0;
                for (int i=(tris-1); i>=0; i--)
                {
                    offset += numSourcePtsArray[i];
                }

                sourceConnectivity = sourceConnectivity + offset;
            }

            // Compute the time index associated with this target
            int timeIndex = int((timept-t0)*((T-1)/(tn-t0))+0.5f);

            if (targetType == SURFACE)
            {
                targetArray[count] = (TargetData*) new SurfaceCurrent(targetPts, targetTris, sourceConnectivity, sigmaW, timept, timeIndex, weight);
            }
            else if (targetType == LANDMARKS)
            {
                targetArray[count] = (TargetData*) new Landmarks(targetPts, targetTris, sourceConnectivity, sigmaW, timept, timeIndex, weight);
            }
            else
            {
                printf("***ERROR*** Unknown target type: \"%s\".  Currently supported types are \"SURFACE\" and \"LANDMARKS\".\n", targetTypeString);
                exit(1);
            }

            targetArray[count]->SetKernelType(kernelType);

            count++;
        }

//        // Create the target array and read target data
//        int curShapeTriIndex = 0;
//        count = 0;
//        TargetData** targetArray = new TargetData*[numTargets];
//        vector<MultiObjectComplex> theTarget;
//        for(node = targetsNode->FirstChild("target"); node; node = node->NextSibling())
//        {
//            // No shape tag found, exit with an error
//            if (node->FirstChild("shape") == 0)
//            {
//                printf("***ERROR*** No <shape> tag found. A target must include shape information.\n\n");
//                exit(1);
//            }

//            MultiObjectComplex curMultiObject;
//            curShapeTriIndex = 0;

//            //--------------------------------------------------------------------------------
//            // <timept> </timept>
//            //--------------------------------------------------------------------------------

//            double timept;
//            if (node->FirstChild("timept") == 0)
//            {
//                printf("***ERROR*** No <timept> tag found.  A target must include a time point.\n\n");
//                exit(1);
//            }
//            char* timeptString = (char*)node->FirstChild("timept")->ToElement()->GetText();
//            sscanf(timeptString, "%lf", &timept);

//            // Compute the time index associated with this target
//            int timeIndex = int((timept-t0)*((T-1)/(tn-t0))+0.5f);

//            TiXmlNode* shapeNode = 0;
//            int numShapesHere = 0;
//            for (shapeNode = node->FirstChild("shape"); shapeNode; shapeNode = node->NextSibling("shape"))
//            {
//                numShapesHere++;
//            }

//            //--------------------------------------------------------------------------------
//            // <shape> </shape>
//            //--------------------------------------------------------------------------------
//            shapeNode = 0;
//            for (shapeNode = node->FirstChild("shape"); shapeNode; shapeNode = shapeNode->NextSibling("shape"))
//            {

//                //--------------------------------------------------------------------------------
//                // <path> </path>
//                //--------------------------------------------------------------------------------
//                if (shapeNode->FirstChild("path") == 0)
//                {
//                    printf("***ERROR*** No <path> tag found.  A target shape must include a path to data.\n\n");
//                    exit(1);
//                }
//                char* targetDataPath = (char*)shapeNode->FirstChild("path")->ToElement()->GetText();

//                //--------------------------------------------------------------------------------
//                // <type> </type>
//                //--------------------------------------------------------------------------------

//                if (shapeNode->FirstChild("type") == 0)
//                {
//                    printf("***ERROR*** No <type> tag found.  A target shape must include a type.  Currently supported types are \"SURFACE\" and \"LANDMARKS\".\n\n");
//                    exit(1);
//                }

//                // Read target type which specifies the data representation and metric for the target
//                char* targetTypeString = (char*)shapeNode->FirstChild("type")->ToElement()->GetText();
//                // Convert the targetTypeString to its #DEFINE value
//                int targetType;
//                if (strcmp(targetTypeString, "SURFACE") == 0)
//                {
//                    targetType = SURFACE;
//                }
//                else if (strcmp(targetTypeString, "LANDMARKS") == 0)
//                {
//                    targetType = LANDMARKS;
//                }
//                else
//                {
//                    printf("***ERROR*** Unknown target type: \"%s\".  Currently supported types are \"SURFACE\" and \"LANDMARKS\".\n\n", targetTypeString);
//                    exit(1);
//                }

//                //--------------------------------------------------------------------------------
//                // <sigmaW> </sigmaW>
//                //--------------------------------------------------------------------------------

//                double sigmaW;
//                if (shapeNode->FirstChild("sigmaW") == 0)
//                {
//                    if (targetType != LANDMARKS)
//                    {
//                        printf("***ERROR*** No <sigmaW> tag found.  A target shape must include a value for sigmaW.\n\n");
//                        exit(1);
//                    }
//                }
//                else
//                {
//                    char* sigmaWString = (char*)shapeNode->FirstChild("sigmaW")->ToElement()->GetText();
//                    sscanf(sigmaWString, "%lf", &sigmaW);
//                }

//                //--------------------------------------------------------------------------------
//                // <weight> </weight>
//                //--------------------------------------------------------------------------------

//                // Check to see if the user included a weight for this target
//                double weight = 1.0;
//                if (shapeNode->FirstChild("weight") != 0)
//                {
//                    char* weightString = (char*)shapeNode->FirstChild("weight")->ToElement()->GetText();
//                    sscanf(weightString, "%lf", &weight);
//                }

//                //--------------------------------------------------------------------------------
//                // We are done reading target information and can now read target data and
//                // build the target object
//                //--------------------------------------------------------------------------------

//                // Read the target data
//                Array2D<double> targetPts;
//                Array2D<int> targetTris;

//                VTKPolyDataReader targetReader(targetDataPath);
//                bool didReadTarget = targetReader.ReadPointsAndTris(targetPts, targetTris);

//                //printf("Num points = %d\n", targetTris.GetWidth());

//                // Did we fail to read the data file
//                if (!didReadTarget)
//                {
//                    printf("***ERROR*** Could not read target VTK file: \"%s\"\n\n", targetDataPath);
//                    exit(1);
//                }

//                // Build the target object and add it to our array
//                Array2D<int> sourceConnectivity = allSourceTris[curShapeTriIndex];
//                // Add the offset necessary since the source array can be a concatination of multiple shape data
//                if (curShapeTriIndex != 0)
//                {
//                    int offset = 0;
//                    for (int i=(curShapeTriIndex-1); i>=0; i--)
//                    {
//                        offset += numSourcePtsArray[i];
//                    }

//                    sourceConnectivity = sourceConnectivity + offset;
//                }

//                if (targetType == SURFACE)
//                {
//                    targetArray[count] = (TargetData*) new SurfaceCurrent(targetPts, targetTris, sourceConnectivity, sigmaW, timept, timeIndex, weight);
//                }
//                else if (targetType == LANDMARKS)
//                {
//                    targetArray[count] = (TargetData*) new Landmarks(targetPts, targetTris, sourceConnectivity, sigmaW, timept, timeIndex, weight);
//                }
//                else
//                {
//                    printf("***ERROR*** Unknown target type: \"%s\".  Currently supported types are \"SURFACE\" and \"LANDMARKS\".\n", targetTypeString);
//                    exit(1);
//                }

//                // Build the new target structures
//                ShapeObject* curShapeObject;

//                if (targetType == SURFACE)
//                {
//                    curShapeObject = (ShapeObject*) new tmpSurfaceCurrent();
//                }
//                else if (targetType == LANDMARKS)
//                {
//                    curShapeObject = (ShapeObject*) new tmpLandmarks();

//                }

//                curShapeObject->SetPoints(targetPts);
//                curShapeObject->SetEdges(targetTris);
//                curShapeObject->SetSigmaW(sigmaW);
//                curShapeObject->SetTimept(timept);
//                curShapeObject->SetTimeIndex(timeIndex);
//                curShapeObject->SetWeight(weight);

//                curMultiObject.AddShape(curShapeObject);

//                count++;
//                curShapeTriIndex++;

//            }

//            theTarget.push_back(curMultiObject);
//        }

        //--------------------------------------------------------------------------------
        // We have built the source and array of target objects and are ready to
        // build and run the specified algorithm
        //--------------------------------------------------------------------------------

        if (warnings)
        {
            printf("\n");
        }

        warnings = false;

        Regression* regression;
        Array3D<double> X;
        Array3D<double> dX;
        Array3D<double> momenta;
        Array3D<double> accel;
        Array3D<double> impulse;

        // Run piecewise geodesic regression parameterized by velocity
        if (strcmp(algorithmType, "RegressionVelocity") == 0)
        {
            // Set up the static object that saves the state and intermediate shapes and vectors during the optimization
            char *fullXMLPath = realpath(pathToFile, NULL);
            // If everything went well getting the full path
            if (fullXMLPath)
            {
                Shape4DState::SetXMLDriver(fullXMLPath);
            }
            else
            {
                // Just use the relative path
                Shape4DState::SetXMLDriver(pathToFile);
            }

            free(fullXMLPath);

            Shape4DState::SetT(T);
            Shape4DState::SetNumSourceShapes(numSourceShapes);
            Shape4DState::SetNumSourcePtsArray(numSourcePtsArray);
            Shape4DState::SetNumSourceTrisArray(numSourceTrisArray);
            Shape4DState::SetSourceTris(allSourceTris);

            Shape4DState::SetOutputDir(outputDir);
            vector<char*> outputNames;
            for (int s=0; s<numSourceShapes; s++)
            {
                outputNames.push_back(outputPathArray[s]);
            }
            Shape4DState::SetOutputNames(outputNames);

            // If we are continuing the experiment instead of starting from the beginning
            if (this->_continueExp)
            {

                printf("Continuing piecewise geodesic shape regression...\n");
                printf("SigmaV = %0.4f   GammaR = %0.4f   SigmaW = %0.4f\n", sigmaV, gammaR, ((TargetData*)targetArray[0])->GetSigmaW());

                vector<double> dataMatching = Shape4DState::GetDataMatchingValues();
                vector<double> regularity = Shape4DState::GetRegularityValues();
                vector<double> stepsizes = Shape4DState::GetStepsizeValues();

                for (unsigned int iters=0; iters<dataMatching.size(); iters++)
                {
                    printf("Iteration %3d   funct = %0.4f   data = %0.4f   reg = %0.4f   step = %0.10lf\n",
                           iters+1, dataMatching[iters] + source.GetGamma()*regularity[iters], dataMatching[iters], regularity[iters], stepsizes[iters]);
                }

            }
            else
            {
                Shape4DState::ResetOptimizationState();
                printf("Starting piecewise geodesic shape regression...\n");
                printf("SigmaV = %0.4f   GammaR = %0.4f   SigmaW = %0.4f\n", sigmaV, gammaR, ((TargetData*)targetArray[0])->GetSigmaW());
            }

            regression = (Regression*) new RegressionVelocity(source, numTargets, targetArray);

            // We set the momenta to continue where we left off
            if (_continueExp)
            {
                ((RegressionVelocity*)regression)->SetMomenta(*(_momenta));
            }

            time_t startTime = time(NULL);
            // Start the PG regression
            X = regression->Run();
            time_t endTime = time(NULL);
            double regressionTime = difftime(endTime, startTime);

            initV0 = ((RegressionVelocity*)regression)->GetV0();
            dX = ((RegressionVelocity*)regression)->GetVelocity();
            momenta = ((RegressionVelocity*)regression)->GetMomenta();

            printf("\nElapsed time is %0.2f seconds\n\n", regressionTime);
        }
        // Run regression parameterized by acceleration
        else if (strcmp(algorithmType, "RegressionAccel") == 0)
        {
            // Set up the static object that saves the state and intermediate shapes and vectors during the optimization
            char *fullXMLPath = realpath(pathToFile, NULL);
            // If everything went well getting the full path
            if (fullXMLPath)
            {
                Shape4DState::SetXMLDriver(fullXMLPath);
            }
            else
            {
                // Just use the relative path
                Shape4DState::SetXMLDriver(pathToFile);
            }

            //free(fullXMLPath);

            Shape4DState::SetT(T);
            Shape4DState::SetNumSourceShapes(numSourceShapes);
            Shape4DState::SetNumSourcePtsArray(numSourcePtsArray);
            Shape4DState::SetNumSourceTrisArray(numSourceTrisArray);
            Shape4DState::SetSourceTris(allSourceTris);

            Shape4DState::SetOutputDir(outputDir);
            vector<char*> outputNames;
            for (int s=0; s<numSourceShapes; s++)
            {
                outputNames.push_back(outputPathArray[s]);
            }
            Shape4DState::SetOutputNames(outputNames);

            if (usePGV0)
            {
                source.SetInitV0(initV0);
                source.SetShouldUseInitV0(true);
                source.SetV0Weight(0.0);
            }
            else
            {
                source.SetV0Weight(v0weight);
            }

            if (estimateBaseline)
            {
                source.SetShouldEstimateBaseline(true);
            }
            else
            {
                source.SetShouldEstimateBaseline(false);
            }

            if (useFista)
            {
                source.SetShouldUseFista(true);
            }
            else
            {
                source.SetShouldUseFista(false);
            }

            // If we are continuing the experiment instead of starting from the beginning
            if (this->_continueExp)
            {
                // We will use only the initial velocity
                if (usePGV0)
                {
                    source.SetInitV0(*(_initV0));
                    source.SetShouldUseInitV0(true);
                    source.SetV0Weight(0.0);
                }
                else
                {
                    source.SetInitV0(*(_initV0));
                    source.SetShouldUseInitV0(false);
                    source.SetV0Weight(v0weight);
                }

                printf("Continuing acceleration controlled shape regression...\n");
                printf("SigmaV = %0.4f   GammaR = %0.4f   SigmaW = %0.4f\n\n", sigmaV, gammaR, ((TargetData*)targetArray[0])->GetSigmaW());

                vector<double> dataMatching = Shape4DState::GetDataMatchingValues();
                vector<double> regularity = Shape4DState::GetRegularityValues();
                vector<double> stepsizes = Shape4DState::GetStepsizeValues();

                for (unsigned int iters=0; iters<dataMatching.size(); iters++)
                {
                    printf("Iteration %3d   funct = %0.4f   data = %0.4f   reg = %0.4f   step = %0.10lf\n",
                           iters+1, dataMatching[iters] + source.GetGamma()*regularity[iters], dataMatching[iters], regularity[iters], stepsizes[iters]);
                }

            }
            else
            {
                Shape4DState::ResetOptimizationState();
                printf("Starting acceleration controlled shape regression...\n");
                printf("SigmaV = %0.4f   GammaR = %0.4f   SigmaW = %0.4f\n", sigmaV, gammaR, ((TargetData*)targetArray[0])->GetSigmaW());
            }

            regression = (Regression*) new RegressionAcceleration(source, numTargets, targetArray);
            //regression->SetTheSource(theSource);
            //regression->SetTargets(theTarget);

            // We set the impulse to continue where we left off
            if (_continueExp)
            {
                ((RegressionAcceleration*)regression)->SetImpulse(*(_impulse));
            }

            time_t startTime = time(NULL);
            // Start the accel based regression
            X = regression->Run();
            time_t endTime = time(NULL);
            double regressionTime = difftime(endTime, startTime);

            dX = ((RegressionAcceleration*)regression)->GetVelocity();
            accel = ((RegressionAcceleration*)regression)->GetAcceleration();
            impulse = ((RegressionAcceleration*)regression)->GetImpulse();

            printf("\nElapsed time is %0.2f seconds\n\n", regressionTime);
        }
        else
        {
            printf("***ERROR*** No algorithm named \"%s\".\n", algorithmType);
            exit(1);
        }

        // Clean up memory
        delete regression;
        delete [] targetArray;
        delete [] sourcePathArray;
        delete [] outputPathArray;
        delete [] numSourcePtsArray;
        delete [] numSourceTrisArray;
    }
}

//----------------------------------------------------------------
// ContinueExperiment
//----------------------------------------------------------------
// Inputs: 
//   pathToFile - path to .exo file
//
// Outputs: 
//----------------------------------------------------------------
// Continues an experiment defined in an .exo file
//----------------------------------------------------------------
void RunExperiment::ContinueExperiment(char* pathToFile)
{
    this->_continueExp = true;

    std::string line;
    char* charLine;

    std::ifstream exoFile(pathToFile);

    // Get the version
    char* version;
    getline(exoFile, line);
    version = new char[line.size()+1];
    version[line.size()] = 0;
    memcpy(version, line.c_str(), line.size());

    // Get the experiment file
    char* xmlFile;
    getline(exoFile, line);								// EXPERIMENT FILE
    getline(exoFile, line);
    xmlFile = new char[line.size()+1];
    xmlFile[line.size()] = 0;
    memcpy(xmlFile, line.c_str(), line.size());

    // Read number of iterations
    getline(exoFile, line);								// COMPLETED ITERATIONS
    getline(exoFile, line);
    charLine = new char[line.size()+1];
    charLine[line.size()] = 0;
    memcpy(charLine, line.c_str(), line.size());
    int iterations;
    sscanf(charLine,"%d", &iterations);
    Shape4DState::SetIterations(iterations);

    vector<double> dataMatching(iterations);
    vector<double> regularity(iterations);
    vector<double> stepsize(iterations);

    // Get iteration information
    getline(exoFile, line);								// MATCHING REGULARITY STEPSIZE
    for (int i=0;i<iterations; i++)
    {
        getline(exoFile, line);
        charLine = new char[line.size()+1];
        charLine[line.size()] = 0;
        memcpy(charLine, line.c_str(), line.size());
        double d,r,s;
        sscanf(charLine,"%lf %lf %lf", &d, &r, &s);

        dataMatching[i] = d;
        regularity[i] = r;
        stepsize[i] = s;
    }

    Shape4DState::SetDataMatchingValues(dataMatching);
    Shape4DState::SetRegularityValues(regularity);
    Shape4DState::SetStepsizeValues(stepsize);

    // Read the vector fields
    getline(exoFile, line);
    charLine = new char[line.size()+1];
    charLine[line.size()] = 0;
    memcpy(charLine, line.c_str(), line.size());

    // This is a regression accel experiment
    if (strcmp(charLine, "INIT VELOCITY") == 0)
    {
        this->_skipRegressionVelocity = true;

        getline(exoFile, line);								// INIT VELOCITY
        delete [] charLine;
        charLine = new char[line.size()+1];
        charLine[line.size()] = 0;
        memcpy(charLine, line.c_str(), line.size());
        int dim, numPts;
        sscanf(charLine,"%d %d", &dim, &numPts);

        _initV0 = new Array2D<double>(dim, numPts);

        for (int i=0; i<numPts; i++)
        {
            getline(exoFile, line);
            delete [] charLine;
            charLine = new char[line.size()+1];
            charLine[line.size()] = 0;
            memcpy(charLine, line.c_str(), line.size());
            double x, y, z;
            sscanf(charLine,"%lf %lf %lf", &x, &y, &z);

            _initV0->SetAt(0,i,x);
            _initV0->SetAt(1,i,y);
            _initV0->SetAt(2,i,z);
        }

        // Read the impulse vectors
        getline(exoFile, line);								// IMPULSE
        getline(exoFile, line);
        delete [] charLine;
        charLine = new char[line.size()+1];
        charLine[line.size()] = 0;
        memcpy(charLine, line.c_str(), line.size());
        int T;
        sscanf(charLine,"%d %d %d", &dim, &numPts, &T);

        _impulse = new Array3D<double>(dim, numPts, T);

        for (int i=0; i<T; i++)
        {
            for (int j=0; j<numPts; j++)
            {
                getline(exoFile, line);
                delete [] charLine;
                charLine = new char[line.size()+1];
                charLine[line.size()] = 0;
                memcpy(charLine, line.c_str(), line.size());
                double x, y, z;
                sscanf(charLine,"%lf %lf %lf", &x, &y, &z);

                _impulse->SetAt(0,j,i,x);
                _impulse->SetAt(1,j,i,y);
                _impulse->SetAt(2,j,i,z);
            }
        }
    }
    // This is a regression velocity experiment
    else
    {
        // Read the momenta vectors
        getline(exoFile, line);								// MOMENTA
        delete [] charLine;
        charLine = new char[line.size()+1];
        charLine[line.size()] = 0;
        memcpy(charLine, line.c_str(), line.size());
        int T, dim, numPts;
        sscanf(charLine,"%d %d %d", &dim, &numPts, &T);

        _momenta = new Array3D<double>(dim, numPts, T);

        for (int i=0; i<T; i++)
        {
            for (int j=0; j<numPts; j++)
            {
                getline(exoFile, line);
                delete [] charLine;
                charLine = new char[line.size()+1];
                charLine[line.size()] = 0;
                memcpy(charLine, line.c_str(), line.size());
                double x, y, z;
                sscanf(charLine,"%lf %lf %lf", &x, &y, &z);

                _momenta->SetAt(0,j,i,x);
                _momenta->SetAt(1,j,i,y);
                _momenta->SetAt(2,j,i,z);
            }
        }
    }

    // We have everything we need to restart the experiment
    this->StartExperiment(xmlFile);

    exoFile.close();

    delete version;
    delete xmlFile;
}
