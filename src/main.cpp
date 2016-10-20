//--------------------------------------------------------------------------
// main.cpp:  Entry point for the application.
//--------------------------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>

#include "runexperiment.h"

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
int main(int argc, char **argv)
{	
	printf("\n");
		
	// Check command line arguments
	if (argc < 2)
	{
		printf("Usage:  exoshapeaccel driver_file.xml    OR\n        exoshapeaccel --continue progress_file.exo\n");
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
