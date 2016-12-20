#ifndef _INPUTFUNCS_H_
#define _INPUTFUNCS_H_
#include <cmath>
#include <string>
#include <cstdlib>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "CommonTypes.h"
#include "XGetopt.h" // Парcинг параметров командной cтроки

struct Args {
	long leftBottomCornerX;
	long leftBottomCornerY;
	long rightTopCornerX;
	long rightTopCornerY;
	long numRows;
	long numCols;
	string resultGridOutputFname;
};

const double EPSILON = 1E-4;
double FunctionF(PointDouble p);
double FunctionPhi(PointDouble p);
Args ParseArgsStr(int argc, char** argv);
#endif