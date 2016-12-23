#include "InputFunctions.h"
using namespace std;

double FunctionF(PointDouble p) {
	double x = p.first;
	double y = p.second;
	return 2*(x*x + y*y)*(1 - 2*x*x*y*y)*exp(1 - x*x*y*y);
}

double FunctionPhi(PointDouble p) {
	double x = p.first;
	double y = p.second;
	return exp(1 - x*x*y*y);
}

Args ParseArgsStr(int argc, char** argv)
{
	int opt;
	char *opt_val = NULL;
	Args args;
	char optstring[15] = "x:y:a:b:r:c:f:";
	opt = ParsingOptions::getopt(argc, argv, optstring);
	while (opt != -1) {
		switch (opt) {
		case 'x':
			opt_val = ParsingOptions::optarg;
			args.leftBottomCornerX = strtol(opt_val, NULL, 10);
			break;
		case 'y':
			opt_val = ParsingOptions::optarg;
			args.leftBottomCornerY = strtol(opt_val, NULL, 10);
			break;
		case 'a':
			opt_val = ParsingOptions::optarg;
			args.rightTopCornerX = strtol(opt_val, NULL, 10);
			break;
		case 'b':
			opt_val = ParsingOptions::optarg;
			args.rightTopCornerY = strtol(opt_val, NULL, 10);
			break;
		case 'r':
			opt_val = ParsingOptions::optarg;
			args.numRows = strtol(opt_val, NULL, 10);
			break;
		case 'c':
			opt_val = ParsingOptions::optarg;
			args.numCols = strtol(opt_val, NULL, 10);
			break;
		case 'f':
			opt_val = ParsingOptions::optarg;
			args.resultGridOutputFname = string(opt_val);
			break;
		}
		opt = ParsingOptions::getopt(argc, argv, "x:y:a:b:r:c:f:");
	}
	return args;
}