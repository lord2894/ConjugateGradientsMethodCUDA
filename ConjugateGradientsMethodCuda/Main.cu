#include "CGMProcessor.h"
using namespace std;


int main(int argc, char **argv) {
    int totalprocCount, rank;
	//Args args = ParseArgsStr(argc, argv);
	//Теcтовый набор параметров--------
	Args args;
	args.leftBottomCornerX = 0;
	args.leftBottomCornerY = 0;
	args.numCols = 1000;
	args.numRows = 1000;
	args.resultGridOutputFname = "outputGrid.txt";
	args.rightTopCornerX = 2;
	args.rightTopCornerY = 2;
	//----------------------------------
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&totalprocCount);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    cudaSetDevice(0);
	CGMProcessor* processor = new CGMProcessor(args.numCols, args.numRows, 
		totalprocCount, rank, PointSizeT(args.leftBottomCornerX, 
			args.leftBottomCornerY), PointSizeT(args.rightTopCornerX, args.rightTopCornerY));
    MPI_Barrier(MPI_COMM_WORLD);
	processor->processGrid();
	MPI_Barrier(MPI_COMM_WORLD);
	processor->getResultGrid(args.resultGridOutputFname);
    MPI_Finalize();
    return 0;
}
