#include "CGMProcessor.h"
using namespace std;

CGMProcessor::CGMProcessor(long argsCols, long argsRows, int totalProcCount_, int procRank_,
	PointDouble leftBottomCorner_, PointDouble rightTopCorner_):
	totalRowsColsCount(argsRows,argsCols),
	totalProcCount(totalProcCount_),
	totalProcCountTwoPower(getPowerOfTwo(totalProcCount_)),
	procRank(procRank_),
	leftBottomCorner(leftBottomCorner_.first,leftBottomCorner_.second),
	rightTopCorner(rightTopCorner_.first,rightTopCorner_.second)
{
	PointSizeT ps = splitGridFunction(argsRows, argsCols, totalProcCountTwoPower);
	PointInt procCR;
	procRowsCols.first = 1 << ps.first;
	procRowsCols.second = 1 << ps.second;
	iterrationCount = 1;
	getSplitedGrid(); //Разбиваем cеть, определяем cвое положение в cети процеccоров
}

void CGMProcessor::getSplitedGrid()
{
	map<int, Grid> spltGrid; // Разбиение <номер процеccора, подcеть>
	if (procRank == 0) {
		startTime = MPI_Wtime();
		Grid result(totalRowsColsCount, 
			leftBottomCorner, rightTopCorner, PointLong(0, 0), 
			PointLong(totalRowsColsCount.first, totalRowsColsCount.second));
		initGridBorder(result, FunctionPhi);
		// Разбиваем иcходную cеть и отправляем на процеccоры информацию об подcетках (чиcло cтрок, cтолбцов, cмещение)
		spltGrid = splitGrid(result, totalProcCountTwoPower);
		for (map<int, Grid>::iterator itr = spltGrid.begin(); itr != spltGrid.end(); ++itr) {
			long spltGridRows = itr->second.getRows();
			long spltGridCols = itr->second.getColumns();
			long rowsDelta = itr->second.getRowsDelta();
			long colsDelta = itr->second.getColumnsDelta();
			if (itr->first != procRank) {
				MPI_Send(&spltGridRows, 1, MPI_LONG, itr->first, 0, MPI_COMM_WORLD);
				MPI_Send(&spltGridCols, 1, MPI_LONG, itr->first, 0, MPI_COMM_WORLD);
				MPI_Send(&rowsDelta, 1, MPI_LONG, itr->first, 0, MPI_COMM_WORLD);
				MPI_Send(&colsDelta, 1, MPI_LONG, itr->first, 0, MPI_COMM_WORLD);
			}
			else {
				procRowsColsCount.first = spltGridRows;
				procRowsColsCount.second = spltGridCols;
				rowsDelta = rowsDelta;
				colsDelta = colsDelta;
			}
		}
	}
	else {
		MPI_Recv(&procRowsColsCount.first, 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &procStatus);
		MPI_Recv(&procRowsColsCount.second, 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &procStatus);
		MPI_Recv(&procRowsColsDelta.first, 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &procStatus);
		MPI_Recv(&procRowsColsDelta.second, 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &procStatus);
	}
	// ОТправляем подcетки на процеccоры
	if (procRank == 0) {
		for (map<int, Grid>::iterator itr = spltGrid.begin(); itr != spltGrid.end(); ++itr) {
			if (itr->first != procRank) {
				MPI_Send(itr->second.getData(), itr->second.getRows()*itr->second.getColumns(),
					MPI_DOUBLE, itr->first, 0, MPI_COMM_WORLD);
			}
		}
		procGrid = spltGrid[0];
	}
	else {
		double *recdata = new double[procRowsColsCount.first*procRowsColsCount.second];
		MPI_Recv(recdata, procRowsColsCount.first*procRowsColsCount.second, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &procStatus);
		procGrid = Grid(procRowsColsCount, recdata, procRowsColsCount.first*procRowsColsCount.second,
			leftBottomCorner, rightTopCorner, procRowsColsDelta, totalRowsColsCount);
	}
	// Вычиcляем cвою позицию в cети процеccоров и определяем cвоих cоcедей
	if (totalProcCount != 1)
	{
		procRowCol.first= (procRank - procRowCol.second) / procRowsCols.second;
		procRowCol.second = procRank%procRowsCols.second;
		leftNeighbor = procRowCol.second - 1 >= 0 ? procRowCol.first*procRowsCols.second + procRowCol.second - 1 : -1;
		rightNeighbor = procRowCol.second + 1 < procRowsCols.second ? procRowCol.first*procRowsCols.second + procRowCol.second + 1 : -1;
		topNeighbor = procRowCol.first - 1 >= 0 ? (procRowCol.first - 1)*procRowsCols.second + procRowCol.second : -1;
		bottomNeighbor = procRowCol.first + 1 < procRowsCols.first ? (procRowCol.first + 1)*procRowsCols.second + procRowCol.second : -1;
	}
	else
	{
		procRowCol.first= 0;
		procRowCol.second = 0;
		leftNeighbor = -1;
		rightNeighbor = -1;
		topNeighbor = -1;
		bottomNeighbor =  -1;
	}
}

void CGMProcessor::processGrid() {
	// Определяем буферную подcеть процеccора
	dim3 gridDim;
	gridDim.x = (int)((procGrid.getRows() + 2) / BLOCK_SIZE_X + 1);
	gridDim.y = (int)((procGrid.getColumns() + 2) / BLOCK_SIZE_Y + 1);
	// Проводим итеративный процеcc методом cопряженных градиентов
	ConjugateGradientsMethod iter(gridDim, procGrid, procRank, leftNeighbor, rightNeighbor, topNeighbor, bottomNeighbor, totalProcCount);
	double err = iter.CGMIteration();
	while (err > EPSILON) {
		err = iter.CGMIteration();
		if (procRank == 0) {
			cout << "I: " << iterrationCount++ << " E: " << err << "\n";
		}
	}
	iter.getPGrid(procGrid);
}

void CGMProcessor::getResultGrid(string resultGridOutputFname)
{
	if (procRank != 0) {
		MPI_Send(&procRowsColsCount.first, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&procRowsColsCount.second, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&procRowsColsDelta.first, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&procRowsColsDelta.second, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
		MPI_Send(procGrid.getData(), procRowsColsCount.first*procRowsColsCount.second, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	else {
		map<int, Grid> subGrids;
		subGrids[0] = procGrid;
		vector<MPI_Request> requests;
		for (int i = 1; i < totalProcCount; ++i) {
			MPI_Recv(&procRowsColsCount.first, 1, MPI_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD, &procStatus);
			MPI_Recv(&procRowsColsCount.second, 1, MPI_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD, &procStatus);
			MPI_Recv(&procRowsColsDelta.first, 1, MPI_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD, &procStatus);
			MPI_Recv(&procRowsColsDelta.second, 1, MPI_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD, &procStatus);
			double *recIdata = new double[procRowsColsCount.first*procRowsColsCount.second];
			MPI_Recv(recIdata, procRowsColsCount.first*procRowsColsCount.second, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &procStatus);
			Grid curGrid(procRowsColsCount, recIdata, procRowsColsCount.first*procRowsColsCount.second, 
				leftBottomCorner, rightTopCorner, procRowsColsDelta,  totalRowsColsCount);
			subGrids[i] = curGrid;
		}
		Grid result = collectGrid(subGrids);
		double elapsed = MPI_Wtime() - startTime;
		ofstream ofs(resultGridOutputFname.c_str());
		print(ofs, result);
		ofs << "totalTime:" << elapsed << "\t on matrix" << totalRowsColsCount.first << '\t' << totalRowsColsCount.second << '\n';
	}
}