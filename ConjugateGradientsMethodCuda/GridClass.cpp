#include "GridClass.h"

const double Grid::Q_COEF = 1.5;
const double Grid::Q_COEF_TWOPOW = 1.82842712474;
Grid::Grid(PointLong subGridRowsCols, 
	PointSizeT leftBottomCorner_, 
	PointSizeT rightTopCorner_, 
	PointLong rowsColsDelta_, 
	PointLong totalRowsColsCount_):
leftBottomCorner(leftBottomCorner_),
	rightTopCorner(rightTopCorner_),
	rowsColsDelta(rowsColsDelta_),
	totalRowsColsCount(totalRowsColsCount_),
	data(subGridRowsCols.first, subGridRowsCols.second, 0)
{
	initCacheForGridPoints();
}

Grid::Grid(PointLong subGridRowsCols, 
	double *data, 
	long dataCount, 
	PointSizeT leftBottomCorner_, 
	PointSizeT rightTopCorner_, 
	PointLong rowsColsDelta_, 
	PointLong totalRowsColsCount_):
leftBottomCorner(leftBottomCorner_),
	rightTopCorner(rightTopCorner_),
	rowsColsDelta(rowsColsDelta_),
	totalRowsColsCount(totalRowsColsCount_),
	data(data, dataCount, subGridRowsCols.first, subGridRowsCols.second) {
		initCacheForGridPoints();
}

Grid:: Grid(const Matrix &data, 
	PointSizeT leftBottomCorner_, 
	PointSizeT rightTopCorner_, 
	PointLong rowsColsDelta_, 
	PointLong totalRowsColsCount_):
leftBottomCorner(leftBottomCorner_),
	rightTopCorner(rightTopCorner_),
	rowsColsDelta(rowsColsDelta_),
	totalRowsColsCount(totalRowsColsCount_),
	data(data) {
		initCacheForGridPoints();
}

void initGridBorder(Grid &grid, Function f) {
	for (size_t i = 0; i < grid.getRows() ; ++i ){
		grid(i,0) = f(grid.getPoint(i,0));
		grid(i, grid.getColumns() - 1) = f(grid.getPoint(i, grid.getColumns() - 1));
	}
	for (size_t j = 0; j < grid.getColumns(); ++j) { grid(0,j) = f(grid.getPoint(0,j));
	grid(grid.getRows() - 1,j) = f(grid.getPoint(grid.getRows() - 1, j));
	}
}

map<int, Grid> splitGrid(const Grid &grid, long procCountTwoPower) {
	PointSizeT size = splitGridFunction(grid.getRows(), grid.getColumns(), procCountTwoPower);
	long rows = 1 << size.first;
	long cols = 1 << size.second;
	map<int, Grid> result;
	map<pair<int,int>, Matrix> submats = split(grid.data, rows, cols);
	int procCount = 0;
	long deltaRows = -1, deltaCols = -1;
	for(map<pair<int,int>,Matrix>::iterator it = submats.begin(); it != submats.end(); ++it) {
		int r = it->first.first;
		int c = it->first.second;
		if(deltaRows == -1){
			deltaRows = it->second.rowsCount();
			deltaCols = it->second.colsCount();
		}
		result[procCount++] = Grid(it->second, grid.leftBottomCorner, 
			grid.rightTopCorner, PointLong(r*deltaRows, c*deltaCols), PointLong(grid.getRows(), grid.getColumns()));
	}
	return result;
}

Grid collectGrid(const map<int, Grid> &subGrids) {
	Grid first = subGrids.at(0);
	Grid result(PointLong(first.getTotalRowsCount(), first.getTotalColsCount()), 
		first.getLeftBottomCorner(), first.getRightTopCorner(), PointLong(0, 0), 
		PointLong(first.getTotalRowsCount(), first.getTotalColsCount()));
	for(map<int, Grid>::const_iterator itr = subGrids.begin(); itr!=subGrids.end(); ++itr){
		for(long i = 0; i < itr->second.getRows(); ++i) {
			for(long j = 0; j < itr->second.getColumns(); ++j) {
				result(i+itr->second.getRowsDelta(), j + itr->second.getColumnsDelta()) = itr->second(i,j);
			}
		}
	}
	return result;
}

int getPowerOfTwo(int val){
	int pwr = 0;
	while(val >>= 1) ++pwr;
	return pwr;
}

PointSizeT splitGridFunction(int rows, int cols, int sizePower) {
	double rows_ = (double)rows;
	double cols_ = (double)cols;
	int sizePower_ = 0;

	for (int i = 0; i < sizePower; i++) {
		if (rows_ > cols_)
		{
			rows_ = rows_ / 2.0;
			++sizePower_;
		}
		else {
			cols_ = cols_ / 2.0;
		}
	}
	return PointSizeT(sizePower_, sizePower - sizePower_);
}
