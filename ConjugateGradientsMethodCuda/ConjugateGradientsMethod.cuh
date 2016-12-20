#ifndef _KERNELS_CUH_
#define _KERNELS_CUH_
#include "GridCLassCUDA.cuh"
#include "GridClass.h"
#include "InputFunctions.cuh"
#include <mpi.h>
typedef PointDouble DividendDivider;
extern const long BLOCK_SIZE_X;
extern const long BLOCK_SIZE_Y;
extern const long BLOCK_SIZE;
extern const dim3 BLOCK_DIM;

static inline int getColsAdd(int left, int right) {
	return (left!=-1) + (right!=-1);
}
static inline int getRowsAdd(int top, int bottom) {
	return (top!=-1) + (bottom!=-1);
}

class ConjugateGradientsMethod {
private:
	long leftNeighbor; // ����� �����cc���-c�c��� c����
	long rightNeighbor; // ����� �����cc���-c�c��� c�����
	long topNeighbor; // ����� �����cc���-c�c��� c����
	long bottomNeighbor; // ����� �����cc���-c�c��� c�����
	int procRank; // ����� �����cc���
	int totalProcCount; // ������c��� �����cc����
	int iterationNum; // ����� ��������
	// ������� �������� r, g, p c�����������, c����c�� ������ c���������� ����������
	CudaGrid rGrid;
	CudaGrid gGrid;
	CudaGrid pGrid;
	// ��cc��� c������c������ ������
	void gGridMatrixCalculation(CudaGrid gGrid, CudaGrid rGrid, double alpha);
	void rGridMatrixCalculation(CudaGrid rGrid, CudaGrid pGrid);
	double pGridMatrixCalculation(CudaGrid pGrid, CudaGrid gGrid);
	// ����������� ���
	double tauCoef;
	// ��cc��� �������������
	DividendDivider alphaCoefCalculation(CudaGrid rGrid, CudaGrid gGrid);
	DividendDivider tauCoefCalculation(CudaGrid gGrid, CudaGrid rGrid);
	dim3 gridDim; // ��������c�� c���� ��� CUDA
	// ����� c��������� c��c�� ����������� � ����c��� 0-� ��������
	double SteepestDescentMethodForZeroIter(CudaGrid &pGrid);
	double SteepestDescentMethodForZeroIter();
	// �������� ������� c����
	void getGridBorders(CudaGrid &grid);
	// C���������� c�������� ����������� ��� ��cc���� ���������� ��� � ����� � ����c�� ������ CUDA
	DividendDivider getScalarProductFromCuda(DividendDivider *numbers, long size) const;
public:
	ConjugateGradientsMethod(dim3 grid, const Grid& m, int procRank_, int leftNeighbor_, int rightNeighbor_, int topNeighbor_, int bottomNeighbor_, int totalProcCount_):
	gridDim(grid), procRank(procRank_),leftNeighbor(leftNeighbor_),
		rightNeighbor(rightNeighbor_),topNeighbor(topNeighbor_),
		bottomNeighbor(bottomNeighbor_),totalProcCount(totalProcCount_), iterationNum(0)
	{
		PointSizeT leftBottomCorner = m.getLeftBottomCorner(); 
		PointSizeT rightTopCorner = m.getRightTopCorner();
		PointLong rowsColsDelta(m.getRowsDelta(),m.getColumnsDelta());
		PointLong totalRowsColsCount(m.getTotalRowsCount(),m.getTotalColsCount());
		PointLong gridRowsColsCount(m.getRows(),m.getColumns());
		// �������������� ������� ���� ��������
		rGrid.setGridOptions(PointLong(gridRowsColsCount.first + getRowsAdd(topNeighbor,bottomNeighbor), gridRowsColsCount.second + getColsAdd(leftNeighbor, rightNeighbor)),
			leftBottomCorner,rightTopCorner,rowsColsDelta,totalRowsColsCount);
		pGrid.setGridOptions(PointLong(gridRowsColsCount.first + getRowsAdd(topNeighbor,bottomNeighbor), gridRowsColsCount.second + getColsAdd(leftNeighbor, rightNeighbor)),
			leftBottomCorner,rightTopCorner,rowsColsDelta,totalRowsColsCount);
		gGrid.setGridOptions(PointLong(gridRowsColsCount.first + getRowsAdd(topNeighbor,bottomNeighbor), gridRowsColsCount.second + getColsAdd(leftNeighbor, rightNeighbor)),
			leftBottomCorner,rightTopCorner,rowsColsDelta,totalRowsColsCount);
		// ���������� ��������� c���� �� �c�� Ox Oy
		initCacheForGridPoints(pGrid, leftNeighbor!=-1, topNeighbor !=-1);
		initCacheForGridPoints(rGrid, leftNeighbor!=-1, topNeighbor !=-1);
		initCacheForGridPoints(gGrid, leftNeighbor!=-1, topNeighbor !=-1);

		if (leftNeighbor == -1) {
			pGrid.setLeftBorder(m.getLeftColBorder(), topNeighbor!=-1);
		}
		if (rightNeighbor == -1) {
			pGrid.setRightBorder(m.getRightColBorder(), topNeighbor!=-1);
		}
		if (topNeighbor == -1) {
			pGrid.setTopBorder(m.getTopRowBorder(),leftNeighbor!=-1);
		}
		if (leftNeighbor == -1){
			pGrid.setBottomBorder(m.getBottomRowBorder(),leftNeighbor!=-1);
		}
	}
	double getCurrentIter() const { return iterationNum; }
	void getPGrid(Grid &grid) const {
		pGrid.tranformToSimpleGrid(PointInt(leftNeighbor!=-1, topNeighbor!=-1), grid);
	}
	double CGMIteration();	
};

#endif
