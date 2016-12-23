#include "ConjugateGradientsMethod.cuh"
using namespace std;

const long BLOCK_SIZE_X = 16;
const long BLOCK_SIZE_Y = 32;
const long BLOCK_SIZE = BLOCK_SIZE_X*BLOCK_SIZE_Y;
const dim3 BLOCK_DIM(BLOCK_SIZE_X, BLOCK_SIZE_Y);

__device__ double CUDAFunctionF(double x, double y){
	return 2*(x*x + y*y)*(1 - 2*x*x*y*y)*exp(1 - x*x*y*y);
}
__device__ double fiveDotScheme(CudaGrid m, long i,long j) {
	double prevX = m.getGridSplitOx(i-1);
	double prevY = m.getGridSplitOy(j-1);
	double curX = m.getGridSplitOx(i);
	double curY = m.getGridSplitOy(j);
	double nextX = m.getGridSplitOx(i+1);
	double nextY = m.getGridSplitOy(j+1);
	double avgStepOx = (nextX - prevX) / 2;
	double avgStepOy = (nextY - prevY) / 2;
	double leftPoint = m(i, j-1);
	double rightPoint = m(i, j+1);
	double bottomPoint = m(i+1, j);
	double topPoint = m(i-1, j);
	double x = ((m(i,j) - topPoint)/(curX - prevX) - (bottomPoint - m(i,j))/(nextX - curX))/avgStepOx;
	double y = ((m(i,j) - leftPoint)/(curY - prevY) - (rightPoint - m(i,j))/(nextY- curY))/avgStepOy;	
	return x + y;
}
void ConjugateGradientsMethod::getGridBorders(CudaGrid &grid) {
	MPI_Status status[4];
	MPI_Request send[4];
	MPI_Request recv[4];
	vector<double> topBuf;
	vector<double> bottomBuf;
	vector<double> rightBuf;
	vector<double> leftBuf;
	if (topNeighbor >= 0 && topNeighbor < totalProcCount) {
		topBuf.resize(grid.colsCount());
		vector<double> topCur;
		topCur.resize(grid.colsCount());
		grid.getTop(topCur);
		MPI_Isend(&topCur[0], topCur.size(), MPI_DOUBLE, topNeighbor, 0, MPI_COMM_WORLD, &send[0]);
		MPI_Irecv(&topBuf[0], grid.colsCount(), MPI_DOUBLE, topNeighbor, MPI_ANY_TAG, MPI_COMM_WORLD, &recv[0]);

		MPI_Wait(&send[0],&status[0]);
		MPI_Wait(&recv[0],&status[0]);
	}
	if(bottomNeighbor >= 0 && bottomNeighbor < totalProcCount) {
		bottomBuf.resize(grid.colsCount());
		vector<double> bottomCur;
		bottomCur.resize(grid.colsCount());
		grid.getBottom(bottomCur);
		MPI_Isend(&bottomCur[0], bottomCur.size(), MPI_DOUBLE, bottomNeighbor, 0, MPI_COMM_WORLD, &send[1]);
		MPI_Irecv(&bottomBuf[0], grid.colsCount(), MPI_DOUBLE, bottomNeighbor, MPI_ANY_TAG, MPI_COMM_WORLD, &recv[1]);

		MPI_Wait(&send[1],&status[1]);
		MPI_Wait(&recv[1],&status[1]);
	}
	if (rightNeighbor >= 0 && rightNeighbor < totalProcCount) {
		rightBuf.resize(grid.rowsCount());
		vector<double> rightCur;
		rightCur.resize(grid.rowsCount());
		grid.getRight(rightCur);
		MPI_Isend(&rightCur[0], rightCur.size(), MPI_DOUBLE, rightNeighbor, 0, MPI_COMM_WORLD, &send[2]);
		MPI_Irecv(&rightBuf[0], grid.rowsCount(), MPI_DOUBLE, rightNeighbor, MPI_ANY_TAG, MPI_COMM_WORLD, &recv[2]);
		MPI_Wait(&send[2],&status[2]);
		MPI_Wait(&recv[2],&status[2]);
	}
	if(leftNeighbor >= 0 && leftNeighbor < totalProcCount) {
		leftBuf.resize(grid.rowsCount());
		vector<double> leftCur;
		leftCur.resize(grid.rowsCount());
		grid.getLeft(leftCur);
		MPI_Isend(&leftCur[0], leftCur.size(), MPI_DOUBLE, leftNeighbor, 0, MPI_COMM_WORLD, &send[3]);
		MPI_Irecv(&leftBuf[0], grid.rowsCount(), MPI_DOUBLE, leftNeighbor, MPI_ANY_TAG, MPI_COMM_WORLD, &recv[3]);
		MPI_Wait(&send[3],&status[3]);
		MPI_Wait(&recv[3],&status[3]);
	}
	if (!topBuf.empty()) {
		grid.setTopBorder(topBuf);
	}
	if (!bottomBuf.empty()) {
		grid.setBottomBorder(bottomBuf);
	}
	if (!rightBuf.empty()) {
		grid.setRightBorder(rightBuf);
	}
	if (!leftBuf.empty()) {
		grid.setLeftBorder(leftBuf);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}


// Ядро cбора значений cкалярных произвдений
template <unsigned int blockSize>
__device__ void CUDA_Allreduce(DividendDivider*sData, DividendDivider *global, long tid, double nom, double denom) {
	sData[tid].first = nom;
	sData[tid].second = denom;
	__syncthreads();
	if (blockSize >= 512) {
		if (tid < 256) {
			sData[tid].first += sData[tid + 256].first;
			sData[tid].second += sData[tid + 256].second;
		}
		__syncthreads();
	}
	if (blockSize >= 256) {
		if (tid < 128) {
			sData[tid].first += sData[tid + 128].first;
			sData[tid].second += sData[tid + 128].second;
		}
		__syncthreads();
	}
	if (blockSize >= 128) {
		if (tid < 64) {
			sData[tid].first += sData[tid + 64].first;
			sData[tid].second += sData[tid + 64].second;
		}
		__syncthreads();
	}
	if (tid < 32) {
		if (blockSize >= 64) {
			sData[tid].first += sData[tid + 32].first;
			sData[tid].second += sData[tid + 32].second;
		}
		__syncthreads();
		if (blockSize >= 32) {
			sData[tid].first += sData[tid + 16].first;
			sData[tid].second += sData[tid + 16].second;
		}
		__syncthreads();
		if (blockSize >= 16) {
			sData[tid].first += sData[tid + 8].first;
			sData[tid].second += sData[tid + 8].second;
		}
		__syncthreads();
		if (blockSize >= 8) {
			sData[tid].first += sData[tid + 4].first;
			sData[tid].second += sData[tid + 4].second;
		}
		__syncthreads();
		if (blockSize >= 4) {
			sData[tid].first += sData[tid + 2].first;
			sData[tid].second += sData[tid + 2].second;
		}
		__syncthreads();
		if (blockSize >= 2) {
			sData[tid].first += sData[tid + 1].first;
			sData[tid].second += sData[tid + 1].second;
		}
		__syncthreads();
	}
	if (tid == 0) {
		long blockIndex = gridDim.x*blockIdx.y + blockIdx.x;
		global[blockIndex].first = sData[0].first;
		global[blockIndex].second = sData[0].second;
	}
}
// Ядра вычиcления матриц
__global__ void rGridCalculationKernel(CudaGrid rGrid, CudaGrid pGrid) {
	long i = blockDim.x*blockIdx.x + threadIdx.x + 1;
	long j = blockDim.y*blockIdx.y + threadIdx.y + 1;
	if (i < rGrid.rowsCount() - 1 && i > 0 && j < rGrid.colsCount() - 1 && j > 0) {
			rGrid(i,j) = fiveDotScheme(pGrid, i, j) - CUDAFunctionF(pGrid.getGridSplitOx(i), pGrid.getGridSplitOy(j));
	}
}
void ConjugateGradientsMethod::rGridMatrixCalculation(CudaGrid rGrid, CudaGrid pGrid) {
	rGridCalculationKernel<<<gridDim, BLOCK_DIM>>>(rGrid, pGrid);
}

__global__ void gGridCalculationKernel(CudaGrid gGrid, CudaGrid rGrid, double alpha ) {
	long i = blockDim.x*blockIdx.x + threadIdx.x + 1;
	long j = blockDim.y*blockIdx.y + threadIdx.y + 1;
	if (i < rGrid.rowsCount() - 1 && i > 0 && j < rGrid.colsCount() - 1 && j > 0) {
			gGrid(i, j) = rGrid(i, j) - alpha*gGrid(i,j);
	}
}
void ConjugateGradientsMethod::gGridMatrixCalculation(CudaGrid gGrid, CudaGrid rGrid, double alpha) {
	gGridCalculationKernel<<<gridDim, BLOCK_DIM>>>(gGrid, rGrid, alpha);
}

__global__ void pGridCalculationKernel(CudaGrid pGrid, CudaGrid gGrid, double tauCoef, DividendDivider *erros) {
	extern __shared__ DividendDivider shared[];
	long i = blockDim.x*blockIdx.x + threadIdx.x + 1;
	long j = blockDim.y*blockIdx.y + threadIdx.y + 1;
	long tid = threadIdx.y*BLOCK_SIZE_X + threadIdx.x;
	double err = 0;
	if (i < pGrid.rowsCount() - 1 && i > 0 && j < pGrid.colsCount() - 1 && j > 0) {
			double val = pGrid(i,j) - tauCoef*gGrid(i,j);
			double valueError = pGrid(i,j) - val;
			double avgStepOx = pGrid.getAvgStepOx(i);
			double avgStepOy = pGrid.getAvgStepOy(j);
			err = valueError*valueError*avgStepOx*avgStepOy;
			pGrid(i, j) = val;
	}
	CUDA_Allreduce<BLOCK_SIZE>(shared, erros, tid, err, 0);
}
double ConjugateGradientsMethod::pGridMatrixCalculation(CudaGrid pGrid, CudaGrid gGrid) {
	DividendDivider*fractions;
	int errAlloc=cudaMalloc(&fractions, gridDim.x*gridDim.y*sizeof(DividendDivider));
	if (errAlloc!=cudaSuccess) 
		cerr<<"CUDA ERROR "<<errAlloc<<" in cudaMalloc at Function pGridMatrixCalculation\n";
	pGridCalculationKernel<<<gridDim, BLOCK_DIM, BLOCK_SIZE*sizeof(DividendDivider)>>>(pGrid, gGrid, tauCoef, fractions);
	return getScalarProductFromCuda(fractions, gridDim.x*gridDim.y).first;
}

// Раccчет коэфициентов тау и альфа
__global__ void kernelCalcTau(CudaGrid gGrid, CudaGrid rGrid, DividendDivider *tauCoefs) {
	extern __shared__ DividendDivider shared[];
	long i = blockDim.x*blockIdx.x + threadIdx.x + 1;
	long j = blockDim.y*blockIdx.y + threadIdx.y + 1;
	long tid = threadIdx.y*BLOCK_SIZE_X + threadIdx.x;
	double tDividend = 0, tDivider = 0;
	if (i < rGrid.rowsCount() - 1 && i > 0 &&
		j < rGrid.colsCount() - 1 && j > 0) {
			tDividend = gGrid(i,j)*rGrid(i,j);
			tDivider = fiveDotScheme(gGrid, i, j)*gGrid(i, j);
	}
	CUDA_Allreduce<BLOCK_SIZE>(shared, tauCoefs, tid, tDividend, tDivider);
}
DividendDivider ConjugateGradientsMethod::tauCoefCalculation(CudaGrid gGrid, CudaGrid rGrid) {
	DividendDivider *fractions;
	int errAlloc=cudaMalloc(&fractions, gridDim.x*gridDim.y*sizeof(DividendDivider));
	if (errAlloc!=cudaSuccess) 
		cerr<<"CUDA ERROR "<<errAlloc<<" in cudaMalloc at Function tauCoefCalculation\n";
	kernelCalcTau<<<gridDim, BLOCK_DIM, BLOCK_SIZE*sizeof(DividendDivider)>>>(gGrid, rGrid, fractions);
	return getScalarProductFromCuda(fractions, gridDim.x*gridDim.y);
}

__global__ void kernelCalcAlpha(CudaGrid rGrid, CudaGrid gGrid, DividendDivider *alphas) {
	extern __shared__ DividendDivider shared[];
	long i = blockDim.x*blockIdx.x + threadIdx.x +1 ;
	long j = blockDim.y*blockIdx.y + threadIdx.y +1;
	double aDividend = 0, aDivider = 0;
	long tid = threadIdx.y*BLOCK_SIZE_X + threadIdx.x;
	if (i < rGrid.rowsCount() - 1 && i > 0 &&
		j < rGrid.colsCount() - 1 && j > 0) {
			aDividend = fiveDotScheme(rGrid, i, j)*gGrid(i, j);
			aDivider = fiveDotScheme(gGrid, i, j)*gGrid(i, j);
	}
	CUDA_Allreduce<BLOCK_SIZE>(shared, alphas, tid, aDividend, aDivider);
}
DividendDivider ConjugateGradientsMethod::alphaCoefCalculation(CudaGrid rGrid, CudaGrid gGrid) {
	DividendDivider *fractions;
	int errAlloc=cudaMalloc(&fractions, gridDim.x*gridDim.y*sizeof(DividendDivider));
	if (errAlloc!=cudaSuccess) 
		cerr<<"CUDA ERROR "<<errAlloc<<" in cudaMalloc at Function alphaCoefCalculation\n";
	kernelCalcAlpha<<<gridDim, BLOCK_DIM, BLOCK_SIZE*sizeof(DividendDivider)>>>(rGrid, gGrid, fractions);
	return getScalarProductFromCuda(fractions, gridDim.x*gridDim.y);
}
// Получить результаты cкалярных умножений c куда
DividendDivider ConjugateGradientsMethod::getScalarProductFromCuda(DividendDivider *numbers, long totalProcCount) const {
	vector<DividendDivider> localData;
	localData.resize(totalProcCount);
	int errCpy=cudaMemcpy(&localData[0], numbers, totalProcCount*sizeof(DividendDivider), cudaMemcpyDeviceToHost);
	if (errCpy!=cudaSuccess) 
		cerr<<"CUDA ERROR "<<errCpy<<" in cudaMemcpy at Function getScalarProductFromCuda\n";
	DividendDivider result;
	result.first = 0;
	result.second = 0;
	for (long i = 0; i < totalProcCount; ++i){
		result.first += localData[i].first;
		result.second += localData[i].second;
	}
	return result;
}

// Итерация Метода cопряженных градиентов
double ConjugateGradientsMethod::CGMIteration() {
	double error;
	if(iterationNum == 0) {
		getGridBorders(pGrid);
		rGridMatrixCalculation(rGrid, pGrid);
		rGridMatrixCalculation(gGrid, pGrid);
		getGridBorders(rGrid);
		getGridBorders(gGrid);
		DividendDivider tFr = tauCoefCalculation(rGrid, rGrid);
		double allNumerator, allDeDividend;
		MPI_Allreduce(&tFr.first, &allNumerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&tFr.second, &allDeDividend, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		tauCoef = allNumerator/allDeDividend;
		error = 100000;
	} else {
		double localErr = pGridMatrixCalculation(pGrid, gGrid);
		MPI_Allreduce(&localErr, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		getGridBorders(pGrid);
		rGridMatrixCalculation(rGrid, pGrid);

		getGridBorders(rGrid);
		DividendDivider aFr = alphaCoefCalculation(rGrid, gGrid);
		double allADividend, allADivider;
		MPI_Allreduce(&aFr.first, &allADividend, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&aFr.second, &allADivider, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		double alpha = allADividend / allADivider;
		gGridMatrixCalculation(gGrid, rGrid, alpha);

		getGridBorders(gGrid);
		DividendDivider tFr = tauCoefCalculation(gGrid, rGrid);
		double allNumerator, allDeDividend;
		MPI_Allreduce(&tFr.first, &allNumerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&tFr.second, &allDeDividend, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		tauCoef = allNumerator/allDeDividend;
	}
	iterationNum++;
	return sqrt(error);
}

