#include "GridCLassCUDA.cuh"

CudaGrid::CudaGrid():
leftBottomCorner(PointDouble(0,0)),
	rightTopCorner(PointDouble(0,0)),
	rowsColsDelta(0,0),
	totalRowsColsCount(0,0)
{
	data.xsize = 0;
	data.ysize = 0;
}

CudaGrid::CudaGrid(PointLong rowsColsCount_, PointSizeT leftBottomCorner_,
	PointSizeT rightTopCorner_, PointLong rowsColsDelta_, PointLong totalRowsColsCount_):
leftBottomCorner(leftBottomCorner_),
	rightTopCorner(rightTopCorner_),
	rowsColsDelta(rowsColsDelta_),
	totalRowsColsCount(totalRowsColsCount_)	
{
	data.xsize = rowsColsCount_.first;
	data.ysize = rowsColsCount_.second;
	int errAllocPitch=cudaMallocPitch(&data.ptr, &data.pitch, data.ysize*sizeof(double), data.xsize);
	if (errAllocPitch!=cudaSuccess) 
		cerr<<"CUDA ERROR "<<errAllocPitch<<" in cudaMallocPitch at Function CudaGrid\n";
	int errMemset2D = cudaMemset2D(data.ptr, data.pitch, 0, data.ysize*sizeof(double), data.xsize);
	if (errMemset2D!=cudaSuccess) 
		cerr<<"CUDA ERROR "<<errMemset2D<<" in cudaMemset2D at Function CudaGrid\n";
}

void CudaGrid::setGridOptions(PointLong rowsColsCount_, PointSizeT leftBottomCorner_,
	PointSizeT rightTopCorner_, PointLong rowsColsDelta_, PointLong totalRowsColsCount_)
{
	leftBottomCorner = leftBottomCorner_;
	rightTopCorner = rightTopCorner_;
	rowsColsDelta = rowsColsDelta_;
	totalRowsColsCount = totalRowsColsCount_,
	data.xsize = rowsColsCount_.first;
	data.ysize = rowsColsCount_.second;
	int errAllocPitch = cudaMallocPitch(&data.ptr, &data.pitch, data.ysize*sizeof(double), data.xsize);
	if (errAllocPitch!=cudaSuccess) 
		cerr<<"CUDA ERROR "<<errAllocPitch<<" in cudaMallocPitch at Function CudaGrid\n";
	int errMemset2D = cudaMemset2D(data.ptr, data.pitch, 0, data.ysize*sizeof(double), data.xsize);
	if (errMemset2D!=cudaSuccess) 
		cerr<<"CUDA ERROR "<<errMemset2D<<" in cudaMemset2D at Function CudaGrid\n";

}

void CudaGrid::setBorder(PointLong from, PointLong size, const double *ptr) {
	int errMemset2D = cudaMemcpy2D((char*)data.ptr + from.second*sizeof(double) + from.first*data.pitch,
		data.pitch, ptr, size.first*sizeof(double), size.first*sizeof(double), size.second, cudaMemcpyHostToDevice);
	if (errMemset2D!=cudaSuccess) 
		cerr<<"CUDA ERROR "<<errMemset2D<<" in cudaMemset2D at Function setBorder\n";
}

void CudaGrid::setLeftBorder(const vector<double> &left, bool delta) {
	setBorder(PointLong(delta, 0), PointLong(1, left.size()), &left[0]);
}
void CudaGrid::setRightBorder(const vector<double> &right, bool delta) {
	setBorder(PointLong(delta, colsCount()-1), PointLong(1, right.size()), &right[0]);
}
void CudaGrid::setTopBorder(const vector<double> &top, bool delta) {
	setBorder(PointLong(0, delta), PointLong(top.size(), 1), &top[0]);
}
void CudaGrid::setBottomBorder(const vector<double> &bottom, bool delta) {
	setBorder(PointLong(rowsCount() - 1, delta), PointLong(bottom.size(), 1), &bottom[0]);
}

void CudaGrid::getBorder(double *ptr, PointLong from, PointLong size) const {
	int errMemset2D = cudaMemcpy2D(ptr, size.first*sizeof(double),
		(char*)data.ptr + data.pitch*from.first + from.second*sizeof(double),
		data.pitch, size.first*sizeof(double), size.second, cudaMemcpyDeviceToHost
		);
	if (errMemset2D!=cudaSuccess) 
		cerr<<"CUDA ERROR "<<errMemset2D<<" in cudaMemset2D at Function getBorder\n";
}

void CudaGrid::getLeft(vector<double> &left) const {
	getBorder(&left[0], PointLong(0, 1), PointLong(1, left.size()));
}
void CudaGrid::getRight(vector<double> &right) const {
	getBorder(&right[0], PointLong(0, colsCount()-2), PointLong(1, right.size()));
}
void CudaGrid::getTop(vector<double> &top) const {
	getBorder(&top[0], PointLong(1, 0), PointLong(top.size(), 1));
}
void CudaGrid::getBottom(vector<double> &bottom) const  {
	getBorder(&bottom[0], PointLong(rowsCount()-2, 0), PointLong(bottom.size(), 1));
}
void initCacheForGridPoints(CudaGrid &grid, int left, int top) {
	cudaMalloc(&grid.cacheForGridPoints.first, grid.rowsCount()*sizeof(double));
	cudaMalloc(&grid.cacheForGridPoints.second, grid.colsCount()*sizeof(double));
	const long rowsCount = grid.rowsCount();
	const long colsCount = grid.colsCount();
	double *localX = new double[rowsCount];
	double *localY = new double[colsCount];
	for(int i = 0; i < grid.rowsCount(); ++i) {
		double I = static_cast<double>(i-top) + grid.rowsColsDelta.first;
		double xQ_COEF = I / (grid.totalRowsColsCount.first - 1);
		localX[i] = grid.rightTopCorner.first*grid.Q_COEF_tPOW(xQ_COEF);
	}
	for(int j = 0; j < grid.colsCount(); ++j){
		double J = static_cast<double>(j-left) + grid.rowsColsDelta.second;
		double yQ_COEF = J / (grid.totalRowsColsCount.second - 1);
		localY[j]= grid.rightTopCorner.second*grid.Q_COEF_tPOW(yQ_COEF);
	}
	int errCpy = cudaMemcpy(grid.cacheForGridPoints.first, localX, grid.rowsCount()*sizeof(double), cudaMemcpyHostToDevice);
	if (errCpy!=cudaSuccess) 
		cerr<<"CUDA ERROR "<<errCpy<<" in cudaMemcpy1 at Function initCacheForGridPoints\n";
	errCpy = cudaMemcpy(grid.cacheForGridPoints.second, localY, grid.colsCount()*sizeof(double), cudaMemcpyHostToDevice);
	if (errCpy!=cudaSuccess) 
		cerr<<"CUDA ERROR "<<errCpy<<" in cudaMemcpy2 at Function initCacheForGridPoints\n";
	delete[] localX;
	delete[] localY;
}

void CudaGrid::tranformToSimpleGrid(PointInt from, Grid &grid) const {
	double *resptr = grid.getData();
	getBorder(resptr, PointLong(from.first, from.second), PointLong(grid.getColumns(), grid.getRows()));
}

