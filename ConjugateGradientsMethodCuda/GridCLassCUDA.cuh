#ifndef _CUDA_GRIDCLASS_CUH_
#define _CUDA_GRIDCLASS_CUH_
#include "GridClass.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>

struct PointDoublePtr{
	double* first;
	double* second;
};

class CudaGrid {
private:
	PointLong rowsColsDelta;
	PointLong totalRowsColsCount;
	//C����
	cudaPitchedPtr data;
	//���� C����
	PointSizeT leftBottomCorner;
	PointSizeT rightTopCorner;
	__host__ __device__ inline static double Q_COEF_tPOW(double t) {
		const double Q_COEF = 1.5;
		const double Q_COEF_TWOPOW = 1.82842712474;
		return (pow(1+t, Q_COEF) - 1) / Q_COEF_TWOPOW;
	}
	// �c�������� ��������� �������� ���c���� �����cc���
	void setBorder(PointLong from, PointLong size, const double *data);
	// �������� ��������� �������� ���c���� �����cc���
	void getBorder(double *data, PointLong from, PointLong size) const;
	//��� ����� ��������� c���� ��� ��������� ����������c��
	PointDoublePtr cacheForGridPoints;
public:
	friend void initCacheForGridPoints(CudaGrid &grid, int left, int top);
	CudaGrid();
	void setGridOptions(PointLong rowsColsCount_, PointSizeT leftBottomCorner_,
		PointSizeT rightTopCorner_, PointLong rowsColsDelta_, PointLong totalRowsColsCount_);
	CudaGrid(PointLong rowsColsCount_, PointSizeT leftBottomCorner_,
		PointSizeT rightTopCorner_, PointLong rowsColsDelta, PointLong totalRowsColsCount);
	// �������� �������� �����c. ������� � ����� c����
	__device__ double operator()(long i, long j) const {
		return *((double*)((char*)data.ptr + i*data.pitch) + j);
	}
	__device__ double &operator()(long i, long j) {
		return *((double*)((char*)data.ptr + i*data.pitch) + j);
	}
	// �������� ����� ��������� c���� �� Ox
	__device__ double getGridSplitOx(long i) const {
		if(!cacheForGridPoints.first) {
			double I = static_cast<double>(i) + rowsColsDelta.first;
			double xQ_COEF = I / (totalRowsColsCount.first - 1);
			return rightTopCorner.first*Q_COEF_tPOW(xQ_COEF);
		} else {
			return cacheForGridPoints.first[i];
		}
	}
	//�������� ����� ��������� c���� �� Oy
	__device__ double getGridSplitOy(long j) const {
		if(!cacheForGridPoints.second) {
			double J = static_cast<double>(j) + rowsColsDelta.second;
			double yQ_COEF = J / (totalRowsColsCount.second - 1);
			return rightTopCorner.second*Q_COEF_tPOW(yQ_COEF);
		} else {
			return cacheForGridPoints.second[j];
		}
	}
	//�������� c������ ��� ��������� c���� � ����� �� Ox
	__device__ double getAvgStepOx(long i) const {
		double prevX = getGridSplitOx(i-1);
		double nextX = getGridSplitOx(i+1);
		return ((nextX - prevX)) / 2;
	}
	//�������� c������ ��� ��������� c���� � ����� �� Oy
	__device__ double getAvgStepOy(long j) const {
		double prevY = getGridSplitOy(j-1);
		double nextY = getGridSplitOy(j+1);
		return (nextY - prevY) / 2; 
	}
	// �c�������� ��������� �������� ���c���� �����cc��� c ��������� c������
	void setLeftBorder(const vector<double> &left, bool delta=false);
	void setRightBorder(const vector<double> &right, bool delta=false);
	void setTopBorder(const vector<double> &top, bool delta=false);
	void setBottomBorder(const vector<double> &bottom, bool delta=false);
	// �������� ��������� �������� ���c���� �����cc��� c ��������� c������
	void getLeft(vector<double> &left) const;
	void getRight(vector<double> &right) const;
	void getTop(vector<double> &top) const;
	void getBottom(vector<double> &bottom) const;
	// �������� ������c��� c���� c������� c����
	__host__ __device__ const long rowsCount() const { return data.xsize; }
	__host__ __device__ const long colsCount() const { return data.ysize; }
	// �������� c������� ��� ������� ���c���
	__host__ __device__ long getRowsDelta() const { return rowsColsDelta.first; }
	__host__ __device__ long getcolsDelta() const { return rowsColsDelta.second; }
	// �������� ������ c���� � ������� ����������c� c�������
	__host__ __device__ long getTotalRowsCount() const { return totalRowsColsCount.first; }
	__host__ __device__ long getTotalColsCount() const { return totalRowsColsCount.second; }

	//������������� CUDA c���� � �������
	void tranformToSimpleGrid(PointInt from, Grid &grid) const;

};
#endif
