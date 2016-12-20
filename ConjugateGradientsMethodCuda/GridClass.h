#ifndef _GRIDCLASS_H_
#define _GRIDCLASS_H_
#include <cmath>
#include <map>
#include <cfloat>
#include <vector>
#include <ostream>
#include <iostream>
#include <iomanip>
#include "MatrixVectorClass.h"
#include "CommonTypes.h"
using namespace std;

int getPowerOfTwo(int val);
PointSizeT splitGridFunction(int rows, int cols, int sizePower);

static const PointDouble NAN_POINT(DBL_MAX, DBL_MAX);

class Grid {
private:
	PointLong rowsColsDelta;
	PointLong totalRowsColsCount;
	//C����
	Matrix data;
	//���� C����
	PointSizeT leftBottomCorner;
	PointSizeT rightTopCorner;
	//��� ����� ��������� c���� ��� ��������� ����������c��
	vector<vector<PointDouble> >  cacheForGridPoints;
	static const double Q_COEF; // q=3/2
	static const double Q_COEF_TWOPOW; // 2^q-1
	inline static double Q_COEF_tPOW(double t) {
		return (pow(1+t, Q_COEF) - 1) / Q_COEF_TWOPOW;
	}
	PointDouble countPoint(long i, long j) const;
	void initCacheForGridPoints() {
		cacheForGridPoints.clear();
		for(long i = 0; i < data.rowsCount()+2; ++i ) {
			cacheForGridPoints.push_back(vector<PointDouble>(data.colsCount()+2));
		}
		for(long i = 0; i < data.rowsCount()+2; ++i){
			for(long j = 0; j < data.colsCount()+2; ++j){
				double I = static_cast<double>(i-1) + rowsColsDelta.first;
				double J = static_cast<double>(j-1) + rowsColsDelta.second;
				double xQ_COEF = I / (getTotalRowsCount() - 1);
				double yQ_COEF = J / (getTotalColsCount() - 1);
				double x = rightTopCorner.first*Q_COEF_tPOW(xQ_COEF);
				double y = rightTopCorner.second*Q_COEF_tPOW(yQ_COEF);
				cacheForGridPoints[i][j] = PointDouble(x,y);
			}
		}
	}
public:
	Grid():leftBottomCorner(0,0), rightTopCorner(0,0), rowsColsDelta(0,0), totalRowsColsCount(0,0){}

	Grid(PointLong subGridRowsCols, 
		PointSizeT leftBottomCorner_, 
		PointSizeT rightTopCorner_, 
		PointLong rowsColsDelta_, 
		PointLong totalRowsColsCount_);

	Grid(PointLong subGridRowsCols, 
		double *data, 
		long dataCount, 
		PointSizeT leftBottomCorner_, 
		PointSizeT rightTopCorner_, 
		PointLong rowsColsDelta_, 
		PointLong totalRowsColsCount_);

	Grid(const Matrix &data, 
		PointSizeT leftBottomCorner_, 
		PointSizeT rightTopCorner_, 
		PointLong rowsColsDelta_, 
		PointLong totalRowsColsCount_);

	Grid(const Grid &other):
	data(other.data),
		leftBottomCorner(other.leftBottomCorner),
		rightTopCorner(other.rightTopCorner),
		rowsColsDelta(other.rowsColsDelta),
		totalRowsColsCount(other.totalRowsColsCount),
		cacheForGridPoints(other.cacheForGridPoints)
	{
	}

	Grid &operator=(const Grid &other) {
		data = other.data;
		leftBottomCorner = other.leftBottomCorner;
		rightTopCorner = other.rightTopCorner;
		rowsColsDelta.first = other.rowsColsDelta.first;
		rowsColsDelta.second = other.rowsColsDelta.second;
		totalRowsColsCount.first = other.totalRowsColsCount.first;
		totalRowsColsCount.second = other.totalRowsColsCount.second;
		cacheForGridPoints = other.cacheForGridPoints;
		return *this;
	}

	// �������� ���� c����
	PointSizeT getLeftBottomCorner() const {return leftBottomCorner; }
	PointSizeT getRightTopCorner() const {return rightTopCorner; }
	// �������� c������� ��� ������� ���c���
	long getRowsDelta() const { return rowsColsDelta.first; }
	long getColumnsDelta() const {return rowsColsDelta.second; }
	// �������� ������ c���� � ������� ����������c� c�������
	long getTotalRowsCount() const {return totalRowsColsCount.first; }
	long getTotalColsCount() const {return totalRowsColsCount.second; }
	// �������� ������ ���c����
	long getRows() const { return data.rowsCount(); }
	long getColumns() const { return data.colsCount(); }
	// �������� ��������� c�����/c������ ���c����
	vector<double> getTopRowBorder() const {return data.getRow(0);}
	vector<double> getBottomRowBorder() const {return data.getRow(data.rowsCount() - 1);}
	vector<double> getLeftColBorder() const {return data.getCol(0);}
	vector<double> getRightColBorder() const {return data.getCol(data.colsCount() - 1);}
	// �������� �� �������� c����
	friend void initGridBorder(Grid &grid, Function FunctionPhi);
	// �������� ���������� ����� ��������� c���� �� ����
	PointDouble getPoint(long i, long j) const{
		return cacheForGridPoints[i+1][j+1];
	}
	// �������� �������� �����c��������� ������� � ����� c����
	double operator()(long i, long j) const {
		return data(i,j);
	}
	double &operator()(long i, long j) { return data(i,j); }
	//�������� ��������� �� ��cc�� �������� �����c. ������� � ���c����
	double *getData() {return data.basePtr();}
	const double *getData() const {return data.basePtr();}

	// ��������� c���� �� �����cc����
	friend map<int, Grid> splitGrid(const Grid &grid, long procCount);
	// C��� c���� c �����cc����
	friend Grid collectGrid(const map<int, Grid> &subGrids);

	friend ostream &operator<< (ostream &os, const Grid& m) {
		for(int i = 0; i<m.getRows(); ++i){
			for(int j = 0; j < m.getColumns(); ++j){
				os <<setw(10) << m(i,j) << "\t";
			}
			os << "\n";
		}
		return os;
	}
	friend ostream &print(ostream& os, const Grid &grid) {
		for (long i = 0; i< grid.getRows(); ++i){
			for(long j = 0; j < grid.getColumns(); ++j) {
				PointDouble p = grid.getPoint(i, j);
				double val = grid(i,j);
				os << p.first << "\t" << p.second << "\t" << val <<"\n";
			}
		}
		return os;
	}
};
#endif

