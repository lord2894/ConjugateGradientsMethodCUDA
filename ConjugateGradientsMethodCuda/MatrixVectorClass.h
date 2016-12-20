#ifndef _MATRIXVECTORCLASS_H_
#define _MATRIXVECTORCLASS_H_
#define __volatile volatile
#include <map>
#include <cmath>
#include <vector>
using namespace std;

class Matrix {
private:
    vector<double> data;
    long rows;
    long cols;
public:
    Matrix(): rows(0), cols(0){}
    Matrix(long r, long c): rows(r), cols(c){
		data.resize(rows*cols);
    }

    Matrix(long r, long c, double val):rows(r), cols(c) {
		data.reserve(rows*cols);
        for(long i = 0; i < rows*cols; ++i) {
            data.push_back(val);
        }
    }

    Matrix(double *data_, long dataCount, long r, long c): rows(r), cols(c) {
		data.reserve(rows*cols);
		for (int i = 0; i < dataCount; ++i)
			data.push_back(data_[i]);
	}

    Matrix(const Matrix& other): rows(other.rows), cols(other.cols) {
		data.resize(rows*cols);
        for(long i = 0; i < rows; ++i) {
            for(long j = 0; j < cols; ++j) {
                this->operator()(i,j) = other(i,j);
            }
        }
    }
    Matrix &operator=(const Matrix &other) {
        rows = other.rows;
        cols = other.cols;
		data.resize(rows*cols);
        for(long i = 0; i < rows; ++i){
            for(long j = 0; j < cols; ++j){
                this->operator()(i,j) = other(i,j);
            }
        }
        return *this;
    }
    double &operator() (long i, long j) {
		return data[i*cols + j];
    }
    double operator() (long i, long j) const {
		return data[i*cols + j];
    }
    double *basePtr() {
        return &data[0];
    }
    const double *basePtr() const {
		return &data[0];
    }
    long rowsCount() const { return rows; }
    long colsCount() const { return cols; }
    vector<double> getRow(long index) const {
		vector<double> result;
		result.reserve(cols);
        for(long i = 0; i < cols; ++i) {
            result.push_back(data[index*cols + i]);
        }
        return result;
    }
	vector<double> getCol(long index) const {
		vector<double> result;
		result.reserve(rows);
        for(long i = 0; i < rows; ++i){
			result.push_back(data[i*cols + index]);
        }
        return result;
    }
    Matrix getSubmat(long startr, long startc, long rws, long cls) const {
        if(startr + rws > rows){
            rws = rows - startr;
        }
        if(startc + cls > cols) {
            cls = cols - startc;
        }
        Matrix result(rws,cls);
        for(long i = 0; i < rws; ++i){
            for(long j = 0; j < cls; ++j) {
                result(i,j) = this->operator()(i + startr, j+startc);
            }
        }
        return result;
    }
    friend map<pair<int,int>,Matrix> split(const Matrix& m, long rows, long cols) {
       long rElem = (long)(m.rowsCount() / (double)rows + 1);
       long cElem = (long)(m.colsCount() / (double)cols + 1);
       map<pair<int,int>, Matrix> result;
       for(long i = 0; i < rows; ++i){
           for(long j = 0; j < cols; ++j){
               result[pair<int,int>(i,j)] = m.getSubmat(i*rElem, j*cElem, rElem, cElem);
           }
       }
       return result;
    }
};
#endif
