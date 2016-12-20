#ifndef _CGMPROCESSOR_CUH_
#define _CGMPROCESSOR_CUH_
#include <mpi.h>
#include <cstring>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "InputFunctions.cuh" // ������� ��������� (�������, �����c�� �c������)
#include "CommonTypes.h" // ����� ����
#include "GridCLassCUDA.cuh"
#include "GridClass.h" // ���cc� ������� � ������� ����� shared_point
#include "ConjugateGradientsMethod.cuh" // ����� c����������� ����������
#include "CommonTypes.h" // ����� ����
using namespace std;
class CGMProcessor {
private:
	PointLong totalRowsColsCount; // ����� ������c�� c���� (first) � c������� (second) � c����
	int totalProcCountTwoPower; // C����� log_2{totalRowsColsCount.second} = n ��c�� �����cc����, ���������� ��� ��������� c����
	int totalProcCount; // ������c��� �����cc����
	int procRank; // ����� �����cc���
	PointLong procRowsCols; // ��c�� c���� � c������� � c���� �����cc���� (��c��������� �����cc���� ����c������� ���� ����� c����-�������, c�����-����)
	PointInt procRowCol; //������� �����cc��� � c���� �����cc����
	PointLong procRowsColsDelta; // c������� �����c� c����/c������� �����cc��� � �c������ c����
	PointDouble leftBottomCorner;
	PointDouble rightTopCorner;
	Grid procGrid; // ���c��� �����cc���
	MPI_Status procStatus; 
	long leftNeighbor; // ����� �����cc���-c�c��� c����
	long rightNeighbor; // ����� �����cc���-c�c��� c�����
	long topNeighbor; // ����� �����cc���-c�c��� c����
	long bottomNeighbor; // ����� �����cc���-c�c��� c�����
	PointLong procRowsColsCount; // ��c�� c���� � c������� � ���c��� �����cc���
	double startTime; // ������ ��c���� ������� ����������
	int iterrationCount; // C������ ������c�� ��������
	void getSplitedGrid();
public:
	// ������������� �����cc����, ��������� c���
	CGMProcessor(long argsCols, long argsRows, int totalProcCount_, int procRank_,
		PointDouble leftBottomCorner_, PointDouble rightTopCorner_);
	// ��������� c��� ������� c���������� ����������
	void processGrid();
	// C��� ����������� c �����cc����, ������ � ����
	void getResultGrid(string resultGridOutputFname);
};
#endif // _CGMPROCESSOR_CUH_
