#pragma once
#include <iostream>
#include <vector>
#include "Objects2D.h"
#include "Funcs.h"

using namespace std;

// �����-�������� ��� ����� ������
class  DllExport Curve
{
public:
	// �������� �������������
	int accuracy;
	// ��������� ������ � ��������������� ���� (� ������� ���������� ����� ���� ���������)
	virtual Point2D MainFunc(double  t) = 0;
	// ��������� ������� ����� (Point2D *) ������
	virtual vector<Point2D> GetCurvePoints() = 0;
	// ��������� ���������� ������� ��������� �����
	virtual vector<vector<double>> GetCurveCoords() = 0;
	// ����� ���������� �������� ������������� ������ 
	virtual vector<Point2D> ImproveAccuracy(double t1, double t2, int accuracy) = 0;
};

// �����, �������������� ������ ����� 2�� �������
class  DllExport Bezier : public Curve
{
	// ������ �����
	vector<Point2D> points;
	// ������� �����
	vector<Point2D> basePoints;
public:
	// ��� ������� ����� � �������� ������������� (�� ��������� 100)
	Bezier(vector<Point2D> basePoints, int accuracy)
	{
		this->basePoints = basePoints;
		this->accuracy = accuracy;
		points.resize(accuracy);
		for (int i = 0; i < accuracy; i++)
		{
			points[i] = MainFunc(1.0 / (accuracy - 1) * i);
		}
	}
	~Bezier()
	{
		points.clear();
	}
	// ��������� ������ � ��������������� ���� 
	Point2D MainFunc(double  t) override
	{
		Point2D to_return = Point2D(0, 0);
		int n = basePoints.size() - 1;
		for (int i = 0; i < basePoints.size(); i++)
		{
			if (i != 0 && n - i != 0)
			{
				to_return = to_return + basePoints[i] * pow(t, i) * pow(1 - t, n - i) * (double)(fact(n)) / (fact(i) * fact(n - i));
			}
			else
			{
				if (i == 0)
					to_return = to_return + basePoints[i] * 1 * pow(1 - t, n - i) * (double)(fact(n)) / (fact(i) * fact(n - i));
				if (n - i == 0)
					to_return = to_return + basePoints[i] * pow(t, i) * 1 * (double)(fact(n)) / (fact(i) * fact(n - i));
			}

		}
		return to_return;
	}
	// ��������� ������� ����� (Point2D *) ������
	vector<Point2D> GetCurvePoints()override
	{
		return points;
	}
	// ��������� ���������� ������� ��������� �����
	vector<vector<double>> GetCurveCoords()override
	{
		vector<vector<double>> coords;
		coords.resize(2);
		coords[0].resize(accuracy);
		coords[1].resize(accuracy);
		for (int i = 0; i < accuracy; i++)
		{
			coords[0][i] = points[i].e1;
			coords[1][i] = points[i].e2;
		}
		return coords;
	}
	// ����� ���������� �������� ������������� ������ 
	vector<Point2D> ImproveAccuracy(double t1, double t2, int accuracy)override
	{
		vector<Point2D> to_return;
		to_return.resize(accuracy);
		for (int i = 0; i < accuracy; i++)
		{
			to_return[i] = MainFunc(t1 + (t2 - t1) / (accuracy - 1) * i);
		}
		return to_return;
	}
};

/// ������� ������ ����������� ���������� ����� ��������� �����, �������������� ����� 2 ��������� ������.
/// ������������ �������� - ������ �� ���� �����(�� 1� � 2� ������ ��������������)
DllExport vector<Point2D> FindClosestPoints(Curve* curve1, Curve* curve2, double eps = 1e-9)
{
	// ��������� ������������� 
	vector<Point2D> to_return;
	to_return.resize(2);
	vector<Point2D> curClosestPoints;
	curClosestPoints.resize(2);
	int idxP1 = 0;
	int idxP2 = 0;
	double tFrom1 = 0, tFrom2 = 0;
	double tTo1 = 1, tTo2 = 1;
	auto p1 = curve1->GetCurvePoints();
	auto p2 = curve2->GetCurvePoints();
	to_return[0] = p1[0];
	to_return[1] = p2[0];
	// �������� ������ ������������ ���������� ����� ������� ����� �������� ���� ��������� ��� �����
	do
	{
		Vector2D v = Vector2D(p1[0].e1, p1[0].e2, p2[0].e1, p2[0].e2);
		double min = v.GetLength();
		// ������� ���� �������� ���� ��� �����
		for (int i = 0; i < curve1->accuracy; i++)
		{
			for (int j = 0; j < curve2->accuracy; j++)
			{
				v = Vector2D(p1[i].e1, p1[i].e2, p2[j].e1, p2[j].e2);
				auto dist = v.GetLength();
				if (dist < min)
				{
					idxP1 = i;
					idxP2 = j;
					curClosestPoints[0] = p1[i];
					curClosestPoints[1] = p2[j];
					min = dist;
				}
			}
		}
		// ������ ������� �����������
		double curEps = abs(Vector2D(curClosestPoints[0].e1, curClosestPoints[0].e2, curClosestPoints[1].e1, curClosestPoints[1].e2).GetLength() -
			Vector2D(to_return[0].e1, to_return[0].e2, to_return[1].e1, to_return[1].e2).GetLength());
		// ������� ������ �� �����
		if (curEps < eps)
		{
			to_return[0] = curClosestPoints[0];
			to_return[1] = curClosestPoints[1];
			return to_return;
		}
		// ���������� �������� ���������� ������������ ����������
		else
		{
			to_return[0] = curClosestPoints[0];
			to_return[1] = curClosestPoints[1];
			if (idxP1 == 0) idxP1++;
			if (idxP1 == curve1->accuracy - 1) idxP1--;
			if (idxP2 == 0) idxP2++;
			if (idxP2 == curve2->accuracy - 1) idxP2--;
			// ���������� �������� ������
			double tFromNew1 = tFrom1 + ((idxP1 - 1) / (double)(curve1->accuracy - 1)) * (tTo1 - tFrom1);
			double tToNew1 = tFrom1 + ((idxP1 + 1) / (double)(curve1->accuracy - 1)) * (tTo1 - tFrom1);
			double tFromNew2 = tFrom2 + ((idxP2 - 1) / (double)(curve2->accuracy - 1)) * (tTo2 - tFrom2);
			double tToNew2 = tFrom2 + ((idxP2 + 1) / (double)(curve2->accuracy - 1)) * (tTo2 - tFrom2);
			tFrom1 = tFromNew1;
			tTo1 = tToNew1;
			tFrom2 = tFromNew2;
			tTo2 = tToNew2;
			p1 = curve1->ImproveAccuracy(tFromNew1, tToNew1, curve1->accuracy);
			p2 = curve2->ImproveAccuracy(tFromNew2, tToNew2, curve2->accuracy);
		}
	} while (true);
	return to_return;
}

// ������� ������ ����� ����������� ����� ��������� �����, �������������� ����� 2 ��������� ������
DllExport vector<Point2D> FindCrossPoints(Curve* curve1, Curve* curve2, double eps = 1e-9)
{
	std::cout << "\n���������� ������ ������ ����� ������� ������:\n";

	// ��������� �������������
	vector<Point2D> to_return;
	auto p1 = curve1->GetCurvePoints();
	auto p2 = curve2->GetCurvePoints();
	int accuracy1 = curve1->accuracy;
	int accuracy2 = curve2->accuracy;
	std::cout.precision(12);

	// ������� ���� ����� ����� �����������
	for (int i = 0; i < accuracy1 - 1; i++)
	{
		Vector2D v1 = Vector2D(p1[i].e1, p1[i].e2, p1[i + 1].e1, p1[i + 1].e2);
		double a1 = 1.0 / v1.e1, b1 = -1.0 / v1.e2, c1 = p1[i].e1 / v1.e1 - p1[i].e2 / v1.e2;
		for (int j = 0; j < accuracy2 - 1; j++)
		{
			Vector2D v2 = Vector2D(p2[j].e1, p2[j].e2, p2[j + 1].e1, p2[j + 1].e2);
			double a2 = 1.0 / v2.e1, b2 = -1.0 / v2.e2, c2 = p2[j].e1 / v2.e1 - p2[j].e2 / v2.e2;
			// ������� ������� 2�2 �� ����� ����� ����������� ���� ��������
			double det = a1 * b2 - a2 * b1;
			double root1 = c1 * b2 - c2 * b1;
			double root2 = a1 * c2 - a2 * c1;
			double e1 = root1 / det;
			double e2 = root2 / det;
			// ������� ������� ����� ����������� 
			if (((e1 > p1[i].e1 && e1 < p1[i + 1].e1) ||
				(e1 < p1[i].e1 && e1 > p1[i + 1].e1)) &&
				((e1 > p2[j].e1 && e1 < p2[j + 1].e1) ||
					(e1 < p2[j].e1 && e1 > p2[j + 1].e1)))
			{
				// ��������� ������������� ��� ��������� ���������� �������� ����� �����������
				int idxP1 = i;
				int idxP2 = j;
				double curEps = 1;
				double curEps1 = 1, curEps2 = 1;
				double tFrom1 = 0, tFrom2 = 0;
				double tTo1 = 1, tTo2 = 1;
				Point2D prevCrossPoint = Point2D(e1, e2);
				Point2D crossPoint;

				// ��������� �������� ������ ����� �����������
				do
				{
					// ���������� �������� ������
					double tFromNew1 = tFrom1 + ((idxP1) / (double)(accuracy1 - 1)) * (tTo1 - tFrom1);
					double tToNew1 = tFrom1 + ((idxP1 + 1) / (double)(accuracy1 - 1)) * (tTo1 - tFrom1);
					double tFromNew2 = tFrom2 + ((idxP2) / (double)(accuracy2 - 1)) * (tTo2 - tFrom2);
					double tToNew2 = tFrom2 + ((idxP2 + 1) / (double)(accuracy2 - 1)) * (tTo2 - tFrom2);
					tFrom1 = tFromNew1;
					tTo1 = tToNew1;
					tFrom2 = tFromNew2;
					tTo2 = tToNew2;
					p1 = curve1->ImproveAccuracy(tFromNew1, tToNew1, accuracy1);
					p2 = curve2->ImproveAccuracy(tFromNew2, tToNew2, accuracy2);

					// ������� ���� ������ ���������� ����� �����������
					for (int ii = 0; ii < accuracy1 - 1; ii++)
					{
						v1 = Vector2D(p1[ii].e1, p1[ii].e2, p1[ii + 1].e1, p1[ii + 1].e2);
						a1 = 1.0 / v1.e1, b1 = -1.0 / v1.e2, c1 = p1[ii].e1 / v1.e1 - p1[ii].e2 / v1.e2;
						for (int jj = 0; jj < accuracy2 - 1; jj++)
						{
							v2 = Vector2D(p2[jj].e1, p2[jj].e2, p2[jj + 1].e1, p2[jj + 1].e2);
							a2 = 1.0 / v2.e1, b2 = -1.0 / v2.e2, c2 = p2[jj].e1 / v2.e1 - p2[jj].e2 / v2.e2;
							// ������� ������� 2�2 �� ����� ���������� ����� ����������� ���� ��������
							det = a1 * b2 - a2 * b1;
							root1 = c1 * b2 - c2 * b1;
							root2 = a1 * c2 - a2 * c1;
							e1 = root1 / det;
							e2 = root2 / det;
							// ������� ������� ���������� ����� ����������� 
							if (((e1 > p1[ii].e1 && e1 < p1[ii + 1].e1) ||
								(e1 < p1[ii].e1 && e1 > p1[ii + 1].e1)) &&
								((e1 > p2[jj].e1 && e1 < p2[jj + 1].e1) ||
									(e1 < p2[jj].e1 && e1 > p2[jj + 1].e1)))
							{
								idxP1 = ii;
								idxP2 = jj;
								crossPoint = Point2D(e1, e2);
								curEps1 = Vector2D(p1[idxP1].e1, p1[idxP1].e2, crossPoint.e1, crossPoint.e2).GetLength();
								curEps2 = Vector2D(p2[idxP2].e1, p2[idxP2].e2, crossPoint.e1, crossPoint.e2).GetLength();
								curEps = curEps1;
								if (curEps1 < curEps2)
									curEps = curEps2;
								//curEps = Vector2D(prevCrossPoint.e1, prevCrossPoint.e2, crossPoint.e1, crossPoint.e2).GetLength();
								// ����� �� �������� ����� (����� ���, ���������� ��-�� goto � ������ ������ �� ��������)
								goto cicleEnd;
							}
						}
					}
					// ����� ������ �� �������� ����� ������ ���������� ����� �����������
				cicleEnd:
					prevCrossPoint = crossPoint;
				} while (curEps > eps);
				p1 = curve1->GetCurvePoints();
				p2 = curve2->GetCurvePoints();
				to_return.push_back(crossPoint);
				auto realPoint1 = curve1->MainFunc(tFrom1 + idxP1 * (tTo1 - tFrom1) / accuracy1);
				auto realPoint2 = curve2->MainFunc(tFrom2 + idxP2 * (tTo2 - tFrom2) / accuracy2);
				std::cout << std::endl;
				std::cout << "����� ����������� #" << to_return.size() << ", ����� ��������������� ��������� ������ 1:  " << realPoint1.e1 << " " << realPoint1.e2 << "\n";
				std::cout << "����� ����������� #" << to_return.size() << ", ����� ��������������� ��������� ������ 2:  " << realPoint2.e1 << " " << realPoint2.e2 << "\n";
				std::cout << "����� ����������� #" << to_return.size() << ", ��������� � ���������� ���������� �������: " << crossPoint.e1 << " " << crossPoint.e2 << "\n\n";
			}
		}
	}
	// ����� �������� ���� �����
	return to_return;
}

// ������� ������ ����� ����������� ����� ��������� �����, �������������� ����� 2 ��������� ������
DllExport vector<Point2D> FindCrossPointsViaEquations(Curve* curve1, Curve* curve2, double eps = 1e-9, int accuracy = 100)

{
	std::cout << "\n���������� ������ ������ ����� ��������� ������:\n";
	// ��������� �������������
	vector<Point2D> to_return;
	Point2D p1, p2, realPoint1, realPoint2;
	double dt1 = 1.0 / accuracy;
	double dt2 = 1.0 / accuracy;
	double curEps, newCurEps;
	double
		tFrom1 = 0,
		tTo1 = 1,
		tFrom2 = 0,
		tTo2 = 1;
	int idxFrom1, idxTo1, idxFrom2, idxTo2;
	// ������� ���� �������� ���� ��������� ������
	for (int i = 0; i < accuracy; i++)
	{
		p1 = curve1->MainFunc(i * dt1);
		for (int j = 0; j < accuracy; j++)
		{
			if (i == 0 && j == 0) continue;
			p2 = curve2->MainFunc(j * dt2);
			curEps = Vector2D(p1.e1, p1.e2, p2.e1, p2.e2).GetLength();
			// ������� ������������� ����� �� ���������� (����� 1�-1 ����� ������, �� ��������)
			// ������������ ��-�� ����� ����� ������������ ������ �����, �� ���� ��������� ��� ������ ��� ��������
			if (curEps < 1e-1)
			{
				// ���������� ������
				idxFrom1 = i - 1;
				idxFrom2 = j - 1;
				idxTo1 = idxFrom1 + 9;
				idxTo2 = idxFrom2 + 9;
				// ������� �������, ����� �� ���� ���������� ��������� � �������������� �������
				i = idxTo1;
				j = idxTo2;
				// ���� ��������� �������� �����
				do
				{
					// ����������������� �������� t1, t2 ������ �����
					double tFrom1new, tTo1new, tFrom2new, tTo2new;
					tFrom1new = tFrom1 + idxFrom1 * dt1;
					tTo1new = tFrom1 + idxTo1 * dt1;
					tFrom2new = tFrom2 + idxFrom2 * dt2;
					tTo2new = tFrom2 + idxTo2 * dt2;

					tFrom1 = tFrom1new;
					tTo1 = tTo1new;
					tFrom2 = tFrom2new;
					tTo2 = tTo2new;

					dt1 = (tTo1 - tFrom1) / accuracy;
					dt2 = (tTo2 - tFrom2) / accuracy;

					// �������� � ������, ������
					if (dt1 == 0 || dt2 == 0)
					{
						goto wrongRootCase;
					}


					// ������� ���� ������ ���������� � ��������� ����� ��������
					p1 = curve1->MainFunc(tFrom1);
					p2 = curve2->MainFunc(tFrom2);
					curEps = Vector2D(p1.e1, p1.e2, p2.e1, p2.e2).GetLength();
					for (int ii = 0; ii < accuracy; ii++)
					{
						p1 = curve1->MainFunc(tFrom1 + ii * dt1);
						for (int jj = 0; jj < accuracy; jj++)
						{
							if (ii == 0 && jj == 0) continue;
							p2 = curve2->MainFunc(tFrom2 + jj * dt2);
							newCurEps = Vector2D(p1.e1, p1.e2, p2.e1, p2.e2).GetLength();
							if (newCurEps < curEps)
							{
								curEps = newCurEps;
								idxFrom1 = ii;
								idxFrom2 = jj;
							}
						}
					}
					idxFrom1 -= 1;
					idxFrom2 -= 1;
					idxTo1 = idxFrom1 + 9;
					idxTo2 = idxFrom2 + 9;
				} while (curEps > eps);
				// ������ ��������� ������ � ������ ������
				to_return.push_back((p1 + p2) / 2);

				// ����������� ����������
				realPoint1 = curve1->MainFunc(tFrom1 + dt1 * (idxFrom1 + 1));
				realPoint2 = curve2->MainFunc(tFrom2 + dt2 * (idxFrom2 + 1));
				std::cout.precision(12);
				std::cout << "\n����� �����������, ��������������� ��������� ������ #1: " << realPoint1.e1 << " " << realPoint1.e2 << "\n";
				std::cout << "����� �����������, ��������������� ��������� ������ #2: " << realPoint2.e1 << " " << realPoint2.e2 << "\n";
				std::cout.precision(8);
				std::cout << "����������� � ������� ��������� f1(t1) - f2(t2) = 0:    " << (realPoint1 - realPoint2).e1 << " " << (realPoint1 - realPoint2).e2 << "\n\n";
				// ����� �������

			wrongRootCase:
				// ���������� �������� ���������� ��� "�������" ������ ������ �� ������
				dt1 = 1.0 / accuracy;
				dt2 = 1.0 / accuracy;
				p1 = curve1->MainFunc(i * dt1);
				p2 = curve2->MainFunc(j * dt2);
				tFrom1 = 0;
				tFrom2 = 0;
				tTo1 = 1;
				tTo2 = 1;
			}
		}
	}
	// ����� ������ ������
	return to_return;
}



