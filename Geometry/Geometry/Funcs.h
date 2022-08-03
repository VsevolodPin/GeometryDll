#pragma once
#include "Curves.h"

// ������� ������ ����������� ���������� ����� ��������� �����, �������������� ����� 2 ��������� ������
template<typename CurveType1, typename CurveType2> DllExport vector<Point2D> FindClosestPoints(CurveType1 curve1, CurveType2 curve2, double eps = 1e-9)
{
	vector<Point2D> to_return;
	to_return.resize(2);
	vector<Point2D> closestPoints;
	closestPoints.resize(2);
	int idxP1 = 0;
	int idxP2 = 0;
	double tFrom1 = 0, tFrom2 = 0;
	double tTo1 = 1, tTo2 = 1;
	auto p1 = curve1.GetCurvePoints();
	auto p2 = curve2.GetCurvePoints();
	to_return[0] = p1[0];
	to_return[1] = p2[0];
	do
	{
		Vector2D v = Vector2D(p1[0].e1, p1[0].e2, p2[0].e1, p2[0].e2);
		double  min = v.GetLength();

		for (int i = 0; i < curve1.accuracy; i++)
		{
			for (int j = 0; j < curve2.accuracy; j++)
			{
				v = Vector2D(p1[i].e1, p1[i].e2, p2[j].e1, p2[j].e2);
				auto dist = v.GetLength();
				if (dist < min)
				{
					idxP1 = i;
					idxP2 = j;
					closestPoints[0] = p1[i];
					closestPoints[1] = p2[j];
					min = dist;
				}
			}
		}

		double curEps = abs(Vector2D(closestPoints[0].e1, closestPoints[0].e2, closestPoints[1].e1, closestPoints[1].e2).GetLength() -
			Vector2D(to_return[0].e1, to_return[0].e2, to_return[1].e1, to_return[1].e2).GetLength());
		if (curEps < eps)
		{
			to_return[0] = closestPoints[0];
			to_return[1] = closestPoints[1];
			return to_return;
		}
		else
		{
			to_return[0] = closestPoints[0];
			to_return[1] = closestPoints[1];
			if (idxP1 == 0) idxP1++;
			if (idxP1 == curve1.accuracy - 1) idxP1--;
			if (idxP2 == 0) idxP2++;
			if (idxP2 == curve2.accuracy - 1) idxP2--;
			double tFromNew1 = tFrom1 + ((idxP1 - 1) / (double)(curve1.accuracy - 1)) * (tTo1 - tFrom1);
			double tToNew1 = tFrom1 + ((idxP1 + 1) / (double)(curve1.accuracy - 1)) * (tTo1 - tFrom1);
			double tFromNew2 = tFrom2 + ((idxP2 - 1) / (double)(curve2.accuracy - 1)) * (tTo2 - tFrom2);
			double tToNew2 = tFrom2 + ((idxP2 + 1) / (double)(curve2.accuracy - 1)) * (tTo2 - tFrom2);
			tFrom1 = tFromNew1;
			tTo1 = tToNew1;
			tFrom2 = tFromNew2;
			tTo2 = tToNew2;
			p1 = curve1.ImproveAccuracy(tFromNew1, tToNew1, curve1.accuracy);
			p2 = curve2.ImproveAccuracy(tFromNew2, tToNew2, curve2.accuracy);
		}

	} while (true);
	return to_return;
}

// ������� ����������� ���������� ����� n
int fact(int n)
{
	if (n <= 1) return 1;
	int to_return = 1;
	for (int i = n; i > 0; i--)
	{
		to_return *= i;
	}
	return to_return;
}

