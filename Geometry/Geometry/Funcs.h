#pragma once
#include "Curves.h"

/// ‘ункци€ поиска наименьшего рассто€ни€ между массивами точек, представл€ющих собой 2 различные кривые.
/// ¬озвращаемое значение - массив из двух точек(на 1й и 2й кривой соответственно)
template<typename CurveType1, typename CurveType2> DllExport vector<Point2D> FindClosestPoints(CurveType1 curve1, CurveType2 curve2, double eps = 1e-9)
{
	// Ќачальна€ инициализаци€ 
	vector<Point2D> to_return;
	to_return.resize(2);
	vector<Point2D> curClosestPoints;
	curClosestPoints.resize(2);
	int idxP1 = 0;
	int idxP2 = 0;
	double tFrom1 = 0, tFrom2 = 0;
	double tTo1 = 1, tTo2 = 1;
	auto p1 = curve1.GetCurvePoints();
	auto p2 = curve2.GetCurvePoints();
	to_return[0] = p1[0];
	to_return[1] = p2[0];
	// јлгоритм поиска минимального рассто€ни€ между кривыми путем перебора всех возможных пар точек
	do
	{
		Vector2D v = Vector2D(p1[0].e1, p1[0].e2, p2[0].e1, p2[0].e2);
		double min = v.GetLength();
		// ƒвойной цикл перебора всех пар точек
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
					curClosestPoints[0] = p1[i];
					curClosestPoints[1] = p2[j];
					min = dist;
				}
			}
		}
		// –асчет текущей погрешности
		double curEps = abs(Vector2D(curClosestPoints[0].e1, curClosestPoints[0].e2, curClosestPoints[1].e1, curClosestPoints[1].e2).GetLength() -
			Vector2D(to_return[0].e1, to_return[0].e2, to_return[1].e1, to_return[1].e2).GetLength());
		// ”словие выхода из цикла
		if (curEps < eps)
		{
			to_return[0] = curClosestPoints[0];
			to_return[1] = curClosestPoints[1];
			return to_return;
		}
		// ”величение точности нахождени€ минимального рассто€ни€
		else
		{
			to_return[0] = curClosestPoints[0];
			to_return[1] = curClosestPoints[1];
			if (idxP1 == 0) idxP1++;
			if (idxP1 == curve1.accuracy - 1) idxP1--;
			if (idxP2 == 0) idxP2++;
			if (idxP2 == curve2.accuracy - 1) idxP2--;
			// ”величение точности кривых
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

// ‘ункци€ поиска наименьшего рассто€ни€ между массивами точек, представл€ющих собой 2 различные кривые
template<typename CurveType1, typename CurveType2> DllExport vector<Point2D> FindCrossPoints(CurveType1 curve1, CurveType2 curve2, double eps = 1e-9)
{
	// Ќачальна€ инициализаци€
	vector<Point2D> to_return;
	auto p1 = curve1.GetCurvePoints();
	auto p2 = curve2.GetCurvePoints();

	// ƒвойной цикл поиск точек пересечени€
	for (int i = 0; i < curve1.accuracy - 1; i++)
	{
		Vector2D v1 = Vector2D(p1[i].e1, p1[i].e2, p1[i + 1].e1, p1[i + 1].e2);
		double a1 = 1.0 / v1.e1, b1 = -1.0 / v1.e2, c1 = p1[i].e1 / v1.e1 - p1[i].e2 / v1.e2;
		for (int j = 0; j < curve2.accuracy - 1; j++)
		{
			Vector2D v2 = Vector2D(p2[j].e1, p2[j].e2, p2[j + 1].e1, p2[j + 1].e2);
			double a2 = 1.0 / v2.e1, b2 = -1.0 / v2.e2, c2 = p2[j].e1 / v2.e1 - p2[j].e2 / v2.e2;
			// –ешение системы 2х2 на поиск точки пересечени€ двух векторов
			double det = a1 * b2 - a2 * b1;
			double root1 = c1 * b2 - c2 * b1;
			double root2 = a1 * c2 - a2 * c1;
			double e1 = root1 / det;
			double e2 = root2 / det;
			// ”словие искомой точки пересечени€ 
			if (((e1 > p1[i].e1 && e1 < p1[i + 1].e1) ||
				(e1 < p1[i].e1 && e1 > p1[i + 1].e1)) &&
				((e1 > p2[j].e1 && e1 < p2[j + 1].e1) ||
					(e1 < p2[j].e1 && e1 > p2[j + 1].e1)))
			{
				// Ќачальна€ инициализаци€ дл€ алгоритма увеличени€ точности точек пересечени€
				int idxP1 = i;
				int idxP2 = j;
				double curEps = 1;
				double tFrom1 = 0, tFrom2 = 0;
				double tTo1 = 1, tTo2 = 1;
				Point2D prevCrossPoint = Point2D(e1, e2);
				Point2D crossPoint;

				// ”точнение точности поиска точек пересечени€
				do
				{
					// ”величение точности кривых
					double tFromNew1 = tFrom1 + ((idxP1) / (double)(curve1.accuracy - 1)) * (tTo1 - tFrom1);
					double tToNew1 = tFrom1 + ((idxP1 + 1) / (double)(curve1.accuracy - 1)) * (tTo1 - tFrom1);
					double tFromNew2 = tFrom2 + ((idxP2) / (double)(curve2.accuracy - 1)) * (tTo2 - tFrom2);
					double tToNew2 = tFrom2 + ((idxP2 + 1) / (double)(curve2.accuracy - 1)) * (tTo2 - tFrom2);
					tFrom1 = tFromNew1;
					tTo1 = tToNew1;
					tFrom2 = tFromNew2;
					tTo2 = tToNew2;
					p1 = curve1.ImproveAccuracy(tFromNew1, tToNew1, curve1.accuracy);
					p2 = curve2.ImproveAccuracy(tFromNew2, tToNew2, curve2.accuracy);

					// ƒвойной цикл поиска уточненной точки пересечени€
					for (int ii = 0; ii < curve1.accuracy - 1; ii++)
					{
						v1 = Vector2D(p1[ii].e1, p1[ii].e2, p1[ii + 1].e1, p1[ii + 1].e2);
						a1 = 1.0 / v1.e1, b1 = -1.0 / v1.e2, c1 = p1[ii].e1 / v1.e1 - p1[ii].e2 / v1.e2;
						for (int jj = 0; jj < curve2.accuracy - 1; jj++)
						{
							v2 = Vector2D(p2[jj].e1, p2[jj].e2, p2[jj + 1].e1, p2[jj + 1].e2);
							a2 = 1.0 / v2.e1, b2 = -1.0 / v2.e2, c2 = p2[jj].e1 / v2.e1 - p2[jj].e2 / v2.e2;
							// –ешение системы 2х2 на поиск уточненной точки пересечени€ двух векторов
							det = a1 * b2 - a2 * b1;
							root1 = c1 * b2 - c2 * b1;
							root2 = a1 * c2 - a2 * c1;
							e1 = root1 / det;
							e2 = root2 / det;
							// ”словие искомой уточненной точки пересечени€ 
							if (((e1 > p1[ii].e1 && e1 < p1[ii + 1].e1) ||
								(e1 < p1[ii].e1 && e1 > p1[ii + 1].e1)) &&
								((e1 > p2[jj].e1 && e1 < p2[jj + 1].e1) ||
									(e1 < p2[jj].e1 && e1 > p2[jj + 1].e1)))
							{
								idxP1 = ii;
								idxP2 = jj;
								crossPoint = Point2D(e1, e2);
								curEps = Vector2D(prevCrossPoint.e1, prevCrossPoint.e2, crossPoint.e1, crossPoint.e2).GetLength();
								// ¬ыход из двойного цикла (лучше так, читаемость из-за goto в данном случае не страдает)
								goto cicleEnd;
							}
						}
					}
					// “очка выхода из двойного цикла поиска уточненной точки пересечени€
				cicleEnd:
					prevCrossPoint = crossPoint;
				} while (curEps > eps);
				p1 = curve1.GetCurvePoints();
				p2 = curve2.GetCurvePoints();
				to_return.push_back(crossPoint);
			}
		}
	}
	//  онец перебора всех точек
	return to_return;
}

// ‘ункци€ определени€ факториала числа n
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

// ћетод поиска градиента функции двух переменных
Vector2D gradient(double f, double fdt1, double fdt2, double dt1, double dt2)
{
	return Vector2D((fdt1 - f) / dt1, (fdt2 - f) / dt2);
}


