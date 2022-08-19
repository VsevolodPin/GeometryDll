#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Objects2D.h"
#include "Funcs.h"

using namespace std;

// Класс-родитель для любой кривой
class  DllExport Curve
{
public:
	// Уравнение кривой в параметрическом виде (у каждого наследника будет свое уравнение)
	virtual const Point2D F(double  t) = 0;
	// Уравнение кривой в параметрическом виде (у каждого наследника будет свое уравнение)
	virtual const Point2D dF(double  t) = 0;
	// Получение массива точек (Point2D *) кривой
	virtual const vector<Point2D> GetCurvePoints(int N) = 0;
	// Получение двумерного массива координат точек
	virtual const vector<vector<double>> GetCurveCoords(int N) = 0;
	// Метод увеличения точности моделирования кривой 
	virtual const vector<Point2D> ImproveAccuracy(double t1, double t2, int accuracy) = 0;
	// Виртуальный деструктор
	virtual ~Curve() = default;
};

// Класс, представляющий кривую Безье 2го порядка
class  DllExport Bezier : public Curve
{
	// Опорные точки
	vector<Point2D> basePoints;
public:
	// Три опорные точки и точность моделирования (по умолчанию 100)
	Bezier(vector<Point2D> basePoints)
	{
		this->basePoints = basePoints;
	}
	/// Конструктор через файловый поток с опорными точками
	/// Структура файла должна иметь следующий вид:
	/// p0.x p0.y
	/// ...
	/// pi.x pi.y
	Bezier(ifstream* stream, bool closeStream = false)
	{
		//int i = 0;
		if (stream->is_open())
			while (!stream->eof())
			{
				double x, y;
				*stream >> x;
				*stream >> y;
				basePoints.push_back(Point2D(x, y));
			}
		if (closeStream)
			stream->close();
	}
	// Деструктор
	~Bezier() = default;
	// Оператор освобождения памяти
	void operator delete (void* ptr)
	{
		free(ptr);
	}
	// Уравнение кривой в параметрическом виде 
	const Point2D F(double  t) override
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
	// Значение производной в точке t
	const Point2D dF(double  t) override
	{
		Point2D to_return = Point2D(0, 0);
		int n = basePoints.size() - 1;
		for (int i = 0; i < basePoints.size(); i++)
		{
			if (i != 0 && n - i != 0)
			{
				to_return = to_return + basePoints[i] * pow(t, i) * pow(1 - t, n - i) * (double)(fact(n)) / (fact(i) * fact(n - i)) * ((i - n * t) / (t - t * t));
			}
			else
			{
				if (i == 0)
					to_return = to_return + basePoints[i] * 1 * pow(1 - t, n - i) * (double)(fact(n)) / (fact(i) * fact(n - i)) * ((i - n * t) / (t - t * t));
				if (n - i == 0)
					to_return = to_return + basePoints[i] * pow(t, i) * 1 * (double)(fact(n)) / (fact(i) * fact(n - i)) * ((i - n * t) / (t - t * t));
			}
		}
		return to_return;
	}

	// Получение массива точек (Point2D *) кривой
	const vector<Point2D> GetCurvePoints(int N)override
	{
		vector<Point2D> points;
		points.resize(N);
		for (int i = 0; i < N; i++)
			points[i] = F(1.0 / (N - 1) * i);
		return points;
	}
	// Получение двумерного массива координат точек
	const vector<vector<double>> GetCurveCoords(int N)override
	{
		vector<vector<double>> coords;
		coords.resize(2);
		coords[0].resize(N);
		coords[1].resize(N);
		for (int i = 0; i < N; i++)
		{
			auto p = this->F(1.0 / (N - 1) * i);
			coords[0][i] = p.e1;
			coords[1][i] = p.e2;
		}
		return coords;
	}
	// Метод увеличения точности моделирования кривой 
	const vector<Point2D> ImproveAccuracy(double t1, double t2, int accuracy)override
	{
		vector<Point2D> to_return;
		to_return.resize(accuracy);
		for (int i = 0; i < accuracy; i++)
		{
			to_return[i] = F(t1 + (t2 - t1) / (accuracy - 1) * i);
		}
		return to_return;
	}
};

/// Функция поиска наименьшего расстояния между массивами точек, представляющих собой 2 различные кривые.
/// Возвращаемое значение - массив из двух точек(на 1й и 2й кривой соответственно)
DllExport vector<Point2D> FindClosestPoints(Curve* curve1, Curve* curve2, double eps = 1e-9)
{
	// Начальная инициализация 
	vector<Point2D> to_return;
	to_return.resize(2);
	vector<Point2D> curClosestPoints;
	curClosestPoints.resize(2);
	int idxP1 = 0;
	int idxP2 = 0;
	double tFrom1 = 0, tFrom2 = 0;
	double tTo1 = 1, tTo2 = 1;
	int accuracy = 100;
	auto p1 = curve1->GetCurvePoints(100);
	auto p2 = curve2->GetCurvePoints(100);
	to_return[0] = p1[0];
	to_return[1] = p2[0];
	// Алгоритм поиска минимального расстояния между кривыми путем перебора всех возможных пар точек
	do
	{
		Vector2D v = Vector2D(p1[0].e1, p1[0].e2, p2[0].e1, p2[0].e2);
		double min = v.GetLength();
		// Двойной цикл перебора всех пар точек
		for (int i = 0; i < accuracy; i++)
		{
			for (int j = 0; j < accuracy; j++)
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
		// Расчет текущей погрешности
		double curEps = abs(Vector2D(curClosestPoints[0].e1, curClosestPoints[0].e2, curClosestPoints[1].e1, curClosestPoints[1].e2).GetLength() -
			Vector2D(to_return[0].e1, to_return[0].e2, to_return[1].e1, to_return[1].e2).GetLength());
		// Условие выхода из цикла
		if (curEps < eps)
		{
			to_return[0] = curClosestPoints[0];
			to_return[1] = curClosestPoints[1];
			return to_return;
		}
		// Увеличение точности нахождения минимального расстояния
		else
		{
			to_return[0] = curClosestPoints[0];
			to_return[1] = curClosestPoints[1];
			if (idxP1 == 0) idxP1++;
			if (idxP1 == accuracy - 1) idxP1--;
			if (idxP2 == 0) idxP2++;
			if (idxP2 == accuracy - 1) idxP2--;
			// Увеличение точности кривых
			double tFromNew1 = tFrom1 + ((idxP1 - 1) / (double)(accuracy - 1)) * (tTo1 - tFrom1);
			double tToNew1 = tFrom1 + ((idxP1 + 1) / (double)(accuracy - 1)) * (tTo1 - tFrom1);
			double tFromNew2 = tFrom2 + ((idxP2 - 1) / (double)(accuracy - 1)) * (tTo2 - tFrom2);
			double tToNew2 = tFrom2 + ((idxP2 + 1) / (double)(accuracy - 1)) * (tTo2 - tFrom2);
			tFrom1 = tFromNew1;
			tTo1 = tToNew1;
			tFrom2 = tFromNew2;
			tTo2 = tToNew2;
			p1 = curve1->ImproveAccuracy(tFromNew1, tToNew1, accuracy);
			p2 = curve2->ImproveAccuracy(tFromNew2, tToNew2, accuracy);
		}
	} while (true);
	return to_return;
}

// Функция поиска точки пересечения между массивами точек, представляющих собой 2 различные кривые
DllExport vector<Point2D> FindCrossPointsViaSegments(Curve* curve1, Curve* curve2, double eps = 1e-9)
{
	std::cout << "\nРезультаты поиска корней через отрезки кривых:\n";

	// Начальная инициализация
	vector<Point2D> to_return;
	auto p1 = curve1->GetCurvePoints(100);
	auto p2 = curve2->GetCurvePoints(100);
	int accuracy1 = 100;
	int accuracy2 = 100;
	std::cout.precision(12);

	// Двойной цикл поиск точек пересечения
	for (int i = 0; i < accuracy1 - 1; i++)
	{
		Vector2D v1 = Vector2D(p1[i].e1, p1[i].e2, p1[i + 1].e1, p1[i + 1].e2);
		double a1 = 1.0 / v1.e1, b1 = -1.0 / v1.e2, c1 = p1[i].e1 / v1.e1 - p1[i].e2 / v1.e2;
		for (int j = 0; j < accuracy2 - 1; j++)
		{
			Vector2D v2 = Vector2D(p2[j].e1, p2[j].e2, p2[j + 1].e1, p2[j + 1].e2);
			double a2 = 1.0 / v2.e1, b2 = -1.0 / v2.e2, c2 = p2[j].e1 / v2.e1 - p2[j].e2 / v2.e2;
			// Решение системы 2х2 на поиск точки пересечения двух векторов
			double det = a1 * b2 - a2 * b1;
			double root1 = c1 * b2 - c2 * b1;
			double root2 = a1 * c2 - a2 * c1;
			double e1 = root1 / det;
			double e2 = root2 / det;
			// Условие искомой точки пересечения 
			if (((e1 > p1[i].e1 && e1 < p1[i + 1].e1) ||
				(e1 < p1[i].e1 && e1 > p1[i + 1].e1)) &&
				((e1 > p2[j].e1 && e1 < p2[j + 1].e1) ||
					(e1 < p2[j].e1 && e1 > p2[j + 1].e1)))
			{
				// Начальная инициализация для алгоритма увеличения точности точек пересечения
				int idxP1 = i;
				int idxP2 = j;
				double curEps = 1;
				double curEps1 = 1, curEps2 = 1;
				double tFrom1 = 0, tFrom2 = 0;
				double tTo1 = 1, tTo2 = 1;
				Point2D prevCrossPoint = Point2D(e1, e2);
				Point2D crossPoint;

				// Уточнение точности поиска точек пересечения
				do
				{
					// Увеличение точности кривых
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

					// Двойной цикл поиска уточненной точки пересечения
					for (int ii = 0; ii < accuracy1 - 1; ii++)
					{
						v1 = Vector2D(p1[ii].e1, p1[ii].e2, p1[ii + 1].e1, p1[ii + 1].e2);
						a1 = 1.0 / v1.e1, b1 = -1.0 / v1.e2, c1 = p1[ii].e1 / v1.e1 - p1[ii].e2 / v1.e2;
						for (int jj = 0; jj < accuracy2 - 1; jj++)
						{
							v2 = Vector2D(p2[jj].e1, p2[jj].e2, p2[jj + 1].e1, p2[jj + 1].e2);
							a2 = 1.0 / v2.e1, b2 = -1.0 / v2.e2, c2 = p2[jj].e1 / v2.e1 - p2[jj].e2 / v2.e2;
							// Решение системы 2х2 на поиск уточненной точки пересечения двух векторов
							det = a1 * b2 - a2 * b1;
							root1 = c1 * b2 - c2 * b1;
							root2 = a1 * c2 - a2 * c1;
							e1 = root1 / det;
							e2 = root2 / det;
							// Условие искомой уточненной точки пересечения 
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
								// Выход из двойного цикла (лучше так, читаемость из-за goto в данном случае не страдает)
								goto cicleEnd;
							}
						}
					}
					// Точка выхода из двойного цикла поиска уточненной точки пересечения
				cicleEnd:
					prevCrossPoint = crossPoint;
				} while (curEps > eps);
				p1 = curve1->GetCurvePoints(100);
				p2 = curve2->GetCurvePoints(100);
				to_return.push_back(crossPoint);
				auto realPoint1 = curve1->F(tFrom1 + idxP1 * (tTo1 - tFrom1) / accuracy1);
				auto realPoint2 = curve2->F(tFrom2 + idxP2 * (tTo2 - tFrom2) / accuracy2);
				std::cout << std::endl;
				std::cout << "Точка пересечения #" << to_return.size() << ", точно соответствующая уравнению кривой 1:  " << realPoint1.e1 << " " << realPoint1.e2 << "\n";
				std::cout << "Точка пересечения #" << to_return.size() << ", точно соответствующая уравнению кривой 2:  " << realPoint2.e1 << " " << realPoint2.e2 << "\n";
				std::cout << "Точка пересечения #" << to_return.size() << ", найденная в результате численного решения: " << crossPoint.e1 << " " << crossPoint.e2 << "\n\n";
			}
		}
	}
	// Конец перебора всех точек
	return to_return;
}

// Функция поиска точки пересечения между массивами точек, представляющих собой 2 различные кривые
DllExport vector<Point2D> FindCrossPointsViaEquations(Curve* curve1, Curve* curve2, double eps = 1e-9, int accuracy = 100)
{
	std::cout << "\nРезультаты поиска корней через уравнения кривых:\n";
	// Начальная инициализация
	int count = 0;
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
	// Двойной цикл перебора всех возможных корней
	for (int i = 0; i < accuracy; i++)
	{
		p1 = curve1->F(i * dt1);
		for (int j = 0; j < accuracy; j++)
		{
			if (i == 0 && j == 0) continue;
			p2 = curve2->F(j * dt2);
			curEps = Vector2D(p1.e1, p1.e2, p2.e1, p2.e2).GetLength();
			// Условие существования корня на промежутке (число 1е-1 взято наобум, но работает)
			// Теоретически из-за этого могут существовать ложные корни, но пока непонятно как решить эту проблему
			if (curEps < 1e-1)
			{
				// Локализуем корень
				idxFrom1 = i - 1;
				idxFrom2 = j - 1;
				idxTo1 = idxFrom1 + 9;
				idxTo2 = idxFrom2 + 9;
				// Скипаем индексы, чтобы не было повторного попадания в локализованный участок
				i = idxTo1;
				j = idxTo2;
				// Цикл уточнения значения корня
				do
				{
					// Переинициализация пределов t1, t2 поиска корня
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

					// Ошиблись с корнем, бывает
					if (dt1 == 0 || dt2 == 0)
					{
						goto wrongRootCase;
					}


					// Двойной цикл поиска ближайшего к реальному корню значения
					p1 = curve1->F(tFrom1);
					p2 = curve2->F(tFrom2);
					curEps = Vector2D(p1.e1, p1.e2, p2.e1, p2.e2).GetLength();
					for (int ii = 0; ii < accuracy; ii++)
					{
						p1 = curve1->F(tFrom1 + ii * dt1);
						for (int jj = 0; jj < accuracy; jj++)
						{
							if (ii == 0 && jj == 0) continue;
							p2 = curve2->F(tFrom2 + jj * dt2);
							newCurEps = Vector2D(p1.e1, p1.e2, p2.e1, p2.e2).GetLength();
							if (newCurEps < curEps)
							{
								curEps = newCurEps;
								idxFrom1 = ii;
								idxFrom2 = jj;
							}
							count++;
						}
					}
					idxFrom1 -= 1;
					idxFrom2 -= 1;
					idxTo1 = idxFrom1 + 9;
					idxTo2 = idxFrom2 + 9;
				} while (curEps > eps);
				// Вносим найденный корень в массив корней
				to_return.push_back((p1 + p2) / 2);

				// Отладкочная информация
				realPoint1 = curve1->F(tFrom1 + dt1 * (idxFrom1 + 1));
				realPoint2 = curve2->F(tFrom2 + dt2 * (idxFrom2 + 1));
				std::cout.precision(12);
				std::cout << "\nТочка пересечения, соответствующая уравнению кривой #1: " << realPoint1.e1 << " " << realPoint1.e2 << "\n";
				std::cout << "Точка пересечения, соответствующая уравнению кривой #2: " << realPoint2.e1 << " " << realPoint2.e2 << "\n";
				std::cout.precision(8);
				std::cout << "Погрешность в решении уравнения f1(" << tFrom1 + dt1 * (idxFrom1 + 1) << ") - f2(" << tFrom2 + dt2 * (idxFrom2 + 1) << ") = 0:    " << (realPoint1 - realPoint2).e1 << " " << (realPoint1 - realPoint2).e2 << "\n\n";
				// Конец отладки

			wrongRootCase:
				// Возвращаем значения переменных для "грубого" поиска корней по кривой
				dt1 = 1.0 / accuracy;
				dt2 = 1.0 / accuracy;
				p1 = curve1->F(i * dt1);
				p2 = curve2->F(j * dt2);
				tFrom1 = 0;
				tFrom2 = 0;
				tTo1 = 1;
				tTo2 = 1;
			}
		}
	}
	std::cout << "Количество итераций: " << count << "\n";
	// Конец поиска корней
	return to_return;
}

// Функция нахождения точного градиента в точке
Vector2D Gradient(Curve* curve1, Curve* curve2, double t1, double t2)
{
	double res1, res2;
	Point2D d1, d2;
	Point2D p1, p2;
	p1 = curve1->F(t1);
	p2 = curve2->F(t2);
	d1 = curve1->dF(t1);
	d2 = curve2->dF(t2);

	res1 = 2 * (d1.e1 * (p1.e1 - p2.e1) + d1.e2 * (p1.e2 - p2.e2));
	res2 = 2 * (d2.e1 * (p2.e1 - p1.e1) + d2.e2 * (p2.e2 - p1.e2));

	return Vector2D(res1, res2);
}

// Функция поиска точки пересечения между массивами точек, представляющих собой 2 различные кривые
DllExport vector<Point2D> FindCrossPointsViaGradient(Curve* curve1, Curve* curve2, double eps = 1e-9, int hypothesisCountOfPoints = 2, double startGradientSpeed = 0.1, double dt = 1e-10)
{
	int count = 0;
	std::cout << "\nРезультаты поиска корней через градиентный спуск:\n";
	// Начальная инициализация
	vector<Point2D> to_return;
	Point2D realPoint1, realPoint2;
	// Функция, которая будет минимизароваться - расстояние между точками кривых
	auto metric = [](Point2D p1, Point2D p2)->double {return Vector2D(p1.e1, p1.e2, p2.e1, p2.e2).GetLength(); };
	// Что-то вроде локализации корней - по-прежнему не совсем понятно, как это сделать
	for (int i = 0; i < hypothesisCountOfPoints; i++)
	{
		double t1 = 1.0 / (hypothesisCountOfPoints + 1) * (i + 1);
		double t2 = 1.0 / (hypothesisCountOfPoints + 1) * (i + 1);
		double a = startGradientSpeed;
		double prevEps, curEps;
		double l1, l2;
		curEps = metric(curve1->F(t1), curve2->F(t2));
		Vector2D next_grad, grad;

		//Vector2D gr;
		// Градиентый спуск к корню
		do
		{
			l1 = metric(curve1->F(t1), curve2->F(t2));
			//grad = gradient(
			//	metric(curve1->F(t1), curve2->F(t2)),
			//	metric(curve1->F(t1 + dt), curve2->F(t2)),
			//	metric(curve1->F(t1), curve2->F(t2 + dt)),
			//	dt, dt);

			grad = Gradient(curve1, curve2, t1, t2);

			//gr = Gradient(curve1, curve2, t1, t2);

			double dt1, dt2;
			dt1 = -grad.e1 * a / grad.GetLength();
			dt2 = -grad.e2 * a / grad.GetLength();

			//next_grad = gradient(
			//	metric(curve1->F(t1 + dt1), curve2->F(t2 + dt2)),
			//	metric(curve1->F(t1 + dt1 + dt), curve2->F(t2 + dt2)),
			//	metric(curve1->F(t1 + dt1), curve2->F(t2 + dt2 + dt)),
			//	dt, dt);

			next_grad = Gradient(curve1, curve2, t1+dt1, t2+dt2);

			// Точка минимума не "перепрыгнута" - делаем шаг побольше
			if (grad * next_grad > 0)
			{
				t1 -= grad.e1 * a / grad.GetLength();
				t2 -= grad.e2 * a / grad.GetLength();
			}
			// Точка минимума "перепрыгнута" - делаем шаг поменьше, и в дальнейшем будем более аккуратны
			else
			{
				a /= 1.25;
				t1 -= grad.e1 * a / grad.GetLength();
				t2 -= grad.e2 * a / grad.GetLength();
			}
			l2 = metric(curve1->F(t1), curve2->F(t2));

			curEps = abs((l2 + l1) / 2);

			count++;
		} while (curEps > eps);
		Point2D p1 = curve1->F(t1);
		Point2D p2 = curve2->F(t2);
		to_return.push_back((p1 + p2) / 2);
		// Отладочная информация
		realPoint1 = curve1->F(t1);
		realPoint2 = curve2->F(t2);
		std::cout.precision(10);
		std::cout << "\nТочка пересечения, соответствующая уравнению кривой #1: " << realPoint1.e1 << " " << realPoint1.e2 << "\n";
		std::cout << "Точка пересечения, соответствующая уравнению кривой #2: " << realPoint2.e1 << " " << realPoint2.e2 << "\n";
		std::cout << "Погрешность f1(" << t1 << ") - f2(" << t2 << ") = 0:   ";
		std::cout.precision(6);
		std::cout << (realPoint1 - realPoint2).e1 << " " << (realPoint1 - realPoint2).e2 << "\n\n";
	iteration_end:;
	}
	std::cout << "Количество итераций: " << count << "\n";
	// Конец поиска корней
	return to_return;
}


