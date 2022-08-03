#pragma once
#include <vector>
#include "Objects2D.h"
#include "Funcs.h"

using namespace std;

// Класс-родитель для любой кривой
class  DllExport Curve
{
protected:
public:
	// Точность моделирования
	int accuracy;
	// Уравнение кривой в параметрическом виде (у каждого наследника будет свое уравнение)
	virtual Point2D MainFunc(double  t) = 0;
	// Получение массива точек (Point2D *) кривой
	virtual vector<Point2D> GetCurvePoints() = 0;
	// Получение двумерного массива координат точек
	virtual vector<vector<double>> GetCurveCoords() = 0;
	// Метод увеличения точности моделирования кривой 
	virtual vector<Point2D> ImproveAccuracy(double t1, double t2, int accuracy) = 0;
};

// Класс, представляющий кривую Безье 2го порядка
class  DllExport Bezier : public Curve
{
	// Массив точек
	vector<Point2D> points;
	// Опорные точки
	vector<Point2D> basePoints;
public:
	// Три опорные точки и точность моделирования (по умолчанию 100)
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
	// Уравнение кривой в параметрическом виде 
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
	// Получение массива точек (Point2D *) кривой
	vector<Point2D> GetCurvePoints()override
	{
		return points;
	}
	// Получение двумерного массива координат точек
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
	// Метод увеличения точности моделирования кривой 
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

