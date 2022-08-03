// Test.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Objects2D.h"
#include "Curves.h"

int main()
{
	setlocale(LC_ALL, "Russian");

	vector<Point2D> basePoints1 = vector<Point2D>{ Point2D(10, 10), Point2D(11, 11), Point2D(12, 13), Point2D(15, 15) };
	vector<Point2D> basePoints2 = vector<Point2D>{ Point2D(1, 15), Point2D(7, 12), Point2D(13, 9) };

	auto curve1 = Bezier(basePoints1, 100);
	auto curve2 = Bezier(basePoints2, 100);
	auto test = curve1.GetCurveCoords();

	auto closestPoints = FindClosestPoints<Bezier, Bezier>(curve1, curve2, 1e-9);
	std::cout.precision(9);
	std::cout << "Точки с наименьшим расстоянием между кривымы: " <<
		closestPoints[0].e1 << " " <<
		closestPoints[0].e2 << " " <<
		closestPoints[1].e1 << " " <<
		closestPoints[1].e2 << "\n";

	auto crossPoint = FindCrossPoints<Bezier, Bezier>(curve1, curve2, 1e-9);
	std::cout << "Точки с наименьшим расстоянием между кривымы: " <<
		crossPoint[0].e1 << " " <<
		crossPoint[0].e2 << "\n";
}