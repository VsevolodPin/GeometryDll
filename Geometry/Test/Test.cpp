// Test.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Object2D.h"

int main()
{
	auto curve1 = Bezier(Point2D(5, 5), Point2D(1, 9), Point2D(15, 15), 100);
	auto curve2 = Bezier(Point2D(1, 0), Point2D(10, 1), Point2D(11, 10), 100);
	auto test = curve1.GetCurveCoords();

	auto closestPoints = FindClosestPoints<Bezier, Bezier>(curve1, curve2, 1e-9);
	std::cout.precision(9);
	std::cout << "Min distance between curve1 and curve2 = " <<
		closestPoints[0].e1 << " " <<
		closestPoints[0].e2 << " " <<
		closestPoints[1].e1 << " " <<
		closestPoints[1].e2 << "\n";
}