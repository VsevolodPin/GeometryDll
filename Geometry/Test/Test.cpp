// Test.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Object2D.h"

int main()
{
	auto curve1 = Bezier(Point2D(0, 0), Point2D(1, 9), Point2D(10, 10), 1000);
	auto curve2 = Bezier(Point2D(1, 0), Point2D(10, 1), Point2D(11, 10), 1000);

	auto minDist = FindMinDistance(curve1.GetCurvePoints(), curve2.GetCurvePoints(), 1000, 1000);

	std::cout << "Min distance between curve1 and curve2 = " << minDist << "\n";
}