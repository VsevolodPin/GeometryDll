// Test.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Objects2D.h"
#include "Curves.h"

int main()
{
	setlocale(LC_ALL, "Russian");

	vector<Point2D> basePoints1 = vector<Point2D>{ Point2D(3, 0),  Point2D(1, 3),  Point2D(5, 5),  Point2D(7, 3) };
	vector<Point2D> basePoints2 = vector<Point2D>{ Point2D(2, 0),  Point2D(6, 3),  Point2D(10, 10) };

	auto curve1 = Bezier(basePoints1);
	auto curve2 = Bezier(basePoints2);

	Curve* p1 = &curve1;
	Curve* p2 = &curve2;

	vector<Point2D> crossPoint;
	vector<Point2D> closestPoints;

	closestPoints = FindClosestPoints(p1, p2, 1e-9);
	std::cout.precision(9);
	std::cout << "Точки с наименьшим расстоянием между кривыми, найденные через функцию без шаблона:\n" <<
		closestPoints[0].e1 << " " <<
		closestPoints[0].e2 << " " <<
		closestPoints[1].e1 << " " <<
		closestPoints[1].e2 << "\n";

	crossPoint = FindCrossPointsViaSegments(p1, p2, 1e-9);

	crossPoint = FindCrossPointsViaEquations(p1, p2, 1e-9, 100);

	crossPoint = FindCrossPointsViaGradient(p1, p2, 1e-9, 2);

	std::cout << "Конец выполнения программы\n";
	char res;
	std::cin >> res;
}