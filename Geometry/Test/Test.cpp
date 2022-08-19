#include <iostream>
#include <fstream>
#include <string>
#include "Curves.h"

int main()
{
	setlocale(LC_ALL, "Russian");

	vector<Point2D> basePoints1;
	vector<Point2D> basePoints2 = vector<Point2D>{ Point2D(2, 0),  Point2D(6, 3),  Point2D(10, 10) };

	std::ifstream reader;
	std::string line;
	reader.open("curve 1 points.txt");
	if (reader.is_open())
	{
		while (!reader.eof())
		{
			double x, y;
			reader >> x;
			reader >> y;
			basePoints1.push_back(Point2D(x, y));
		}
	}
	reader.close();

	reader.open("curve 2 points.txt");

	auto curve1 = Bezier(basePoints1);
	auto curve2 = Bezier(&reader);
	reader.close();

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

	// Проверка на утечку памяти
	for (int i = 0; i < 1000000; i++)
	{
		if (i % 100000 == 0)
		{
			int abc = 1;
		}
		Curve* curve;
		curve = new Bezier(basePoints1);
		delete curve;
	}

	std::cout << "Конец выполнения программы\n";
	char res;
	std::cin >> res;

	return 0;
}