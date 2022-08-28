#include <iostream>
#include <fstream>
#include <string>
#include "Curves.h"

int main()
{
	setlocale(LC_ALL, "Russian");

	vector<Point2D> basePoints1;
	vector<Point2D> crossPoint;
	vector<Point2D> closestPoints;

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
	crossPoint = FindCrossPoints(p1, p2, 1e-9, true);

	reader.open("curve 3 points.txt");
	auto curve3 = Bezier(&reader, true);
	reader.open("curve 4 points.txt");
	auto curve4 = Bezier(&reader, true);

	p1 = &curve3;
	p2 = &curve4;
	closestPoints = FindClosestPoints(p1, p2, 1e-9, true);

	std::cout << "Конец выполнения программы\n";
	char res;
	std::cin >> res;

	return 0;
}