#pragma once
#include "Objects2D.h"

// ‘ункци€ определени€ факториала числа n
static int fact(int n)
{
	if (n <= 1) return 1;
	int to_return = 1;
	for (int i = n; i > 0; i--)
	{
		to_return *= i;
	}
	return to_return;
}

// ‘ункци€, позвол€юща€ сравнивать точки с определенной точностью
static bool EqualPoints(Point2D p1, Point2D p2, double eps = 1e-9)
{
	int precision = -log10(eps);
	int mult = 1;
	for (int i = 0; i < precision; i++)
		mult *= 10;

	if ((int)(p1.e1 * mult) == (int)(p2.e1 * mult) &&
		(int)(p1.e2 * mult) == (int)(p2.e2 * mult))
		return true;
	else return false;
}


