#pragma once
#include "Objects2D.h"

// Метод поиска градиента функции двух переменных
Vector2D gradient(double f, double fdt1, double fdt2, double dt1, double dt2)
{
	return Vector2D((fdt1 - f) / dt1, (fdt2 - f) / dt2);
}

// Функция определения факториала числа n
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



