#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Objects2D.h"
#include "Funcs.h"

// Класс-родитель для любой кривой
class  DllExport Curve
{
public:
	// Уравнение кривой в параметрическом виде (у каждого наследника будет свое уравнение)
	virtual const Point2D F(double  t) = 0;
	// Уравнение кривой в параметрическом виде (у каждого наследника будет свое уравнение)
	virtual const Point2D dF(double  t) = 0;
	// Получение степени уравнения кривой 
	virtual const int EquationPow() = 0;
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
		int n = basePoints.size() - 1;
		if (t == 0)
			return basePoints[0];
		if (t == 1)
			return basePoints[n];
		Point2D to_return = Point2D(0, 0);
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
		int n = basePoints.size() - 1;
		if (t == 0)
			return basePoints[0] * (-n);
		if (t == 1)
			return basePoints[n] * n;
		Point2D to_return = Point2D(0, 0);
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
	// Получение степени уравнения кривой 
	const int EquationPow() override
	{
		return basePoints.size();
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

// Функция нахождения чисто математического градиента 
const Vector2D gradient(Curve* curve1, Curve* curve2, double t1, double t2)
{
	double gr_e1, gr_e2;
	Point2D d1, d2;
	Point2D p1, p2;
	p1 = curve1->F(t1);
	p2 = curve2->F(t2);
	d1 = curve1->dF(t1);
	d2 = curve2->dF(t2);

	gr_e1 = 2 * (d1.e1 * (p1.e1 - p2.e1) + d1.e2 * (p1.e2 - p2.e2));
	gr_e2 = 2 * (d2.e1 * (p2.e1 - p1.e1) + d2.e2 * (p2.e2 - p1.e2));

	return Vector2D(gr_e1, gr_e2);
}

/// Функция поиска наименьшего расстояния между массивами точек, представляющих собой 2 различные кривые.
/// Возвращаемое значение - массив из двух точек(на 1й и 2й кривой соответственно)
DllExport const vector<Point2D> FindClosestPoints(Curve* curve1, Curve* curve2, double eps = 1e-9, bool debug = false)
{
	if (curve1 == NULL)
	{
		throw std::invalid_argument("Pointer to the 1st curve is NULL");
		return vector<Point2D>();
	}
	if (curve2 == NULL)
	{
		throw std::invalid_argument("Pointer to the 2nd curve is NULL");
		return vector<Point2D>();
	}
	if (debug)
		std::cout << "\nРезультаты поиска ближайших точек через градиентный спуск:\n";
	// Начальная инициализация
	vector<Point2D> to_return;
	Point2D realPoint1, realPoint2;
	// Функция, которая будет минимизароваться - расстояние между точками кривых
	auto metric = [](Point2D p1, Point2D p2)->double {return Vector2D(p1.e1, p1.e2, p2.e1, p2.e2).GetLength(); };
	double t1 = 0.5;
	double t2 = 0.5;
	double dt1, dt2;
	double a = 0.1;
	double prevEps, curEps;
	double l1, l2;
	curEps = metric(curve1->F(t1), curve2->F(t2));
	Vector2D next_grad, grad;

	// Градиентый спуск к корню
	do
	{
		l1 = metric(curve1->F(t1), curve2->F(t2));
		grad = gradient(curve1, curve2, t1, t2);

		if (grad.GetLength() == 0)
		{
			throw std::invalid_argument("Math error: Attempted to divide by zero");
			break;
		}
		else
		{
			dt1 = -grad.e1 * a / grad.GetLength();
			dt2 = -grad.e2 * a / grad.GetLength();
		}

		next_grad = gradient(curve1, curve2, t1 + dt1, t2 + dt2);

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

		curEps = abs(l2 - l1);

	} while (curEps > eps);
	Point2D p1 = curve1->F(t1);
	Point2D p2 = curve2->F(t2);
	to_return.push_back((p1 + p2) / 2);

	// Отладочная информация
	if (debug)
	{
		realPoint1 = curve1->F(t1);
		realPoint2 = curve2->F(t2);
		std::cout.precision(10);
		std::cout << "\nБлижайшая точка, соответствующая уравнению кривой #1: " << realPoint1.e1 << " " << realPoint1.e2 << "\n";
		std::cout << "Ближайшая точка, соответствующая уравнению кривой #2: " << realPoint2.e1 << " " << realPoint2.e2 << "\n";
		std::cout << "Наименьшее расстояние между кривыми равняется:        " << l2 << "\n";
	}
	// Конец поиска корней
	return to_return;
}

/// Функция поиска точки пересечения между двумя различными кривыми
/// curve1 - указатель на первую кривую
/// curve2 - указатель на вторую кривую
/// eps - необходимая точность нахождения точек
/// Возвращаемое значение - массив всех найденных точек пересечения
DllExport const vector<Point2D> FindCrossPoints(Curve* curve1, Curve* curve2, double eps = 1e-9, bool debug = false)
{
	if (curve1 == NULL)
	{
		throw std::invalid_argument("Pointer to the 1st curve is NULL");
		return vector<Point2D>();
	}
	if (curve2 == NULL)
	{
		throw std::invalid_argument("Pointer to the 2nd curve is NULL");
		return vector<Point2D>();
	}
	if (debug)
	{
		std::cout << "\nРезультаты поиска корней через градиентный спуск:\n";
	}
	// Начальная инициализация
	vector<Point2D> to_return;
	Point2D realPoint1, realPoint2;
	// Функция, которая будет минимизароваться - расстояние между точками кривых
	auto metric = [](Point2D p1, Point2D p2)->double {return Vector2D(p1.e1, p1.e2, p2.e1, p2.e2).GetLength(); };
	// Количество точек пересечения не может превышать значение максимальной степени уравнения - 1 (если только кривые не совпадают)
	int hypothesisCountOfPoints = curve1->EquationPow() > curve2->EquationPow() ? curve1->EquationPow() - 1 : curve2->EquationPow() - 1;
	// Что-то вроде локализации корней - по-прежнему не совсем понятно, как это сделать
	for (int i = 0; i < hypothesisCountOfPoints; i++)
	{
		double t1 = 1.0 / (hypothesisCountOfPoints + 1) * (i + 1);
		double t2 = 1.0 / (hypothesisCountOfPoints + 1) * (i + 1);
		double dt1, dt2;
		double a = 0.1;
		double prevEps, curEps;
		double l1, l2;
		curEps = metric(curve1->F(t1), curve2->F(t2));
		Vector2D next_grad, grad;

		// Градиентый спуск к корню
		do
		{
			l1 = metric(curve1->F(t1), curve2->F(t2));
			grad = gradient(curve1, curve2, t1, t2);

			if (grad.GetLength() == 0)
			{
				throw std::invalid_argument("Math error: Attempted to divide by zero");
				break;
			}
			else
			{
				dt1 = -grad.e1 * a / grad.GetLength();
				dt2 = -grad.e2 * a / grad.GetLength();
			}

			next_grad = gradient(curve1, curve2, t1 + dt1, t2 + dt2);

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

		} while (curEps > eps);
		Point2D p1 = curve1->F(t1);
		Point2D p2 = curve2->F(t2);
		Point2D avgPoint = (p1 + p2) / 2;

		// Проверка на уникальность корня
		bool unique = true;
		for (auto point : to_return)
		{
			if (EqualPoints(avgPoint, point, eps))
			{
				unique = false;
				break;
			}
		}
		if (unique)
		{
			to_return.push_back(avgPoint);

			// Отладочная информация
			if (debug)
			{
				realPoint1 = curve1->F(t1);
				realPoint2 = curve2->F(t2);
				std::cout.precision(10);
				std::cout << "\nТочка пересечения, соответствующая уравнению кривой #1: " << realPoint1.e1 << " " << realPoint1.e2 << "\n";
				std::cout << "Точка пересечения, соответствующая уравнению кривой #2: " << realPoint2.e1 << " " << realPoint2.e2 << "\n";
				std::cout << "Погрешность f1(" << t1 << ") - f2(" << t2 << ") = 0:   ";
				std::cout.precision(6);
				std::cout << (realPoint1 - realPoint2).e1 << " " << (realPoint1 - realPoint2).e2 << "\n\n";
			}
		}
	}
	// Конец поиска корней
	return to_return;
}


