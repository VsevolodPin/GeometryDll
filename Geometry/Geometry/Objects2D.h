#pragma once
#define DllExport __declspec(dllexport)
#include <math.h>
#include <vector>

using namespace std;

// Класс, описывающий вектор в 2D пространстве
class DllExport Vector2D
{
public:
	double e1, e2;
	Vector2D()
	{
		e1 = 0;
		e2 = 0;
	}
	// Конструктор через компоненты вектора
	Vector2D(double e1, double  e2, bool normalize = false)
	{
		this->e1 = e1;
		this->e2 = e2;
		if (normalize)
			this->Normalize();
	}
	// Конструктор через координаты 2 точек: x-y-x-y
	Vector2D(double p1e1, double p1e2, double p2e1, double p2e2, bool normalize = false)
	{
		this->e1 = p1e1 - p2e1;
		this->e2 = p1e2 - p2e2;
		if (normalize)
			this->Normalize();
	}
	// Нормализация вектора
	void Normalize()
	{
		double l = this->GetLength();
		if (l != 0)
		{
			this->e1 /= l;
			this->e2 /= l;
		}
		else
		{
			throw std::invalid_argument("Math error: Attempted to divide by zero");
		}
	}

	// Перегрузка основных операторов
	Vector2D operator + (Vector2D v) const
	{
		return Vector2D(e1 + v.e1, e2 + v.e2);
	}
	Vector2D operator -(Vector2D v) const
	{
		return Vector2D(e1 - v.e1, e2 - v.e2);
	}
	Vector2D operator * (double  v) const
	{
		return Vector2D(e1 * v, e2 * v);
	}
	double operator * (Vector2D v) const
	{
		return e1 * v.e1 + e2 * v.e2;
	}
	double GetLength() const
	{
		return sqrt(e1 * e1 + e2 * e2);
	}
};

// Класс, описывающий точку в 2D пространстве
class  DllExport Point2D
{
public:
	double e1, e2;
	Point2D()
	{
		e1 = 0;
		e2 = 0;
	}
	// Конструктор через 2 координаты
	Point2D(double  e1, double  e2)
	{
		this->e1 = e1;
		this->e2 = e2;
	}

	// Перегрузка основных операторов
	Point2D operator + (Point2D v) const
	{
		return Point2D(e1 + v.e1, e2 + v.e2);
	}
	Point2D operator -(Point2D v) const
	{
		return Point2D(e1 - v.e1, e2 - v.e2);
	}
	Point2D  operator  * (double  v) const
	{
		return Point2D(e1 * v, e2 * v);
	}
	Point2D operator / (double  v) const
	{
		if (v == 0)
		{
			throw std::invalid_argument("Math error: Attempted to divide by zero");
			return  Point2D(e1, e2);
		}
		return Point2D(e1 / v, e2 / v);
	}
	Point2D operator * (Point2D  v) const
	{
		return Point2D(e1 * v.e1, e2 * v.e2);
	}
	bool operator == (Point2D  v1) const
	{
		if (this->e1 == v1.e1 && this->e2 == v1.e2)
			return true;
		else return false;
	}
};

