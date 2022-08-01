#pragma once
#define DllExport __declspec(dllexport)
#include <math.h>


class DllExport Vector2D 
{
public:
	double e1, e2;
	Vector2D()
	{
		e1 = 0;
		e2 = 0;
	}
	Vector2D(double e1, double  e2)
	{
		this->e1 = e1;
		this->e2 = e2;
	}
	Vector2D(double p1e1, double p1e2, double p2e1, double p2e2)
	{
		this->e1 = p1e1 - p2e1;
		this->e2 = p1e2 - p2e2;
	}
	//Vector2D(Point2D p1, Point2D p2)
	//{
	//	this->e1 = p1.e1 - p2.e1;
	//	this->e2 = p1.e2 - p2.e2;
	//}


	Vector2D operator + (Vector2D v)
	{
		return Vector2D(e1 + v.e1, e2 + v.e2);
	}
	Vector2D operator -(Vector2D v)
	{
		return Vector2D(e1 - v.e1, e2 - v.e2);
	}
	Vector2D operator * (double  v)
	{
		return Vector2D(e1 * v, e2 * v);
	}
	double  operator * (Vector2D v)
	{
		return e1 * v.e1 + e2 * v.e2;
	}
	double  GetLength()
	{
		return sqrt(e1 * e1 + e2 * e2);
	}
};

class  DllExport Point2D 
{
public:
	double e1, e2;
	Point2D()
	{
		e1 = 0;
		e2 = 0;
	}
	Point2D(double  e1, double  e2)
	{
		this->e1 = e1;
		this->e2 = e2;
	}


	Point2D operator + (Point2D v)
	{
		return Point2D(e1 + v.e1, e2 + v.e2);
	}
	Point2D operator -(Point2D v)
	{
		return Point2D(e1 - v.e1, e2 - v.e2);
	}
	Point2D operator * (double  v)
	{
		return Point2D(e1 * v, e2 * v);
	}
	Point2D operator / (double  v)
	{
		return Point2D(e1 / v, e2 / v);
	}
};

class  DllExport Curve
{
protected:
	Point2D* arr;
public:
	int accuracy;
	virtual Point2D MainFunc(double  t) = 0;
	virtual Point2D* GetCurvePoints() = 0;
	virtual double** GetCurveCoords() = 0;
};

class  DllExport Bezier : public Curve
{
	Point2D p1, p2, p3;
public:
	Bezier(Point2D p1, Point2D p2, Point2D p3, int accuracy = 10)
	{
		this->p1 = p1;
		this->p2 = p2;
		this->p3 = p3;
		this->accuracy = accuracy;
		arr = new Point2D[accuracy];
		for (int i = 0; i < accuracy; i++)
		{
			arr[i] = MainFunc(1.0 / (accuracy - 1) * i);
		}
	}
	~Bezier()
	{
		delete[] arr;
	}
	Point2D MainFunc(double  t) override
	{
		return p1 * (1 - t) * (1 - t) + p2 * 2 * t * (1 - t) + p3 * t * t;
	}
	Point2D* GetCurvePoints()override
	{
		return arr;
	}
	double** GetCurveCoords()override
	{
		double** coords = new double* [2];
		coords[0] = new double[accuracy];
		coords[1] = new double[accuracy];
		for (int i = 0; i < accuracy; i++)
		{
			coords[0][i] = arr[i].e1;
			coords[1][i] = arr[i].e2;
		}
		return coords;
	}
};

DllExport double  FindMinDistance(Point2D* p1, Point2D* p2, int size1, int size2)
{
	Vector2D v = Vector2D(p1[0].e1, p1[0].e2, p2[0].e1, p2[0].e2);
	double  min = v.GetLength();

	for (int i = 1; i < size1; i++)
	{
		for (int j = 0; j < size2; j++)
		{
			v = Vector2D(p1[0].e1, p1[0].e2, p2[0].e1, p2[0].e2);
			auto dist = v.GetLength();
			if (dist < min)
				min = dist;
		}
	}
	return min;
}
