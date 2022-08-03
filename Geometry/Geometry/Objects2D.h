#pragma once
#define DllExport __declspec(dllexport)
#include <math.h>
#include <vector>

using namespace std;

// �����, ����������� ������ � 2D ������������
class DllExport Vector2D
{
public:
	double e1, e2;
	Vector2D()
	{
		e1 = 0;
		e2 = 0;
	}
	// ����������� ����� ���������� �������
	Vector2D(double e1, double  e2)
	{
		this->e1 = e1;
		this->e2 = e2;
	}
	// ����������� ����� ���������� 2 �����: x-y-x-y
	Vector2D(double p1e1, double p1e2, double p2e1, double p2e2)
	{
		this->e1 = p1e1 - p2e1;
		this->e2 = p1e2 - p2e2;
	}
	// �� �����, ������ �� �������� :�
	//Vector2D(Point2D p1, Point2D p2)
	//{
	//	this->e1 = p1.e1 - p2.e1;
	//	this->e2 = p1.e2 - p2.e2;
	//}

	// ���������� �������� ����������
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

// �����, ����������� ����� � 2D ������������
class  DllExport Point2D
{
public:
	double e1, e2;
	Point2D()
	{
		e1 = 0;
		e2 = 0;
	}
	// ����������� ����� 2 ����������
	Point2D(double  e1, double  e2)
	{
		this->e1 = e1;
		this->e2 = e2;
	}

	// ���������� �������� ����������
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
