#pragma once
#define DllExport __declspec(dllexport)
#include <math.h>
#include <vector>

using namespace std;

//  ласс, описывающий вектор в 2D пространстве
class DllExport Vector2D
{
public:
	double e1, e2;
	Vector2D()
	{
		e1 = 0;
		e2 = 0;
	}
	//  онструктор через компоненты вектора
	Vector2D(double e1, double  e2)
	{
		this->e1 = e1;
		this->e2 = e2;
	}
	//  онструктор через координаты 2 точек: x-y-x-y
	Vector2D(double p1e1, double p1e2, double p2e1, double p2e2)
	{
		this->e1 = p1e1 - p2e1;
		this->e2 = p1e2 - p2e2;
	}
	// Ќе пон€л, почему не работает :с
	//Vector2D(Point2D p1, Point2D p2)
	//{
	//	this->e1 = p1.e1 - p2.e1;
	//	this->e2 = p1.e2 - p2.e2;
	//}

	// ѕерегрузка основных операторов
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

//  ласс, описывающий точку в 2D пространстве
class  DllExport Point2D
{
public:
	double e1, e2;
	Point2D()
	{
		e1 = 0;
		e2 = 0;
	}
	//  онструктор через 2 координаты
	Point2D(double  e1, double  e2)
	{
		this->e1 = e1;
		this->e2 = e2;
	}

	// ѕерегрузка основных операторов
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

//  ласс-родитель дл€ любой кривой
class  DllExport Curve
{
protected:
	// ћассив точек
	vector<Point2D> points;
public:
	// “очность моделировани€
	int accuracy;
	// ”равнение кривой в параметрическом виде (у каждого наследника будет свое уравнение)
	virtual Point2D MainFunc(double  t) = 0;
	// ѕолучение массива точек (Point2D *) кривой
	virtual vector<Point2D> GetCurvePoints() = 0;
	// ѕолучение двумерного массива координат точек
	virtual vector<vector<double>> GetCurveCoords() = 0;
	// ћетод увеличени€ точности моделировани€ кривой 
	virtual vector<Point2D> ImproveAccuracy(double t1, double t2, int accuracy) = 0;
};

//  ласс, представл€ющий кривую Ѕезье 2го пор€дка
class  DllExport Bezier : public Curve
{
	Point2D p1, p2, p3;
public:
	// “ри опорные точки и точность моделировани€ (по умолчанию 100)
	Bezier(Point2D p1, Point2D p2, Point2D p3, int accuracy = 100)
	{
		this->p1 = p1;
		this->p2 = p2;
		this->p3 = p3;
		this->accuracy = accuracy;
		points.resize(accuracy);
		for (int i = 0; i < accuracy; i++)
		{
			points[i] = MainFunc(1.0 / (accuracy - 1) * i);
		}
	}
	~Bezier()
	{
		points.clear();
	}
	// ”равнение кривой в параметрическом виде 
	Point2D MainFunc(double  t) override
	{
		return p1 * (1 - t) * (1 - t) + p2 * 2 * t * (1 - t) + p3 * t * t;
	}
	// ѕолучение массива точек (Point2D *) кривой
	vector<Point2D> GetCurvePoints()override
	{
		return points;
	}
	// ѕолучение двумерного массива координат точек
	vector<vector<double>> GetCurveCoords()override
	{
		vector<vector<double>> coords;
		coords.resize(2);
		coords[0].resize(accuracy);
		coords[1].resize(accuracy);
		for (int i = 0; i < accuracy; i++)
		{
			coords[0][i] = points[i].e1;
			coords[1][i] = points[i].e2;
		}
		return coords;
	}
	// ћетод увеличени€ точности моделировани€ кривой 
	vector<Point2D> ImproveAccuracy(double t1, double t2, int accuracy)override
	{
		vector<Point2D> to_return;
		to_return.resize(accuracy);
		for (int i = 0; i < accuracy; i++)
		{
			to_return[i] = MainFunc(t1 + (t2 - t1) / (accuracy - 1) * i);
		}
		return to_return;
	}
};

//class DllExport Hermite : public Curve
//{
//	Point2D p1, p2;
//public:
//	// ƒве опорные точки и точность моделировани€ (по умолчанию 100)
//	Hermite(Point2D p1, Point2D p2, int accuracy = 100)
//	{
//		this->p1 = p1;
//		this->p2 = p2;
//		this->accuracy = accuracy;
//		points = new Point2D[accuracy];
//		for (int i = 0; i < accuracy; i++)
//		{
//			points[i] = MainFunc(1.0 / (accuracy - 1) * i);
//		}
//	}
//	~Hermite()
//	{
//		delete[] points;
//	}
//	// ”равнение кривой в параметрическом виде 
//	Point2D MainFunc(double  t) override
//	{
//		Point2D m1 = ((p2 - p1) / (p2.e1 - p1.e1) + (p1 / p1.e1)) * 0.5;
//		Point2D m2 = (p2 / p2.e1) + ((p2 - p1) / (p2.e1 - p1.e1)) * 0.5;
//		return p1 * (2 * t * t * t - 3 * t * t + 1) + m1 * (t * t * t - 2 * t * t - 1) * (p2.e1 - p1.e1) + p2 * (-2 * t * t * t + 3 * t * t) + m2 * (t * t * t - t * t) * (p2.e1 - p1.e1);
//	}
//	// ѕолучение массива точек (Point2D *) кривой
//	Point2D* GetCurvePoints()override
//	{
//		return points;
//	}
//	// ѕолучение двумерного массива координат точек
//	double** GetCurveCoords()override
//	{
//		double** coords = new double* [2];
//		coords[0] = new double[accuracy];
//		coords[1] = new double[accuracy];
//		for (int i = 0; i < accuracy; i++)
//		{
//			coords[0][i] = points[i].e1;
//			coords[1][i] = points[i].e2;
//		}
//		return coords;
//	}
//};

// ‘ункци€ поиска наименьшего рассто€ни€ между массивами точек, представл€ющих собой 2 различные кривые
template<typename CurveType1, typename CurveType2> DllExport vector<Point2D> FindClosestPoints(CurveType1 curve1, CurveType2 curve2, double eps = 1e-9)
{
	vector<Point2D> to_return;
	to_return.resize(2);
	vector<Point2D> closestPoints;
	closestPoints.resize(2);
	int idxP1 = 0;
	int idxP2 = 0;
	double tFrom1 = 0, tFrom2 = 0;
	double tTo1 = 1, tTo2 = 1;
	auto p1 = curve1.GetCurvePoints();
	auto p2 = curve2.GetCurvePoints();
	to_return[0] = p1[0];
	to_return[1] = p2[0];
	do
	{
		Vector2D v = Vector2D(p1[0].e1, p1[0].e2, p2[0].e1, p2[0].e2);
		double  min = v.GetLength();

		for (int i = 0; i < curve1.accuracy; i++)
		{
			for (int j = 0; j < curve2.accuracy; j++)
			{
				v = Vector2D(p1[i].e1, p1[i].e2, p2[j].e1, p2[j].e2);
				auto dist = v.GetLength();
				if (dist < min)
				{
					idxP1 = i;
					idxP2 = j;
					closestPoints[0] = p1[i];
					closestPoints[1] = p2[j];
					min = dist;
				}
			}
		}

		double curEps = abs(Vector2D(closestPoints[0].e1, closestPoints[0].e2, closestPoints[1].e1, closestPoints[1].e2).GetLength() - 
			Vector2D(to_return[0].e1, to_return[0].e2, to_return[1].e1, to_return[1].e2).GetLength());
		if (curEps < eps)
		{
			to_return[0] = closestPoints[0];
			to_return[1] = closestPoints[1];
			return to_return;
		}
		else
		{
			to_return[0] = closestPoints[0];
			to_return[1] = closestPoints[1];
			if (idxP1 == 0) idxP1++;
			if (idxP1 == curve1.accuracy - 1) idxP1--;
			if (idxP2 == 0) idxP2++;
			if (idxP2 == curve2.accuracy - 1) idxP2--;
			double tFromNew1 = ((idxP1 - 1) / (double)(curve1.accuracy-1)) * (tTo1 - tFrom1);
			double tToNew1 = ((idxP1 + 1) / (double)(curve1.accuracy-1)) * (tTo1 - tFrom1);
			double tFromNew2 = ((idxP2 - 1) / (double)(curve2.accuracy-1)) * (tTo2 - tFrom2);
			double tToNew2 = ((idxP2 + 1) / (double)(curve2.accuracy-1)) * (tTo2 - tFrom2);
			p1 = curve1.ImproveAccuracy(tFromNew1, tToNew1, curve1.accuracy);
			p2 = curve2.ImproveAccuracy(tFromNew2, tToNew2, curve2.accuracy);
		}

	} while (true);
	return to_return;
}
