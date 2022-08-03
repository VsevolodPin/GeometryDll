using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geom
{
    public class Point2D
    {
        public double e1, e2;
        public Point2D()
        {
            e1 = 0;
            e2 = 0;
        }
        public Point2D(double e1, double e2)
        {
            this.e1 = e1;
            this.e2 = e2;
        }


        public static Point2D operator +(Point2D p1, Point2D p2)
        {
            return new Point2D(p1.e1 + p2.e1, p1.e2 + p2.e2);
        }
        public static Point2D operator -(Point2D p1, Point2D p2)
        {
            return new Point2D(p1.e1 - p2.e1, p1.e2 - p2.e2);
        }
        public static Point2D operator *(Point2D p, double a)
        {
            return new Point2D(p.e1 * a, p.e2 * a);
        }
        public static Point2D operator /(Point2D p, double a)
        {
            return new Point2D(p.e1 / a, p.e2 / a);
        }
    }

    public class Vector2D
    {
        double e1, e2;
        public Vector2D()
        {
            e1 = 0;
            e2 = 0;
        }
        public Vector2D(double e1, double e2)
        {
            this.e1 = e1;
            this.e2 = e2;
        }

        public Vector2D(Point2D p1, Point2D p2)
        {
            this.e1 = p1.e1 - p2.e1;
            this.e2 = p1.e2 - p2.e2;
        }

        public static Vector2D operator +(Vector2D v1, Vector2D v2)
        {
            return new Vector2D(v1.e1 + v2.e1, v1.e2 + v2.e2);
        }
        public static Vector2D operator -(Vector2D v1, Vector2D v2)
        {
            return new Vector2D(v1.e1 - v2.e1, v1.e2 - v2.e2);
        }
        public static Vector2D operator *(Vector2D v, double a)
        {
            return new Vector2D(v.e1 * a, v.e2 * a);
        }
        public static double operator *(Vector2D v1, Vector2D v2)
        {
            return v1.e1 * v2.e1 + v1.e2 * v2.e2;
        }
        public double GetLength()
        {
            return Math.Sqrt(e1 * e1 + e2 * e2);
        }
    };

    public abstract class Curve
    {
        public int accuracy;
        public List<Point2D> points;
        public abstract Point2D MainFunc(double t);
        public abstract List<Point2D> ImproveAccuracy(double t1, double t2, int accuracy);
    }

    public class Bezier : Curve
    {
        public Point2D p1, p2, p3;
        public Point2D[] basePoints;

        public Bezier(Point2D p1, Point2D p2, Point2D p3, int accuracy = 10)
        {
            this.p1 = p1;
            this.p2 = p2;
            this.p3 = p3;
            this.accuracy = accuracy;
            points = new List<Point2D>();
            for (int i = 0; i < accuracy; i++)
            {
                var res = MainFunc(1.0 / (accuracy - 1) * i);
                points.Add(res);
            }
        }
        public Bezier(Point2D[] basePoints, int accuracy = 10)
        {
            this.basePoints = basePoints;
            this.accuracy = accuracy;
            points = new List<Point2D>();
            for (int i = 0; i < accuracy; i++)
            {
                var res = MainFunc2(1.0 / (accuracy - 1) * i);
                points.Add(res);
            }
        }

        public override Point2D MainFunc(double t)
        {
            var P1 = p1 * (1 - t) * (1 - t);
            var P2 = p2 * 2 * t * (1 - t);
            var P3 = p3 * t * t;
            var to_return = P1 + P2 + P3;
            return to_return;
        }

        public int fact(int n)
        {
            if (n <= 1) return 1;
            int to_return = 1;
            for (int i = n; i > 0; i--)
            {
                to_return *= i;
            }
            return to_return;
        }

        public Point2D MainFunc2(double t)
        {
            Point2D to_return = new Point2D(0, 0);
            int n = basePoints.Length - 1;
            for (int i = 0; i < basePoints.Length; i++)
            {
                if (i != 0 && n - i != 0)
                {
                    double N = fact(n);
                    double K = fact(i);
                    double N_K = fact(n - i);
                    to_return = to_return + basePoints[i] * Math.Pow(t, i) * Math.Pow(1 - t, n - i) * (double)(fact(n)) / (fact(i) * fact(n - i));
                }
                else
                {
                    if (i == 0)
                        to_return = to_return + basePoints[i] * 1 * Math.Pow(1 - t, n - i) * (double)(fact(n)) / (fact(i) * fact(n - i));
                    if (n - i == 0)
                        to_return = to_return + basePoints[i] * Math.Pow(t, i) * 1 * (double)(fact(n)) / (fact(i) * fact(n - i));
                }

            }
            return to_return;
        }

        public List<Point2D> GetCurvePoints()
        {
            return points;
        }
        public double[][] GetCurveCoords()
        {
            double[][] coords = new double[2][];
            coords[0] = new double[accuracy];
            coords[1] = new double[accuracy];
            for (int i = 0; i < accuracy; i++)
            {
                coords[0][i] = points[i].e1;
                coords[1][i] = points[i].e2;
            }
            return coords;
        }
        public override List<Point2D> ImproveAccuracy(double t1, double t2, int accuracy)
        {
            List<Point2D> to_return = new List<Point2D>();
            for (int i = 0; i < accuracy; i++)
            {
                to_return.Add(MainFunc(t1 + (t2 - t1) / (accuracy - 1) * i));
            }
            return to_return;
        }
    }
}

