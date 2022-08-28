using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;


namespace Geom
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        public double FindMinDistance(List<Point2D> p1, List<Point2D> p2, out int idxP1, out int idxP2)
        {
            double min = new Vector2D(p1[0], p2[0]).GetLength();
            idxP1 = 0;
            idxP2 = 0;

            for (int i = 0; i < p1.Count; i++)
            {
                for (int j = 0; j < p1.Count; j++)
                {
                    var dist = new Vector2D(p1[i], p2[j]).GetLength();
                    if (dist < min)
                    {
                        idxP1 = i;
                        idxP2 = j;
                        min = dist;
                    }
                }
            }
            return min;
        }

        public Point2D FindCrossPoint(Curve curve1, Curve curve2, double eps)
        {
            double prevLength = double.MaxValue, curLength = 0;
            double curEps = 0.1;
            List<Point2D> arr1 = curve1.points;
            List<Point2D> arr2 = curve2.points;
            do
            {
                int idxP1, idxP2;
                curLength = FindMinDistance(arr1, arr2, out idxP1, out idxP2);
                if (idxP1 == 0) idxP1++;
                if (idxP2 == 0) idxP2++;
                if (idxP1 == arr1.Count - 1) idxP1--;
                if (idxP2 == arr2.Count - 1) idxP2--;

                if (curLength < prevLength)
                {
                    if ((curve1.points[idxP1 - 1].e2 - curve2.points[idxP2 - 1].e2) * (curve1.points[idxP1 + 1].e2 - curve2.points[idxP2 + 1].e2) < 0)
                    {
                        curEps = prevLength - curLength;
                        if (curEps < 0) curEps *= -1;
                        if (curEps < eps)
                        {
                            return (curve1.points[idxP1] + curve2.points[idxP2]) / 2;
                        }
                        else
                        {
                            prevLength = curLength;
                            arr1 = curve1.ImproveAccuracy(((double)(idxP1 - 1)) / curve1.accuracy, ((double)(idxP1 + 1)) / curve1.accuracy, 100);
                            arr2 = curve2.ImproveAccuracy(((double)(idxP2 - 1)) / curve2.accuracy, ((double)(idxP2 + 1)) / curve2.accuracy, 100);
                        }
                    }
                }
            } while (true);
            return new Point2D(0, 0);
        }

        private void Button1_Click(object sender, EventArgs e)
        {
            //var curve1 = new Bezier(new Point2D(5, 5), new Point2D(1, 9), new Point2D(15, 15), 100);
            //var curve2 = new Bezier(new Point2D(1, 0), new Point2D(10, 1), new Point2D(11, 10), 100);
            var curve1 = new Bezier(new Point2D[] { new Point2D(3, 0), new Point2D(1, 3), new Point2D(5, 5), new Point2D(7, 3) }, 100);
            var curve2 = new Bezier(new Point2D[] { new Point2D(2, 0), new Point2D(6, 3), new Point2D(10, 10) }, 100);

            var curve3 = new Bezier(new Point2D[] { new Point2D(0, 4), new Point2D(4, 6), new Point2D(7, 7), new Point2D(4, 10) }, 100);
            var curve4 = new Bezier(new Point2D[] { new Point2D(0, 0), new Point2D(3, 3), new Point2D(7, 0) }, 100);

            var points1 = curve1.GetCurveCoords();
            var points2 = curve2.GetCurveCoords();

            int i, j;

            var minDist = FindMinDistance(curve1.GetCurvePoints(), curve2.GetCurvePoints(), out i, out j);
            //var crossPoint = FindCrossPoint(curve1, curve2, 1e-3);

            //var p1 = curve1.ImproveAccuracy((double)(i - 2) / curve1.accuracy, (double)(i + 2) / curve1.accuracy, 100);
            //var p2 = curve2.ImproveAccuracy((double)(j - 2) / curve1.accuracy, (double)(j + 2) / curve1.accuracy, 100);

            //double[][] coords1 = new double[2][];
            //double[][] coords2 = new double[2][];
            //coords1[0] = new double[p1.Count];
            //coords1[1] = new double[p1.Count];
            //coords2[0] = new double[p2.Count];
            //coords2[1] = new double[p2.Count];
            //for (int k = 0; k < p1.Count; k++)
            //{
            //    coords1[0][k] = p1[k].e1;
            //    coords1[1][k] = p1[k].e2;
            //    coords2[0][k] = p2[k].e1;
            //    coords2[1][k] = p2[k].e2;
            //}

            DrawerCollection d = new DrawerCollection();

            DRWDrawer drawer = d.DRWCreator(pictureBox1, true, true, "");

            drawer.axis_countX = 11;
            drawer.axis_countY = 11;
            drawer.AxisConfig(new int[] { 11, 11 });

            drawer.Resize(new double[] { 0, 10, 0, 10 });
            //drawer.DrawGraph(coords1[0], coords1[1], false);
            //drawer.AddGraph(coords2[0], coords2[1]);
            drawer.DrawGraph(curve1.GetCurveCoords()[0], curve1.GetCurveCoords()[1], false);
            drawer.AddGraph(curve2.GetCurveCoords()[0], curve2.GetCurveCoords()[1]);
        }
    }
}
