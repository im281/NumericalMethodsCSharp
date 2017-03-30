using System;
using XuMath;

namespace InterpolationTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestLinear();
            //TestLagrangian();
            //TestBarycentric();
            //TestNewtonDividedDifference();
            TestSpline();
            //TestBilinear();
            Console.ReadLine();
        }

        static void TestLinear()
        {
            double[] xarray = new double[] {0, 2, 4, 6, 8 };
            double[] yarray = new double[] { 0, 4, 16, 36, 64 };
            double[] x = new double[] { 1, 3, 5, 7 };
            double[] y = Interpolation.Linear(xarray,yarray,x);
            VectorR vx = new VectorR(x);
            VectorR vy = new VectorR(y);
            Console.WriteLine(" x = " + vx.ToString());
            Console.WriteLine(" y=" + vy.ToString());
        }

        static void TestLagrangian()
        {
            double[] xarray = new double[5] { 1, 2, 3, 4, 5 };
            double[] yarray = new double[5] { 1, 4, 9, 16, 25 };
            double[] x = new double[3] { 2.5, 3.5, 1.5 };

            double[] y = Interpolation.Lagrangian(xarray, yarray, x);
            VectorR vx = new VectorR(x);
            VectorR vy = new VectorR(y);
            Console.WriteLine(" x = " + vx.ToString());
            Console.WriteLine(" y=" + vy.ToString());
        }

        static void TestBarycentric()
        {
            double[] xarray = new double[] { 0, 2, 4, 6, 8 };
            double[] yarray = new double[] { 0, 4, 16, 36, 64 };
            double[] x = new double[] { 1, 3, 5, 7 };
            double[] y = Interpolation.Barycentric(xarray, yarray, x);
            VectorR vx = new VectorR(x);
            VectorR vy = new VectorR(y);
            Console.WriteLine(" x = " + vx.ToString());
            Console.WriteLine(" y=" + vy.ToString());
        }

        static void TestNewtonDividedDifference()
        {
            double[] xarray = new double[] {1950, 1960, 1970, 1980, 1990 };
            double[] yarray = new double[] {150.697, 179.323, 203.212, 226.505, 249.633 };
            double[] x = new double[] { 1955, 1965, 1975, 1985 };
            double[] y = Interpolation.NewtonDividedDifference(xarray, yarray, x);
            VectorR vx = new VectorR(x);
            VectorR vy = new VectorR(y);
            Console.WriteLine(" x = " + vx.ToString());
            Console.WriteLine(" y=" + vy.ToString());
        }

        static void TestSpline()
        {
            double[] xarray = new double[] { 0, 2, 4, 6, 8 };
            double[] yarray = new double[] { 0, 4, 16, 36, 64 };
            double[] x = new double[] { 3, 5, 7 };
            double[] y = Interpolation.Spline(xarray, yarray, x);
            VectorR vx = new VectorR(x);
            VectorR vy = new VectorR(y);
            Console.WriteLine(" x = " + vx.ToString());
            Console.WriteLine(" y=" + vy.ToString());
        }

        static void TestBilinear()
        {
            double[] xarray = new double[] { 0, 1 };
            double[] yarray = new double[] { 0, 1 };
            double[,] zarray = new double[,] { { 0, 10 }, { 10, 5 } };
            double[] x = new double[9];
            double[] y = new double[9];
            MatrixR z = new MatrixR(9, 9);
            for (int i = 0; i < 9; i++)
            {
                x[i] = (i + 1.0) / 10.0;
                y[i] = (i + 1.0) / 10.0;
            }
            VectorR vx = new VectorR(x);
            VectorR vy = new VectorR(y);

            for (int i = 0; i < 9; i++)
            {
                for (int j = 0; j < 9; j++)
                {
                    z[i, j] = Interpolation.Bilinear(xarray, yarray, zarray, x[i], y[j]);
                }
            }
            Console.WriteLine("x = " + vx.ToString());
            Console.WriteLine("y = " + vy.ToString());
            Console.WriteLine("\nResults z = \n" + z.ToString());
        }
    }
}
