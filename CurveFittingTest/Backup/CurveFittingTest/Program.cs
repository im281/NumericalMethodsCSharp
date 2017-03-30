using System;
using XuMath;

namespace CurveFittingTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestStraightLineFit();
            //TestLinearRegression();
            //TestPolynomialFit();
            //TestWeightedLinearRegression();
            //TestSimpleMovingAverage();
            //TestWeightedMovingAverage();
            TestExponentialMovingAverage();
            Console.ReadLine();
        }

        private static void TestStraightLineFit()
        {
            double[] xarray = new double[] { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
            double[] yarray = new double[] { 1.9, 2.7, 3.3, 4.4, 5.5, 6.5 };
            double[] results = CurveFitting.StraightLineFit(xarray, yarray);
            VectorR v = new VectorR(results);
            Console.WriteLine(v.ToString());
        }

        private static void TestLinearRegression()
        {
            double[] xarray = new double[] { 0, 1, 2, 3, 4, 5 };
            double[] yarray = new double[] { 2, 1, 4, 4, 3, 2 };

            // First order polynomial (m = 1):
            CurveFitting.ModelFunction[] f = new CurveFitting.ModelFunction[] { f0, f1};
            double sigma = 0.0;
            VectorR results = CurveFitting.LinearRegression(xarray, yarray,f, out sigma);
            Console.WriteLine("Order of polynomial m = 1" + ", Standard deviation = " + sigma.ToString());
            Console.WriteLine("Ceofficients = " + results.ToString() +"\n");

            //Second order polynomial (m = 2):
            f = new CurveFitting.ModelFunction[] { f0, f1, f2 };
            results = CurveFitting.LinearRegression(xarray, yarray, f, out sigma);
            Console.WriteLine("Order of polynomial m = 2" + ", Standard deviation = " + sigma.ToString());
            Console.WriteLine("Ceofficients = " + results.ToString() + "\n");

            //Third order polynomial (m = 3):
            f = new CurveFitting.ModelFunction[] { f0, f1, f2, f3 };
            results = CurveFitting.LinearRegression(xarray, yarray, f, out sigma);
            Console.WriteLine("Order of polynomial m = 3" + ", Standard deviation = " + sigma.ToString());
            Console.WriteLine("Ceofficients = " + results.ToString() + "\n");

        }

        private static double f0(double x)
        {
            return 1.0;
        }

        private static double f1(double x)
        {
            return x;
        }

        private static double f2(double x)
        {
            return x * x;
        }

        private static double f3(double x)
        {
            return x * x * x;
        }

        private static void TestPolynomialFit()
        {
            double[] xarray = new double[] { 0, 1, 2, 3, 4, 5 };
            double[] yarray = new double[] {2, 1, 4, 4, 3, 2 };
            for (int m = 1; m < 4; m++)
            {
                double sigma = 0.0;
                VectorR results = CurveFitting.PolynomialFit(xarray, yarray, m, out sigma);
                Console.WriteLine("\nOrder of polynomial m = " + m.ToString() + ", Standard deviation = " + sigma.ToString());
                Console.WriteLine("Ceofficients = " + results.ToString());
            }
        }

        private static void TestWeightedLinearRegression()
        {
            double[] x = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            double[] y = new double[] { 1.9398, 2.9836, 5.9890, 10.2000, 20.7414, 23.2320, 69.5855, 82.5836, 98.1779, 339.3256 };
            double[] ylog = new double[] { 0.6626, 1.0931, 1.7899, 2.3224, 3.0321, 3.1455, 4.2426, 4.4138, 4.5868, 5.8270 };
            double[] results = CurveFitting.WeightedLinearRegression(x, ylog, y);
            VectorR v = new VectorR(results);
            Console.WriteLine(v.ToString());
        }

        private static void TestSimpleMovingAverage()
        {
            double[] data = new double[] {45.375, 45.500, 45.000, 43.625, 43.375, 43.125, 43.125, 44.250,
                                          43.500, 44.375, 45.875, 46.750, 47.625, 48.000, 49.125, 48.750,
                                          46.125, 46.750, 46.625, 46.000};
            VectorR sma = CurveFitting.SimpleMovingAverage(data, 5);
            Console.WriteLine(sma.ToString());
        }

        private static void TestWeightedMovingAverage()
        {
            double[] data = new double[] {45.375, 45.500, 45.000, 43.625, 43.375, 43.125, 43.125, 44.250,
                                          43.500, 44.375, 45.875, 46.750, 47.625, 48.000, 49.125, 48.750,
                                          46.125, 46.750, 46.625, 46.000};
            VectorR sma = CurveFitting.WeightedMovingAverage(data, 5);
            Console.WriteLine(sma.ToString());
        }

        private static void TestExponentialMovingAverage()
        {
            double[] data = new double[] {45.375, 45.500, 45.000, 43.625, 43.375, 43.125, 43.125, 44.250,
                                          43.500, 44.375, 45.875, 46.750, 47.625, 48.000, 49.125, 48.750,
                                          46.125, 46.750, 46.625, 46.000};
            VectorR sma = CurveFitting.ExponentialMovingAverage(data, 5);
            Console.WriteLine(sma.ToString());
        }
    }
}
