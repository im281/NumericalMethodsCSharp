using System;
using XuMath;

namespace OptimizationTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestBisection();
            //TestGoldenSearch();
            //TestNewton();
            //TestBrent();
            //TestMultiNewton();
            //TestSimplex();
            //TestPeaks();
            //TestPeaks1();
            TestDifferentialEvolution();
            
            Console.ReadLine();
        }

        static void TestBisection()
        {
            double result = Optimization.Bisection(f, 0.0, 1.0, 1.0e-5);
            Console.WriteLine("x = " + result.ToString() + ", f(x) = " + f(result).ToString());
        }

        static void TestGoldenSearch()
        {
            double result = Optimization.GoldenSearch(f, 0.0, 1.0, 1.0e-5);
            Console.WriteLine("x = " + result.ToString() + ", f(x) = " + f(result).ToString());
        }

        static void TestNewton()
        {
            double result = Optimization.Newton(f, 0.0, 1.0e-5);
            Console.WriteLine("x = " + result.ToString() + ", f(x) = " + f(result).ToString());
        }

        static void TestBrent()
        {
            double result = Optimization.Brent(f, 0.0, 1.0, 1.0e-5);
            Console.WriteLine("x = " + result.ToString() + ", f(x) = " + f(result).ToString());
        }

        static void TestMultiNewton()
        {
            double[] xarray = new double[] { 0, 0 };
            VectorR result = Optimization.multiNewton(f1, xarray, 1.0e-5);
            Console.WriteLine("x = " + result.ToString());
            Console.WriteLine("f1(x) = " + f1(result).ToString());
        }

        static void TestSimplex()
        {
            MatrixR x = new MatrixR(3, 2);

            x[0, 0] = 0;
            x[0, 1] = 0;
            x[1, 0] = 1;
            x[1, 1] = 0;
            x[2, 0] = 0;
            x[2, 1] = 1;

            VectorR result = Optimization.Simplex(f2, x, 1000);
            Console.WriteLine("x = " + result.ToString());
            Console.WriteLine("f(x) = " + f2(result));
        }

        static void TestAnneal()
        {
            double[] xarray = new double[] { 0, 1};
            VectorR result = Optimization.Anneal(f3, xarray, 1e-15, 200);
            Console.WriteLine(result.ToString());
            Console.WriteLine(" f3 = " + f3(result).ToString());
        }

        static void TestPeaks()
        {
            double[] xarray = new double[] { 0, -2 };
            VectorR result = Optimization.Anneal(Peaks, xarray, 1e-15, 200);
            Console.WriteLine(result.ToString());
            Console.WriteLine(" Peaks = " + Peaks(result).ToString());
        }

        static void TestPeaks1()
        {
            double[] xarray = new double[] { 0, 2 };
            VectorR result = Optimization.Anneal(Peaks1, xarray, 1e-15, 200);
            Console.WriteLine(result.ToString());
            Console.WriteLine(" Peaks = " + Peaks1(result).ToString());
        }

        static void TestDifferentialEvolution()
        {
            VectorR bestmember;
            double bestvalue;
            Optimization.MaxIterations = 201;
            Optimization.MinCost = -50;
            Optimization.Refresh = 50;
            Optimization.DifferentialEvolution(Peaks, out bestmember, out bestvalue);
            Console.WriteLine("\n Minimum = {0}\n Location = {1}\n\n", bestvalue, bestmember); 
            Optimization.DifferentialEvolution(Peaks1, out bestmember, out bestvalue);
            Console.WriteLine("\n Maximum = {0}\n Location = {1}", -bestvalue, bestmember); 
        }

        static double f(double x)
        {
            return 1.6 * x * x * x + 3 * x * x - 2 * x;
        }

        static double f1(VectorR x)
        {
            return 1.5 * (x[0] - 0.5) * (x[0] - 0.5) + 3.4 * (x[1] + 1.2) * (x[1] + 1.2) + 2.5;
        }

        static double f2(VectorR x)
        {
            return x[0] * x[0] - 4 * x[0] + x[1] * x[1] - x[1] - x[0] * x[1];
        }

        static double f3(VectorR x)
        {
            return (4.0 - 2.1 * x[0] * x[0] + Math.Pow(Math.Abs(x[0]), 4.0 / 3.0)) * x[0] * x[0] + x[0] * x[1] + 
                    4.0 * (x[1] * x[1] - 1) * x[1] * x[1];
        }

        static double Peaks(VectorR x)
        {
            double z = 3 * (1 - x[0]) * (1 - x[0]) * Math.Exp(-x[0] * x[0] - (x[1] + 1) * (x[1] + 1))
                   - 10 * (x[0] / 5 - Math.Pow(x[0], 3) - Math.Pow(x[1], 5)) * Math.Exp(-x[0] * x[0] - x[1] * x[1])
                   - 1 / 3 * Math.Exp(-(x[0] + 1) * (x[0] + 1) - x[1] * x[1]);
            return z;
        }

        static double Peaks1(VectorR x)
        {
            return -Peaks(x);
        }
    }
}
