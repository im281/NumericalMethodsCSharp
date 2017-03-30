using System;
using XuMath;

namespace NonlinearSystemTest
{
    class Program
    {
        static void Main(string[] args)
        {
            TestIncrementSearch();
            //TestFixedPoint();
            //TestBisection();
            //TestFalsePosition();
            //TestNewtonRaphson();
            //TestSecant();
            //TestNewtonMultiRoots();
            //TestbirgeVieta();
            //TestNewtonMultiEquations();
            Console.ReadLine();
        }

        static void TestIncrementSearch()
        {
            double h = 0.01;
            int n = 500;
            double x = -3;
            for (int i = 1; i <= 3; i++)
            {
                x = NonlinearSystem.IncrementSearch(F, x, h, n);
                Console.WriteLine("\n Solution " + i.ToString() + " = " + x.ToString());
                Console.WriteLine(" Solution confirmation: f(x) = " + F(x).ToString());
            }
        }

        static void TestFixedPoint()
        {
            double tol = 0.0001;
            int n = 10000;
            double x0 = 1.6;
            double x = NonlinearSystem.FixedPoint(G, x0, tol, n);
            Console.WriteLine("solution from the fixed point method: " + x.ToString());
        }


        static void TestBisection()
        {
            double x = NonlinearSystem.Bisection(F, 1.0, 2.0, 0.0001);
            Console.WriteLine("Solution from the bisection method: " + x.ToString());
        }

        static void TestFalsePosition()
        {
            double x = NonlinearSystem.FalsePosition(F, 1.0, 2.0, 0.0001);
            Console.WriteLine("Solution from the false position method: " + x.ToString());
        }

        static void TestNewtonRaphson()
        {
            double x = NonlinearSystem.NewtonRaphson(FF, FF1, 1.0, 0.0001);
            Console.WriteLine("Solution from the Newton-Raphson method: " + x.ToString());
        }

        static void TestSecant()
        {
            double x = NonlinearSystem.Secant(F, 1.0, 1.5, 0.0001);
            Console.WriteLine("Solution from the secant method: " + x.ToString());
        }

        static void TestNewtonMultiRoots()
        {
            double[] x = NonlinearSystem.NewtonMultiRoots(F, 0.0, 3, 1000, 0.0001);
            Console.WriteLine(x[0].ToString());
            Console.WriteLine(x[1].ToString());
            Console.WriteLine(x[2].ToString());
        }

        static void TestbirgeVieta()
        {
            double[] x = NonlinearSystem.BirgeVieta(new double[4] { 1, -3, 0, 1.0 }, 0.0, 3, 3, 1000, 0.0001);
            Console.WriteLine(x[0].ToString());
            Console.WriteLine(x[1].ToString());
            Console.WriteLine(x[2].ToString());
        }

        static void TestNewtonMultiEquations()
        {
            VectorR x0 = new VectorR(new double[] { 0.0, 0.0 });
            VectorR x = NonlinearSystem.NewtonMultiEquations(FV, x0, 1e-5);
            Console.WriteLine("\n x[0] =  {0,8:n6}, x[1] = {1,8:n6} \n   f1 = {2,8:n6},   f2 = {3,8:n6}", x[0], x[1],FV(x)[0], FV(x)[1]);
        }

        static VectorR FV(VectorR x)
        {
            VectorR result = new VectorR(2);
            result[0] = 2 * x[0] - x[1] - Math.Exp(-2 * x[0]);
            result[1] = -x[0] + 2 * x[1] - Math.Exp(-x[1]);
            return result;
        }

        static double F(double x)
        {
            return x * x * x - 3.0 * x + 1.0;
        }

        static double FF(double x)
        {
            return Math.Cos(x) - x * x * x;
        }

        static double FF1(double x)
        {
            return -Math.Sin(x) - 3 * x * x;
        }

        static double G(double x)
        {
            return Math.Sqrt(2 * x + 3);
        }
    }
}
