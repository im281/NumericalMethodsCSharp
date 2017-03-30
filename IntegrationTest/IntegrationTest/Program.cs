using System;
using XuMath;

namespace IntegrationTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestTrapezoidal();
            //TestSimpson();
            //TestRomberg();
            //TestGaussLegendre();
            //TestGaussLaguerre();
            //TestGaussHermite();
            TestChebyshev();

            Console.ReadLine();
        }

        static void TestTrapezoidal()
        {

            int n = 101;
            double result;

            result = f1(1) - f1(0);
            Console.WriteLine("\n Analytic result = " + result.ToString());

            result = Integration.Trapezoidal(f, 0, 1, n);
            Console.WriteLine(" Result for function = " + result.ToString());

            double[] ya = new double[n];
            double h = 1.0 / (n - 1);
            for (int i = 0; i < n; i++)
            {
                double x = i * h;
                ya[i] = f(x);
            }
            result = Integration.Trapezoidal(ya, h);
            Console.WriteLine(" Result for data array = " + result.ToString());
        }

        static void TestSimpson()
        {
            int n = 101;
            double result;

            result = f1(1) - f1(0);
            Console.WriteLine("\n Analytic result = " + result.ToString());

            result = Integration.Simpson(f, 0, 1, n);
            Console.WriteLine(" Result for function = " + result.ToString());

            double[] ya = new double[n];
            double h = 1.0 / (n - 1);
            for (int i = 0; i < n; i++)
            {
                double x = i * h;
                ya[i] = f(x);
            }
            result = Integration.Simpson(ya, h);
            Console.WriteLine(" Result for data array = " + result.ToString());
        }

        static void TestRomberg()
        {
            double result;

            result = f1(1) - f1(0);
            Console.WriteLine("\n Analytic result = " + result.ToString());

            result = Integration.Romberg(f, 0, 1, 15, 1e-9);
            Console.WriteLine(" Result from Romberg method = " + result.ToString());
        }

        static void TestGaussLegendre()
        {
            Console.WriteLine("\n Result from Gauss-Legendre method:\n");
            double result;
            for (int n = 1; n < 9; n++)
            {
                result = Integration.GaussLegendre(f2, 1, 2, n);
                Console.WriteLine(" n = {0}, result = {1}", n, result);
            }
        }

        static void TestGaussLaguerre()
        {
            Console.WriteLine("\n Result from Gauss-Laguerre method:\n");
            double result;
            for (int n = 1; n < 9; n++)
            {
                result = Integration.GaussLaguerre(f3, n);
                Console.WriteLine(" n = {0}, result = {1}", n, result);
            }
        }

        static void TestGaussHermite()
        {
            Console.WriteLine("\n Result from Gauss-Hermit method:\n");
            double exact = Math.Sqrt(Math.PI) / 2;
            Console.WriteLine(" Exact result = {0} \n", exact);
            double result;
            for (int n = 1; n < 9; n++)
            {
                result = Integration.GaussHermite(f4, n);
                Console.WriteLine(" n = {0}, result = {1}", n, result);
            }
        }

        static void TestChebyshev()
        {
            Console.WriteLine("\n Result from Gauss-Chebyshev method:\n");
            double result;
            for (int n = 1; n < 9; n++)
            {
                result = Integration.GaussChebyshev(f5, n);
                Console.WriteLine(" n = {0}, result = {1}", n, result);
            }
        }

        static double f(double x)
        {
            return Math.Exp(x) - 3 * x * x;
        }

        static double f1(double x)
        {
            return Math.Exp(x) - x * x * x;
        }

        static double f2(double x)
        {
            return Math.Exp(-0.5 * x * x) / Math.Sqrt(2.0 * Math.PI);
        }

        static double f3(double x)
        {
            return Math.Sin(x);
        }

        static double f4(double x)
        {
            return x * x;
        }

        static double f5(double x)
        {
            return (1 - x * x) * (1 - x * x);
        }
    }
}
