using System;
using XuMath;

namespace SpecialFunctionsTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestGamma();
            //TestBeta();
            //TestError();
            //TestSiCi();
            //TestLaguerre();
            //TestHermite();
            //TestChebyshev();
            //TestLegendre();
            TestBessel();
            Console.ReadLine();
        }

        static void TestGamma()
        {
            Console.WriteLine("Gamma(5) = " + SpecialFunctions.Gamma(5).ToString());
            Console.WriteLine("Gamma(1/2) = " + SpecialFunctions.Gamma(0.5).ToString());
        }

        static void TestBeta()
        {
            Console.WriteLine("Beta(2,3) = " + SpecialFunctions.Beta(2, 3).ToString());
        }

        static void TestError()
        {
            for (int i = 0; i < 21; i++)
            {
                double x = (i - 10.0) / 2.0;
                Console.WriteLine("x = {0,5:n2}, erf(x) = {1,20:e12}, erfc(x) = {2:e12}", x,
                    SpecialFunctions.Erf(x), SpecialFunctions.Erfc(x));
            }
        }

        static void TestSiCi()
        {
            for (int i = 1; i < 21; i++)
            {
                double x = i / 2.0;
                Console.WriteLine("x = {0,5:n2}, Si(x) = {1,18:e12}, Ci(x) = {2,20:e12}", x,
                    SpecialFunctions.Si(x), SpecialFunctions.Ci(x));
            }
        }

        static void TestLaguerre()
        {
            for (int i = 0; i < 7; i++)
            {
                double x = 1.0 * i - 1.0;
                Console.WriteLine("x = {0,5:n2}, L20(x) = {1,20:e12}", x, SpecialFunctions.Laguerre(x, 20));
            }
        }

        static void TestHermite()
        {
            for (int i = 0; i < 7; i++)
            {
                double x = 1.0 * i - 1.0;
                Console.WriteLine("x = {0,5:n2}, H10(x) = {1,20:e12}", x, SpecialFunctions.Hermite(x, 10));
            }
        }

        static void TestChebyshev()
        {
            for (int i = 0; i < 11; i++)
            {
                double x = 0.25 * (i - 5.0);
                Console.WriteLine("x = {0,5:n2}, T15(x) = {1,20:e12}, U15(x) = {2,20:e12}",
                    x, SpecialFunctions.ChebyshevT(x, 15), SpecialFunctions.ChebyshevU(x, 15));
            }
        }

        static void TestLegendre()
        {
            for (int i = 0; i < 9; i++)
            {
                double x = 0.25 * (i - 4.0);
                Console.WriteLine("x = {0,5:n2}, P10(x) = {1,20:e12}", x, SpecialFunctions.Legendre(x, 10));
            }
        }

        static void TestBessel()
        {
            for (int i = 1; i < 21; i++)
            {
                double x = 1.0 * i;
                Console.WriteLine("x = {0,2:n0}, J0(x) = {1,18:e10}, Y0(x) = {2,18:e10}",
                    x, SpecialFunctions.BesselJ(x, 0), SpecialFunctions.BesselY(x, 0.0+1.0e-5));
            }
        }

    }
}
