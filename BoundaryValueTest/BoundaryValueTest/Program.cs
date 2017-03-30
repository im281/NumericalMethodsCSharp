using System;
using XuMath;

namespace BoundaryValueTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestShooting2();
            //TestFiniteDifferenceLinear2();
            //TestfiniteDifferenceNonlinear2();
            TestFiniteDifferenceLinear4();
            Console.ReadLine();
        }

        static void TestShooting2()
        {
            BoundaryValue bv = new BoundaryValue();
            bv.xa = 0.0;
            bv.xb = 1.0;
            bv.ya = 0.0;
            bv.yb = 2.0;
            bv.u1 = 0.0;
            bv.u2 = 2.0;
            bv.StepSize = 0.05;
            bv.F1 = f1;

            Console.WriteLine("\n Results from the shooting method:\n");
            Console.WriteLine(" x        y           y'         Exact y     Exact y'");
            for (int i = 0; i < 11; i++)
            {
                bv.xOut = 0.1 * i;
                VectorR y = bv.Shooting2();
                double yexact = bv.xOut * bv.xOut * bv.xOut + bv.xOut;
                double y1exact = 3 * bv.xOut * bv.xOut + 1.0;
                Console.WriteLine(" {0:n3}, {1,10:n6}, {2,10:n6} {3,10:n6}, {4,10:n6}", bv.xOut, y[0], y[1], yexact, y1exact);
            }
        }

        static void TestFiniteDifferenceLinear2()
        {
            BoundaryValue bv = new BoundaryValue();
            /*bv.xa = 0.0;
            bv.xb = 1.0;
            bv.ya = 0.0;
            bv.yb = 0.0;
            bv.n = 8;*/
            bv.BoundaryFlag = new int[] { 0, 1 };
            bv.xa = 0.0;
            bv.xb = 0.5 * Math.PI;
            bv.ya = 0.0;
            bv.vb = 0;
            bv.n = 10;
            double[] x;
            VectorR y = bv.FiniteDifferenceLinear2(f3, out x);

            Console.WriteLine("\n  x           y            Exact y");
            for (int i = 0; i < x.Length; i++)
            {
                //double exact = x[i] - x[i] * x[i] * x[i];
                double exact = 2 * x[i] + Math.Sin(2 * x[i]);
                Console.WriteLine(" {0,8:n5}   {1,10:n6}   {2,10:n6}", x[i], y[i], exact);
            }
        }

        static void TestfiniteDifferenceNonlinear2()
        {
            BoundaryValue bv = new BoundaryValue();
            bv.xa = 0.0;
            bv.xb = 1.0;
            bv.ya = 0.0;
            bv.yb = 1.0;
            bv.n = 10;
            bv.fd = f4;
            double[] x;
            VectorR y = bv.FiniteDifferenceNonlinear2(out x);

            Console.WriteLine("\n  x        y");
            for (int i = 0; i < x.Length; i++)
            {
                Console.WriteLine(" {0,5:n2}   {1,10:n6}", x[i], y[i]);
            }
        }

        static void TestFiniteDifferenceLinear4()
        {
            BoundaryValue bv = new BoundaryValue();
            bv.xa = 0.0;
            bv.xb = 1.0;
            bv.ya = 0.0;
            bv.va = 0.0;
            bv.vb = 0.0;
            bv.gb = -5;
            bv.n = 10;
            bv.BoundaryFlag = new int[] { 0, 1 };
            
            double[] x;
            VectorR y = bv.FiniteDifferenceLinear4(f5, out x);

            Console.WriteLine("\n  x           y");
            for (int i = 0; i < x.Length; i++)
            {
                Console.WriteLine(" {0,8:n5}   {1,10:n6}", x[i], y[i]);
            }
        }

        static VectorR f5(double x)
        {
            VectorR result = new VectorR(2);
            result[0] = Math.Sin(x);
            result[1] = x;
            return result;
        }

        static double f4(double x, double y, double yprime)
        {
            return x * y * yprime + Math.Sin(2 * x * y);
        }

        static VectorR f3(double x)
        {
            VectorR result = new VectorR(4);
            result[0] = 1;
            result[1] = 0;
            result[2] = 4;
            result[3] = -8 * x;
            return result;
        }

        static VectorR f2(double x)
        {
            VectorR result = new VectorR(4);
            result[0] = x;
            result[1] = -2;
            result[2] = 0;
            result[3] = 2;
            return result;
        }

        static VectorR f1(double x, VectorR y)
        {
            VectorR result = new VectorR(2);
            result[0] = y[1];
            result[1] = -2 * x * y[1] + 6 * y[0] + 2 * x;
            return result;
        }

    }
}
