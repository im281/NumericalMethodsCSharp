using System;
using XuMath;

namespace ODETest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestEuler();
            //TestRungeKutta2();
            //TestRungeKutta4();
            //TestRungeKuttaFehlberg();  
            //TestMultiRungeKutta4();
            TestMultiRungeKuttaFehlberg();
            Console.ReadLine();
        }

        static void TestEuler()
        {
            double h = 0.001;
            double x0 = 0;
            double y0 = 1.0;
            Console.WriteLine("\n Results from the Euler's method with h = {0}\n", h);
            double result = y0;
            for (int i = 0; i < 11; i++)
            {
                double x = 0.1 * i;
                result = ODE.Euler(f, x0, result, h, x);
                double exact = 2.0 * Math.Exp(x) - x - 1.0;
                Console.WriteLine(" x = {0:n1}, y = {1:e12}, exact = {2:e12}", x, result, exact);
                x0 = x;
            }
        }

        static void TestRungeKutta2()
        {
            double h = 0.001;
            double x0 = 0;
            double y0 = 1.0;
            Console.WriteLine("\n Results from the second-order Runge-Kutta method with h = {0}\n", h);
            double result = y0;
            for (int i = 0; i < 11; i++)
            {
                double x = 0.1 * i;
                result = ODE.RungeKutta2(f, x0, result, h, x);
                double exact = 2.0 * Math.Exp(x) - x - 1.0;
                Console.WriteLine(" x = {0:n1}, y = {1:e12}, exact = {2:e12}", x, result, exact);
                x0 = x;
            }
        }

        static void TestRungeKutta4()
        {
            double h = 0.001;
            double x0 = 0;
            double y0 = 1.0;
            Console.WriteLine("\n Results from the fourth-order Runge-Kutta method with h = {0}\n", h);
            double result = y0;
            for (int i = 0; i < 11; i++)
            {
                double x = 0.1 * i;
                result = ODE.RungeKutta4(f, x0, result, h, x);
                double exact = 2.0 * Math.Exp(x) - x - 1.0;
                Console.WriteLine(" x = {0:n1}, y = {1:e12}, exact = {2:e12}", x, result, exact);
                x0 = x;
            }
        }

        static void TestRungeKuttaFehlberg()
        {
            double h = 0.2;
            double x0 = 0;
            double y0 = 1.0;
            Console.WriteLine("\n Results from the fourth-order Runge-Kutta-Fehlberg method with h = {0}\n", h);
            double result = y0;
            for (int i = 0; i < 11; i++)
            {
                double x = 0.1 * i;
                result = ODE.RungeKuttaFehlberg(f, x0, result, x, h, 1e-8);
                double exact = 2.0 * Math.Exp(x) - x - 1.0;
                Console.WriteLine(" x = {0:n1}, y = {1:e12}, exact = {2:e12}", x, result, exact);
                x0 = x;
            }
        }

        static double m1 = 0.2;
        static double m2 = 0.2;
        static double k1 = 10.0;
        static double k2 = 1.0;
        static double k3 = 10.0;
        static double b1 = 0.01;
        static double b2 = 0.01;
        static double b3 = 0.01;
        static double x10 = 1.0;
        static double x20 = 0.0;
        static double v10 = 0.0;
        static double v20 = 0.0;

        static void TestMultiRungeKutta4()
        {
            double dt = 0.02;
            double t0 = 0.0;
            VectorR x0 = new VectorR(new double[] { x10, x20, v10, v20 });

            Console.WriteLine("\n Results for a coupled spring system with h = {0}:\n", dt);
            Console.WriteLine(" t       x1         x2         v1         v2");

            VectorR x = x0;
            for (int i = 0; i < 21; i++)
            {
                double t = 0.1 * i;
                x = ODE.MultiRungeKutta4(f1, t0, x, dt, t);
                Console.WriteLine(" {0:n2}   {1,8:n5}   {2,8:n5}   {3,8:n5}   {4,8:n5}", t, x[0], x[1], x[2], x[3]);
                t0 = t;
            }
        }

        static void TestMultiRungeKuttaFehlberg()
        {
            double dt = 0.5;
            double t0 = 0.0;
            VectorR x0 = new VectorR(new double[] { x10, x20, v10, v20 });

            Console.WriteLine("\n Results for a coupled spring system with h = {0}:\n", dt);
            Console.WriteLine(" t       x1         x2         v1         v2");

            VectorR x = x0;
            for (int i = 0; i < 21; i++)
            {
                double t = 0.1 * i;
                x = ODE.MultiRungeKuttaFehlberg(f1, t0, x, t, dt, 1e-5);
                Console.WriteLine(" {0:n2}   {1,8:n5}   {2,8:n5}   {3,8:n5}   {4,8:n5}", t, x[0], x[1], x[2], x[3]);
                t0 = t;
            }
        }

        static double f(double x, double y)
        {
            return x + y;
        }

        static VectorR f1(double t, VectorR x)
        {
            VectorR result = new VectorR(4);
            result[0] = x[2];
            result[1] = x[3];
            result[2] = -(k1 + k2) * x[0] / m1 + k2 * x[1] / m1 - (b1 + b2) * x[2] / m1 + b2 * x[3] / m1;
            result[3] = -(k2 + k3) * x[1] / m2 + k2 * x[0] / m2 - (b2 + b3) * x[3] / m2 + b2 * x[2] / m2;
            return result;
        }
    }
}
