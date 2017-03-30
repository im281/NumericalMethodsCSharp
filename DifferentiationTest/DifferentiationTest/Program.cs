using System;
using XuMath;

namespace DifferentiationTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestForwardMethod();
            //TestBackwardMethod();
            //TestCentralMethod();
            //TestExtendedCentralMethod();
            TestRichardson();
            //TestInterpolation();

            Console.ReadLine();
        }

        static void TestForwardMethod()
        {
            double h = 0.1;
            double[] yarray=new double[10];
            // create yarray:
            for (int i = 0; i < 10; i++)
                yarray[i] = f(i * h);

            // Calculate derivatives for function:
            Console.WriteLine("\n Derivatives for f(x) = sin(x):\n");
            double dy = Differentiation.Forward1(f, 0.3, h);
            Console.WriteLine(" f'(x) = " + dy.ToString());
            dy = Differentiation.Forward2(f, 0.3, h);
            Console.WriteLine(" f''(x) = " + dy.ToString());
            dy = Differentiation.Forward3(f, 0.3, h);
            Console.WriteLine(" f'''(x) = " + dy.ToString());
            dy = Differentiation.Forward4(f, 0.3, h);
            Console.WriteLine(" f''''(x) = " + dy.ToString());

            // Calculate derivatives for array values:
            Console.WriteLine("\n Derivatives for array values:");
            Console.WriteLine("\n yarray ="); 
            foreach (double y1 in yarray)
                Console.WriteLine(" " + y1.ToString());
            dy = Differentiation.Forward1(yarray, 3, h);
            Console.WriteLine("\n y' = " + dy.ToString());
            dy = Differentiation.Forward2(yarray, 3, h);
            Console.WriteLine(" y'' = " + dy.ToString());
            dy = Differentiation.Forward3(yarray, 3, h);
            Console.WriteLine(" y''' = " + dy.ToString());
            dy = Differentiation.Forward4(yarray, 3, h);
            Console.WriteLine(" y'''' = " + dy.ToString());
        }

        static void TestBackwardMethod()
        {
            double h = 0.1;
            double[] yarray = new double[10];
            // create yarray:
            for (int i = 0; i < 10; i++)
                yarray[i] = f(i * h);

            // Calculate derivatives for function at 0.7:
            Console.WriteLine("\n Derivatives for f(x) = sin(x):\n");
            double dy = Differentiation.Backward1(f, 0.7, h);
            Console.WriteLine(" f'(x) = " + dy.ToString());
            dy = Differentiation.Backward2(f, 0.7, h);
            Console.WriteLine(" f''(x) = " + dy.ToString());
            dy = Differentiation.Backward3(f, 0.7, h);
            Console.WriteLine(" f'''(x) = " + dy.ToString());
            dy = Differentiation.Backward4(f, 0.7, h);
            Console.WriteLine(" f''''(x) = " + dy.ToString());

            // Calculate derivatives for array values at yindex = 7:
            Console.WriteLine("\n Derivatives for array values:");
            Console.WriteLine("\n yarray =");
            foreach (double y1 in yarray)
                Console.WriteLine(" " + y1.ToString());
            dy = Differentiation.Backward1(yarray, 7, h);
            Console.WriteLine("\n y' = " + dy.ToString());
            dy = Differentiation.Backward2(yarray, 7, h);
            Console.WriteLine(" y'' = " + dy.ToString());
            dy = Differentiation.Backward3(yarray, 7, h);
            Console.WriteLine(" y''' = " + dy.ToString());
            dy = Differentiation.Backward4(yarray, 7, h);
            Console.WriteLine(" y'''' = " + dy.ToString());
        }

        static void TestCentralMethod()
        {
            double h = 0.1;
            double[] yarray = new double[10];
            // create yarray:
            for (int i = 0; i < 10; i++)
                yarray[i] = f(i * h);

            // Calculate derivatives for function at 0.5:
            Console.WriteLine("\n Derivatives for f(x) = sin(x):\n");
            double dy = Differentiation.Central1(f, 0.5, h);
            Console.WriteLine(" f'(x) = " + dy.ToString());
            dy = Differentiation.Central2(f, 0.5, h);
            Console.WriteLine(" f''(x) = " + dy.ToString());
            dy = Differentiation.Central3(f, 0.5, h);
            Console.WriteLine(" f'''(x) = " + dy.ToString());
            dy = Differentiation.Central4(f, 0.5, h);
            Console.WriteLine(" f''''(x) = " + dy.ToString());

            // Calculate derivatives for array values at yindex = 5:
            Console.WriteLine("\n Derivatives for array values:");
            Console.WriteLine("\n yarray =");
            foreach (double y1 in yarray)
                Console.WriteLine(" " + y1.ToString());
            dy = Differentiation.Central1(yarray, 5, h);
            Console.WriteLine("\n y' = " + dy.ToString());
            dy = Differentiation.Central2(yarray, 5, h);
            Console.WriteLine(" y'' = " + dy.ToString());
            dy = Differentiation.Central3(yarray, 5, h);
            Console.WriteLine(" y''' = " + dy.ToString());
            dy = Differentiation.Central4(yarray, 5, h);
            Console.WriteLine(" y'''' = " + dy.ToString());
        }

        static void TestExtendedCentralMethod()
        {
            double h = 0.1;
            double[] yarray = new double[10];
            // create yarray:
            for (int i = 0; i < 10; i++)
                yarray[i] = f(i * h);

            // Calculate derivatives for function at 0.5:
            Console.WriteLine("\n Derivatives for f(x) = sin(x):\n");
            double dy = Differentiation.ExtendedCentral1(f, 0.5, h);
            Console.WriteLine(" f'(x) = " + dy.ToString());
            dy = Differentiation.ExtendedCentral2(f, 0.5, h);
            Console.WriteLine(" f''(x) = " + dy.ToString());
            dy = Differentiation.ExtendedCentral3(f, 0.5, h);
            Console.WriteLine(" f'''(x) = " + dy.ToString());
            dy = Differentiation.ExtendedCentral4(f, 0.5, h);
            Console.WriteLine(" f''''(x) = " + dy.ToString());

            // Calculate derivatives for array values at yindex = 5:
            Console.WriteLine("\n Derivatives for array values:");
            dy = Differentiation.ExtendedCentral1(yarray, 5, h);
            Console.WriteLine("\n y' = " + dy.ToString());
            dy = Differentiation.ExtendedCentral2(yarray, 5, h);
            Console.WriteLine(" y'' = " + dy.ToString());
            dy = Differentiation.ExtendedCentral3(yarray, 5, h);
            Console.WriteLine(" y''' = " + dy.ToString());
            dy = Differentiation.ExtendedCentral4(yarray, 5, h);
            Console.WriteLine(" y'''' = " + dy.ToString());

            // Analytic results:
            Console.WriteLine("\n Analytic results:");
            dy = Math.Cos(0.5);
            Console.WriteLine("\n y' = " + dy.ToString());
            dy = -Math.Sin(0.5);
            Console.WriteLine(" y'' = " + dy.ToString());
            dy = -Math.Cos(0.5);
            Console.WriteLine(" y''' = " + dy.ToString());
            dy = Math.Sin(0.5);
            Console.WriteLine(" y'''' = " + dy.ToString());
        }

        static void TestRichardson()
        {
            double h = 0.1;
            double x = 0.5;
            Console.WriteLine("\n Comparison results:\n");
            double exact1 = 4 * x * x * x + 1 - Math.Exp(-x);
            double c1 = Differentiation.Central1(f1, x, h);
            double r1 = Differentiation.Richardson1(f1, x, h, "Central");
            Console.WriteLine(" exact1 = {0,11:n8}, c1 = {1,11:n8}, r1 = {2,11:n8}", exact1, c1, r1);
            double exact2 = 12 * x * x + Math.Exp(-x);
            double c2 = Differentiation.Central2(f1, x, h);
            double r2 = Differentiation.Richardson2(f1, x, h, "Central");
            Console.WriteLine(" exact2 = {0,11:n8}, c2 = {1,11:n8}, r2 = {2,11:n8}", exact2, c2, r2);
            double exact3 = 24 * x - Math.Exp(-x);
            double c3 = Differentiation.Central3(f1, x, h);
            double r3 = Differentiation.Richardson3(f1, x, h, "Central");
            Console.WriteLine(" exact3 = {0,11:n8}, c3 = {1,11:n8}, r3 = {2,11:n8}", exact3, c3, r3);
            double exact4 = 24 + Math.Exp(-x);
            double c4 = Differentiation.Central4(f1, x, h);
            double r4 = Differentiation.Richardson4(f1, x, h, "Central");
            Console.WriteLine(" exact4 = {0,11:n8}, c4 = {1,11:n8}, r4 = {2,11:n8}", exact4, c4, r4); 
        }

        static void TestInterpolation()
        {
            double x = 1.8;
            double[] xa = new double[] { 1.0, 1.1, 1.3, 1.6, 1.7, 2.0};
            double[] ya = new double[] { 1.2772, 1.1414, 0.7880, 0.1435, -0.0729, -0.6215 };
            double y0 = Differentiation.Interpolation0(xa, ya, x);
            double y1 = Differentiation.Interpolation1(xa, ya, x);
            double y2 = Differentiation.Interpolation2(xa, ya, x);
            Console.WriteLine(" y = {0:n6}, y' = {1:n6}, y'' = {2:n6}", y0, y1, y2);
        }

        static double f(double x)
        {
            return Math.Sin(x);
        }

        static double f1(double x)
        {
            return x * x * x * x + x + Math.Exp(-x);
        }
    }
}
