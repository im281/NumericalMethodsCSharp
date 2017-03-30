using System;
using System.Collections;

namespace XuMath
{
    public class Integration
    {
        public delegate double Function(double x);
        public const double badResult = double.NaN;

        public static double Trapezoidal(Function f, double a, double b, int n)
        {
            double sum = 0.0;
            double h = (b - a) / (n - 1);
            for (int i = 0; i < n - 1; i++)
            {
                sum += 0.5 * h * (f(a + i * h) + f(a + (i + 1) * h));
            }
            return sum;
        }

        public static double Trapezoidal(double[] yarray, double h)
        {
            int n = yarray.Length;
            double sum = 0.0;
            for (int i = 0; i < n - 1; i++)
            {
                sum += 0.5 * h * (yarray[i] + yarray[i + 1]);
            }
            return sum;
        }

        public static double Simpson(Function f, double a, double b, int n)
        {
            if (n < 3)
                return badResult;
            double sum = 0.0;
            double h = (b - a) / (n - 1);
            if (n % 2 != 0)
            {
                for (int i = 0; i < n - 2; i += 2)
                {
                    sum += h * (f(a + i * h) + 4 * f(a + (i + 1) * h) + f(a + (i + 2) * h)) / 3;
                }
            }
            else
            {
                sum = 3 * h * (f(a) + 3 * f(a + h) + 3 * f(a + 2 * h) + f(a + 3 * h)) / 8;
                for (int i = 3; i < n - 2; i += 2)
                {
                    sum += h * (f(a + i * h) + 4 * f(a + (i + 1) * h) + f(a + (i + 2) * h)) / 3;
                }
            }
            return sum;
        }

        public static double Simpson(double[] yarray, double h)
        {
            int n = yarray.Length;
            if (n < 3 || h == 0)
                return badResult;

            double sum = 0.0;
            if (n % 2 != 0)
            {
                for (int i = 0; i < n - 2; i += 2)
                {
                    sum += h * (yarray[i] + 4 * yarray[i + 1] + yarray[i + 2]) / 3;
                }
            }
            else
            {
                sum = 3 * h * (yarray[0] + 3 * yarray[1] + 3 * yarray[2] + yarray[3]) / 8;
                for (int i = 3; i < n - 2; i += 2)
                {
                    sum += h * (yarray[i] + 4 * yarray[i + 1] + yarray[i + 2]) / 3;
                }
            }
            return sum;
        }

        public static double Romberg(Function f, double a, double b, int maxIterations, double tolerance)
        {
            int n = (int)Math.Pow(2,maxIterations) + 1;
            double[] T = new double[n];
            double[] fn = new double[n - 2];

            double c;
            double h = b - a;
            double fa = 0.5 * (f(a) + f(b));
            T[0] = h * fa;
            double result = T[0];
            int i, j;
            int nIterations = 0;
            int nSteps = 1;
            do
            {
                result = T[nIterations];
                nIterations++;
                h /= 2;
                nSteps *= 2;
                c = T[0];

                j = 0;
                i = nSteps - 1;
                do
                {
                    j++;
                    fn[i - 1] = f(a + i * h);
                    if (i > 1)
                        fn[i - 2] = fn[nSteps / 2 - j - 1];
                    i -= 2;
                }
                while (i >= 1);

                T[0] = fa;
                for (i = 1; i < nSteps; i++)
                {
                    T[0] += fn[i - 1];
                }
                T[0] *= h;
                for (i = 2; i < nIterations + 2; i++)
                {
                    T[i - 1] = (Math.Pow(4, i - 1) * T[i - 2] - c) / (Math.Pow(4, i - 1) - 1);
                }
            }
            while (nIterations < maxIterations && Math.Abs(T[nIterations] - result) > tolerance) ;
            return T[nIterations];
        }

        public static double GaussLegendre(Function f, double a, double b, int n)
        {
            double[] x, w;
            LegendreNodesWeights(n, out x, out w);

            double sum = 0.0;
            for (int i = 0; i < n; i++)
            {
                sum += 0.5 * (b - a) * w[i] * f(0.5 * (a + b) + 0.5 * (b - a) * x[i]);
            }
            return sum;
        }

        public static void LegendreNodesWeights(int n, out double[] x, out double[] w)
        {
            double c, d, p1, p2, p3, dp;

            x = new double[n];
            w = new double[n];

            for (int i = 0; i < (n + 1) / 2; i++)
            {
                c = Math.Cos(Math.PI * (4 * i + 3) / (4 * n + 2));
                do
                {
                    p2 = 0;
                    p3 = 1;
                    for (int j = 0; j < n; j++)
                    {
                        p1 = p2;
                        p2 = p3;
                        p3 = ((2 * j + 1) * c * p2 - j * p1) / (j + 1);
                    }
                    dp = n * (c * p3 - p2) / (c * c - 1);
                    d = c;
                    c -= p3 / dp;
                }
                while (Math.Abs(c - d) > 1e-12);
                x[i] = c;
                x[n - 1 - i] = -c;
                w[i] = 2 * (1 - x[i] * x[i]) / (n + 1) / (n + 1) / SpecialFunctions.Legendre(x[i], n + 1) / 
                       SpecialFunctions.Legendre(x[i], n + 1);
                w[n - 1 - i] = w[i];
            }
        }

        public static double GaussLaguerre(Function f, int n)
        {
            double[] x, w;
            LaguerreNodesWeights(n, out x, out w);

            double sum = 0.0;
            for (int i = 0; i < n; i++)
            {
                sum += w[i] * f(x[i]);
            }
            return sum;
        }

        public static void LaguerreNodesWeights(int n, out double[] x, out double[] w)
        {
            double c = 0.0;
            double d, p1, p2, p3, dp;

            x = new double[n];
            w = new double[n];
            for (int i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    c = 3 / (1 + 2.4 * n);
                }
                else
                {
                    if (i == 1)
                    {
                        c += 15 / (1 + 2.5 * n);
                    }
                    else
                    {
                        c+= (1 + 2.55 * (i - 1)) / (1.9 * (i - 1))  * (c - x[i - 2]);
                    }
                }
                do
                {
                    p2 = 0;
                    p3 = 1;
                    for (int j = 0; j < n; j++)
                    {
                        p1 = p2;
                        p2 = p3;
                        p3 = ((-c + 2 * j + 1) * p2 - j * p1) / (j + 1);
                    }
                    dp = (n * p3 - n * p2) / c;
                    d = c;
                    c = c - p3 / dp;
                }
                while (Math.Abs(c- d) > 1e-12);
                x[i] = c;
                w[i] = x[i] / (n + 1) / (n + 1) / SpecialFunctions.Laguerre(x[i], n + 1) / SpecialFunctions.Laguerre(x[i], n + 1);
            }
        }

        public static double GaussHermite(Function f, int n)
        {
            double[] x, w;
            HermiteNodesWeights(n, out x, out w);

            double sum = 0.0;
            for (int i = 0; i < n; i++)
            {
                sum += w[i] * f(x[i]);
            }
            return sum;
        }

        public static void HermiteNodesWeights(int n, out double[] x, out double[] w)
        {
            double c = 0.0;
            double d, p1, p2, p3, dp;

            x = new double[n];
            w = new double[n];
            for (int i = 0; i < (n + 1) / 2; i++)
            {
                if (i == 0)
                {
                    c = Math.Sqrt(2 * n + 1) - 1.85575 * Math.Pow(2 * n + 1, -((double)(1) / (double)(6)));
                }
                else
                {
                    if (i == 1)
                    {
                        c = c - 1.14 * Math.Pow(n, 0.426) / c;
                    }
                    else
                    {
                        if (i == 2)
                        {
                            c = 1.86 * c - 0.86 * x[0];
                        }
                        else
                        {
                            if (i == 3)
                            {
                                c = 1.91 * c - 0.91 * x[1];
                            }
                            else
                            {
                                c = 2 * c - x[i - 2];
                            }
                        }
                    }
                }
                do
                {
                    p2 = 0;
                    p3 = Math.Pow(Math.PI, -0.25);
                    for (int j = 0; j < n; j++)
                    {
                        p1 = p2;
                        p2 = p3;
                        p3 = p2 * c * Math.Sqrt((double)(2) / ((double)(j + 1))) - p1 * Math.Sqrt((double)(j) / ((double)(j + 1)));
                    }
                    dp = Math.Sqrt(2 * n) * p2;
                    d = c;
                     c -=  p3 / dp;
                }
                while (Math.Abs(c - d) > 1e-12);
                x[i] = c;
                w[i] = Math.Pow(2, n + 1) * SpecialFunctions.Gamma(n + 1) * Math.Sqrt(Math.PI) /
                       SpecialFunctions.Hermite(x[i], n + 1) / SpecialFunctions.Hermite(x[i], n + 1);
                x[n - 1 - i] = -x[i];
                w[n - 1 - i] = w[i];
            }
        }

        public static double GaussChebyshev(Function f, int n)
        {
            double x;
            double w = Math.PI / n;

            double sum = 0.0;
            for (int i = 0; i < n; i++)
            {
                x = Math.Cos(Math.PI * (i + 0.5) / n);
                sum += f(x);
            }
            return sum * w;
        }
    }
}
