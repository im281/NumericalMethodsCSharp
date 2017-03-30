using System;
using System.Collections;

namespace XuMath
{
    public class CurveFitting
    {
        public static double[] StraightLineFit(double[] xarray, double[] yarray)
        {
            int n = xarray.Length;
            double xm = 0.0;
            double ym = 0.0;
            double b1 = 0.0;
            double b2 = 0.0;
            double a = 0.0;
            double b = 0.0;
            double s = 0.0;
            double sigma = 0.0;

            for (int i = 0; i < n; i++)
            {
                xm += xarray[i] / n;
                ym += yarray[i] / n;
            }

            for (int i = 0; i < n; i++)
            {
                b1 += yarray[i] * (xarray[i] - xm);
                b2 += xarray[i] * (xarray[i] - xm);
            }
            b = b1 / b2;
            a = ym - xm * b;

            for (int i = 0; i < n; i++)
            {
                s += (yarray[i] - a - b * xarray[i]) * (yarray[i] - a - b * xarray[i]);
            }
            sigma = Math.Sqrt(s / (n - 2));

            return new double[] {a, b, sigma };
        }

        public delegate double ModelFunction(double x);
        public static VectorR LinearRegression(double[] xarray, double[] yarray, ModelFunction[] f, out double sigma)
        {
            int m = f.Length;
            MatrixR A = new MatrixR(m, m);
            VectorR b = new VectorR(m);
            int n = xarray.Length;

            for (int k = 0; k < m; k++)
            {
                b[k] = 0.0;
                for (int i = 0; i < n; i++)
                {
                    b[k] += f[k](xarray[i]) * yarray[i];
                }
            }

            for (int j = 0; j < m; j++)
            {
                for (int k = 0; k < m; k++)
                {
                    A[j, k] = 0.0;
                    for (int i = 0; i < n; i++)
                    {
                        A[j, k] += f[j](xarray[i]) * f[k](xarray[i]);

                    }
                }
            }

            LinearSystem ls = new LinearSystem();
            VectorR coef = ls.GaussJordan(A, b);

            // Calculate the standard deviation:
            double s = 0.0;

            for (int i = 0; i < n; i++)
            {
                double s1 = 0.0;
                for (int j = 0; j < m; j++)
                {
                    s1 += coef[j] * f[j](xarray[i]);
                }
                s += (yarray[i] - s1) * (yarray[i] - s1);
            }
            sigma = Math.Sqrt(s / (n - m));

            return coef;
        }

        public static VectorR PolynomialFit(double[] xarray, double[] yarray, int m, out double sigma)
        {
            m++;
            MatrixR A = new MatrixR(m, m);
            VectorR b = new VectorR(m);
            int n = xarray.Length;

            for (int k = 0; k < m; k++)
            {
                b[k] = 0.0;
                for (int i = 0; i < n; i++)
                {
                    b[k] += Math.Pow(xarray[i], k) * yarray[i];
                }
            }

            for (int j = 0; j < m; j++)
            {
                for (int k = 0; k < m; k++)
                {
                    A[j, k] = 0.0;
                    for (int i = 0; i < n; i++)
                    {
                        A[j, k] += Math.Pow(xarray[i], j + k);
                    }
                }
            }

            LinearSystem ls = new LinearSystem();
            VectorR coef = ls.GaussJordan(A, b);

            // Calculate the standard deviation:
            double s = 0.0;

            for (int i = 0; i < n; i++)
            {
                double s1 = 0.0;
                for (int j = 0; j < m; j++)
                {
                    s1 += coef[j] * Math.Pow(xarray[i], j);
                }
                s += (yarray[i] - s1) * (yarray[i] - s1);
            }
            sigma = Math.Sqrt(s / (n - m));

            return coef;
        }

        public static double[] WeightedLinearRegression(double[] xarray, double[] yarray, double[] warray)
        {
            int n = xarray.Length;
            double xw = 0.0;
            double yw = 0.0;
            double b1 = 0.0;
            double b2 = 0.0;
            double a = 0.0;
            double b = 0.0;

            for (int i = 0; i < n; i++)
            {
                xw += xarray[i] / n;
                yw += yarray[i] / n;
            }

            for (int i = 0; i < n; i++)
            {
                b1 += warray[i] * warray[i] * yarray[i] * (xarray[i] - xw);
                b2 += warray[i] * warray[i] * xarray[i] * (xarray[i] - xw);
            }
            b = b1 / b2;
            a = yw - xw * b;

            return new double[] { a, b };
        }

        public static VectorR SimpleMovingAverage(double[] data, int n)
        {
            int m = data.Length;
            double[] sma = new double[m - n + 1];
            if (m > n)
            {
                double sum = 0.0;
                for (int i = 0; i < n; i++)
                {
                    sum += data[i];
                }

                sma[0] = sum / n;

                for (int i = 1; i <= m - n; i++)
                {
                    sma[i] = sma[i - 1] + (data[n + i - 1] - data[i - 1]) / n;
                }
            }
            return new VectorR(sma);
        }

        public static VectorR WeightedMovingAverage(double[] data, int n)
        {
            int m = data.Length;
            double[] wma = new double[m - n + 1];
            double psum = 0.0;
            double numerator = 0.0;
            double[] numerator1 = new double[m - n  + 1];
            double[] psum1 = new double[m - n + 1];

            if (m > n)
            {
                for (int i = 0; i < n; i++)
                {
                    psum += data[i];
                    numerator += (i + 1) * data[i];
                }
                psum1[0] = psum;
                numerator1[0] = numerator;
                wma[0] = 2 * numerator / n / (n + 1);

                for (int i = 1; i <= m - n; i++)
                {
                    numerator1[i] = numerator1[i - 1] + n * data[i + n - 1] - psum1[i - 1];
                    psum1[i] = psum1[i - 1] + data[i + n - 1] - data[i - 1];
                    wma[i] = 2 * numerator1[i] / n / (n + 1);
                }
            }
            return new VectorR(wma);
        }

        public static VectorR ExponentialMovingAverage(double[] data, int n)
        {
            int m = data.Length;
            double[] ema = new double[m - n + 1];
            double psum = 0.0;
            double alpha = 2.0 / n;

            if (m > n)
            {
                for (int i = 0; i < n; i++)
                {
                    psum += data[i];
                }
                ema[0] = psum / n + alpha * (data[n - 1] - psum / n);

                for (int i = 1; i <= m - n; i++)
                {
                    ema[i] = ema[i - 1] + alpha * (data[i + n - 1] - ema[i - 1]);
                }
            }
            return new VectorR(ema);
        }
    }
}
