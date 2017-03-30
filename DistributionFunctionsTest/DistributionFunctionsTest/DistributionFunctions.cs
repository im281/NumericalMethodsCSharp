using System;
using System.Collections.Generic;

namespace XuMath
{
    public static class DistributionFunctions
    {
        public static double Binomial(int x, int n, double p)
        {
            return SpecialFunctions.Gamma(n + 1) * Math.Pow(p, x) * Math.Pow(1 - p, n - x) /
                SpecialFunctions.Gamma(x + 1) / SpecialFunctions.Gamma(n - x + 1);
        }

        public static double[] Binomial(int[] x, int n, double p)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Binomial(x[i], n, p);
            }
            return y;
        }

        public static double[] Beta(double[] x, double alpha, double beta)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Beta(x[i], alpha, beta);
            }
            return y;
        }

        public static double Beta(double x, double alpha, double beta)
        {
            double b = SpecialFunctions.Beta(alpha, beta);
            return Math.Pow(x, alpha - 1) * Math.Pow(1 - x, beta - 1) / b;
        }

        public static double[] Gamma(double[] x, int r, double alpha)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Gamma(x[i], r, alpha);
            }
            return y;
        }

        public static double Gamma(double x, int r, double alpha)
        {
            return Math.Pow(alpha, r) * Math.Pow(x, r - 1) * 
                Math.Exp(-alpha * x) / SpecialFunctions.Gamma(r);
        }

        public static double[] Exponential(double[] x, double alpha)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Exponential(x[i], alpha);
            }
            return y;
        }

        public static double Exponential(double x, double alpha)
        {
            return alpha * Math.Exp(-alpha * x);
        }

        public static double[] Cauchy(double[] x)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Cauchy(x[i]);
            }
            return y;
        }

        public static double Cauchy(double x)
        {
            return 1 / (Math.PI * (1.0 + x  * x ));
        }

        public static double[] Cauchy(double[] x, double a, double b)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Cauchy(x[i], a, b);
            }
            return y;
        }

        public static double Cauchy(double x, double a, double b)
        {
            return b / (Math.PI * (b * b + (x - a) * (x - a)));
        }

        public static double[] Chi(double[] x, int n)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Chi(x[i], n);
            }
            return y;
        }

        public static double Chi(double x, int n)
        {
            double gamma = SpecialFunctions.Gamma(n / 2.0);
            double exp = Math.Exp(-n * x * x / 2.0);
            return 2.0 * Math.Pow(n / 2.0, n / 2) * Math.Pow(x, n - 1) *
                exp / Math.Pow(2, n / 2) / gamma;
        }

        public static double[] Chi(double[] x, int n, double sigma)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Chi(x[i], n, sigma);
            }
            return y;
        }

        public static double Chi(double x, int n, double sigma)
        {
            double gamma = SpecialFunctions.Gamma(n / 2.0);
            double exp = Math.Exp(-n*x*x / 2.0 / sigma / sigma);
            return 2.0*Math.Pow(n/2.0, n / 2)*Math.Pow(x,n-1) * 
                exp / Math.Pow(2, n / 2) / gamma / Math.Pow(sigma, n);
        }
        
        public static double[] ChiSquare(double[] x, int n)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = ChiSquare(x[i], n);
            }
            return y;
        }

        public static double ChiSquare(double x, int n)
        {
            double gamma = SpecialFunctions.Gamma(n / 2.0);
            double exp = Math.Exp(-x / 2.0);
            return Math.Pow(x, n / 2 - 1) * exp / Math.Pow(2, n / 2) / gamma;
        }

        public static double[] ChiSquare(double[] x, int n, double sigma)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = ChiSquare(x[i], n, sigma);
            }
            return y;
        }

        public static double ChiSquare(double x, int n, double sigma)
        {
            double gamma = SpecialFunctions.Gamma(n / 2.0);
            double exp = Math.Exp(-x / 2.0 / sigma / sigma);
            return Math.Pow(x, n / 2 - 1) * exp / Math.Pow(2, n / 2) / 
                gamma / Math.Pow(sigma, n);
        }
        
        public static double[] StudentT(double[] x, int n)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = StudentT(x[i], n);
            }
            return y;
        }

        public static double StudentT(double x, int n)
        {
            double gamma1 = SpecialFunctions.Gamma((n + 1.0) / 2.0);
            double gamma2 = SpecialFunctions.Gamma(1.0 / 2.0);
            double gamma3 = SpecialFunctions.Gamma(n / 2.0);
            return Math.Pow(n, -0.5) * Math.Pow(1 + x * x / n, -(n + 1) / 2) * 
                gamma1 / gamma2 / gamma3;
        }

        public static double Normal(double x, double mu, double sigma)
        {
            double x1, x2;
            x1 = 1 / sigma / Math.Sqrt(2 * Math.PI);
            x2 = (x - mu) * (x - mu) / (2 * sigma * sigma);
            return x1 * Math.Exp(-x2);
        }

        public static double[] Normal(double[] x, double mu, double sigma)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Normal(x[i], mu, sigma);
            }
            return y;
        }

        public static double Poisson(int x, double lambda)
        {
            double y = Math.Exp(-lambda) * Math.Pow(lambda, x);
            return y / SpecialFunctions.Gamma(x + 1);

        }

        public static double[] Poisson(int[] x, double lambda)
        {
            double[] y = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Poisson(x[i], lambda);
            }
            return y;
        }

    }
}
