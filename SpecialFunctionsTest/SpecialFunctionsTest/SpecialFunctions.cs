using System;

namespace XuMath
{
    public class SpecialFunctions
    {
        public SpecialFunctions()
        {
        }

        public static double Gamma(double x)
        {
            const int g = 7;
            double[] coef = new double[9]{0.99999999999980993, 676.5203681218851,
                                          -1259.1392167224028, 771.32342877765313, 
                                          -176.61502916214059, 12.507343278686905, 
                                          -0.13857109526572012, 9.9843695780195716e-6,
                                          1.5056327351493116e-7};
            if (x < 0.5)
            {
                return Math.PI / (Math.Sin(Math.PI * x) * Gamma(1.0 - x));
            }
            x -= 1.0;
            double y = coef[0];
            for (int i = 1; i < g + 2; i++)
            {
                y += coef[i] / (x + 1.0 * i);
            }
            double z = x + (g + 0.5);
            return Math.Sqrt(2 * Math.PI) * Math.Pow(z, x + 0.5) * Math.Exp(-z) * y;
        }

        public static double Beta(double x, double y)
        {
            return Gamma(x) * Gamma(y) / Gamma(x + y);
        }

        public static double Erf(double x)
        {
            if (Math.Abs(x) > 2.2)
                return 1.0 - Erfc(x);
            double sum = 0.0;
            double sum0 = 0.0;
            int i = 0;
            do
            {
                sum0 = sum;
                sum += Math.Pow(-1, i) * Math.Pow(x, 2 * i + 1) / Gamma(i + 1) / (2 * i + 1);
                i++;
            }
            while (Math.Abs((sum - sum0) / sum0) > 1e-12);
            return 2.0 * sum / Math.Sqrt(Math.PI);
        }

        public static double Erfc(double x)
        {
            if (Math.Abs(x) <= 2.2)
                return 1.0 - Erf(x);
            if (x < 0)
                return 2.0 - Erfc(-x);
            double x1 = Math.Exp(-x * x) / Math.Sqrt(Math.PI);
            double c1 = 1.0;
            double c2 = x;
            double c3 = x;
            double c4 = x * x + 0.5;
            double c5 = 0.0;
            double c6 = c2 / c4;
            double c7 = 1.0;
            double c8 = 0.0;

            do
            {
                c8 = c1 * c7 + c2 * x;
                c1 = c2;
                c2 = c8;
                c8 = c3 * c7 + c4 * x;
                c3 = c4;
                c4 = c8;
                c7 += 0.5;
                c5 = c6;
                c6 = c2 / c4;
            }
            while (Math.Abs(c5 - c6) / c6 > 1.0e-12);
            return c6 * x1;
        }

        public static double Si(double x)
        {
            double sum = 0.0;
            double term = 0.0;
            int i = 0;
            do
            {
                term = Math.Pow(-1, i) * Math.Pow(x, 2 * i + 1) / (2 * i + 1) / Gamma(2 * i + 2);
                sum += term;
                i++;
            }
            while (Math.Abs(term) > 1.0e-12);
            return sum;
        }

        public static double Ci(double x)
        {
            double sum = 0.0;
            double term = 0.0;
            int i = 1;
            do
            {
                term = Math.Pow(-1, i) * Math.Pow(x, 2 * i) / (2 * i) / Gamma(2 * i + 1);
                sum += term;
                i++;
            }
            while (Math.Abs(term) > 1.0e-12);
            return 0.5772156649 + Math.Log(x) + sum;
        }

        public static double Laguerre(double x, int n)
        {
            double L0 = 1;
            double L1 = -x + 1;
            double L2 = (x * x - 4 * x + 2) / 2;
            int i = 1;
            if (n < 0)
                return -1;
            if (n == 0)
                return L0;
            else if (n == 1)
                return L1;
            else
            {
                while (i < n)
                {
                    L2 = ((2.0 * i + 1.0 - x) * L1 - i * L0) / (i + 1);
                    L0 = L1;
                    L1 = L2;
                    i++;
                }
                return L2;
            }
        }

        public static double Hermite(double x, int n)
        {
            double H0 = 1.0;
            double H1 = 2 * x;
            double H2 = 4 * x * x - 2;
            int i = 1;
            if (n < 0)
                return -1;
            if (n == 0)
                return H0;
            else if (n == 1)
                return H1;
            else
            {
                while (i < n)
                {
                    H2 = 2.0 * x * H1 - 2.0 * i * H0;
                    H0 = H1;
                    H1 = H2;
                    i++;
                }
                return H2;
            }
        }

        public static double ChebyshevT(double x, int n)
        {
            double T0 = 1.0;
            double T1 = x;
            double T2 = 2 * x * x - 1;
            int i = 1;
            if (n < 0)
                return -1;
            if (n == 0)
                return T0;
            else if (n == 1)
                return T1;
            else
            {
                while (i < n)
                {
                    T2 = 2.0 * x * T1 - T0;
                    T0 = T1;
                    T1 = T2;
                    i++;
                }
                return T2;
            }
        }

        public static double ChebyshevU(double x, int n)
        {
            double U0 = 1.0;
            double U1 = 2 * x;
            double U2 = 4 * x * x - 1;
            int i = 1;
            if (n < 0)
                return -1;
            if (n == 0)
                return U0;
            else if (n == 1)
                return U1;
            else
            {
                while (i < n)
                {
                    U2 = 2.0 * x * U1 - U0;
                    U0 = U1;
                    U1 = U2;
                    i++;
                }
                return U2;
            }
        }

        public static double Legendre(double x, int n)
        {
            double P0 = 1.0;
            double P1 = x;
            double P2 = (3.0 * x * x - 1) / 2.0;
            int i = 1;
            if (n < 0)
                return -1;
            if (n == 0)
                return P0;
            else if (n == 1)
                return P1;
            else
            {
                while (i < n)
                {
                    P2 = 2.0 * x * P1 - P0 - (x * P1 - P0) / (i + 1);
                    P0 = P1;
                    P1 = P2;
                    i++;
                }
                return P2;
            }
        }

        public static double BesselJ(double x, double v)
        {
            double sum = 0.0;
            double term = 0.0;
            int i = 0;
            do
            {
                term = Math.Pow(-1, i) * Math.Pow(0.5 * x, 2 * i + v) / Gamma(i + 1) / Gamma(i + v + 1);
                sum += term;
                i++;
            }
            while (Math.Abs(term) > 1.0e-12);
            return sum;
        }

        public static double BesselY(double x, double v)
        {
            return (BesselJ(x, v) * Math.Cos(v * Math.PI) - BesselJ(x, -v)) / Math.Sin(v * Math.PI);
        }
    }
}
