using System;

namespace XuMath
{
    public class ODE
    {
        public delegate double Function(double x, double y);

        public static double Euler(Function f, double x0, double y0, double h, double x)
        {
            double y;
            if (x < x0)
            {
                return double.NaN;
            }
            if (x == x0)
            {
                return y0;
            }

            do
            {
                if (h > x - x0)
                    h = x - x0;
                y = y0 + h * f(x0, y0);
                y0 = y;
                x0 += h;
            }
            while (x0 < x);
            return y;
        }

        public static double RungeKutta2(Function f, double x0, double y0, double h, double x)
        {
            double y;
            if (x < x0)
            {
                return double.NaN;
            }
            if (x == x0)
            {
                return y0;
            }

            do
            {
                if (h > x - x0)
                    h = x - x0;
                double k1 = h * f(x0, y0);
                double k2 = h * f(x0 + 0.5 * h, y0 + 0.5 * k1);
                y = y0 + k2;

                y0 = y;
                x0 += h;
            }
            while (x0 < x);
            return y;
        }

        public static double RungeKutta4(Function f, double x0, double y0, double h, double x)
        {
            double y;
            if (x < x0)
            {
                return double.NaN;
            }
            if (x == x0)
            {
                return y0;
            }

            do
            {
                if (h > x - x0)
                    h = x - x0;
                double k1 = h * f(x0, y0);
                double k2 = h * f(x0 + 0.5 * h, y0 + 0.5 * k1);
                double k3 = h * f(x0 + 0.5 * h, y0 + 0.5 * k2);
                double k4 = h * f(x0 + h, y0 + k3);
                y = y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                y0 = y;
                x0 += h;
            }
            while (x0 < x);
            return y;
        }

        public static double RungeKuttaFehlberg(Function f, double x0, double y0, double x, double h, double tolerance)
        {
            double hmin = 0.0001;
            double hmax = 0.5;
            if (h > hmax)
                h = hmax;
            if (h < hmin)
                h = hmin;

            while (x0 < x)
            {
                double k1 = h * f(x0, y0);
                double k2 = h * f(x0 + 0.25 * h, y0 + 0.25 * k1);
                double k3 = h * f(x0 + 3 * h / 8, y0 + 3 * k1 / 32 + 9 * k2 / 32);
                double k4 = h * f(x0 + 12 * h / 13, y0 + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197);
                double k5 = h * f(x0 + h, y0 + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104);
                double k6 = h * f(x0 + 0.5 * h, y0 - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40);
                double error = Math.Abs(k1 / 360 - 128 * k3 / 4275 - 2197 * k4 / 75240 + k5 / 50 + 2 * k6 / 55) / h;
                double s = Math.Pow(0.5 * tolerance / error, 0.25);

                if (error < tolerance)
                {
                    y0 += 25 * k1 / 216 + 1408 * k3 / 2565 + 2197 * k4 / 4104 - 0.2 * k5;
                    x0 += h;
                }

                if (s < 0.1)
                    s = 0.1;
                if (s > 4)
                    s = 4;
                h *= s;
                if (h > hmax)
                    h = hmax;
                if (h < hmin)
                    h = hmin;
                if (h > x - x0)
                    h = x - x0;
            }
            return y0;
        }

        public delegate VectorR MultiFunction(double x, VectorR y);
        public static VectorR MultiRungeKutta4(MultiFunction f, double x0, VectorR y0, double h, double x)
        {
            if (x <= x0)
                return y0;

            int n = y0.GetSize();
            VectorR k1 = new VectorR(n);
            VectorR k2 = new VectorR(n);
            VectorR k3 = new VectorR(n);
            VectorR k4 = new VectorR(n);
            VectorR y = y0;

            do
            {
                if (h > x - x0)
                    h = x - x0;
                k1 = h * f(x0, y);
                k2 = h * f(x0 + h / 2, y + k1 / 2);
                k3 = h * f(x0 + h / 2, y + k2 / 2);
                k4 = h * f(x0 + h, y + k3);
                y = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                x0 += h;
            }
            while (x0 < x);

            return y;
        }

        public static VectorR MultiRungeKuttaFehlberg(MultiFunction f, double x0, VectorR y0, double x, double h, double tolerance)
        {
            if (x <= x0)
                return y0;

            double hmin = 0.0001;
            double hmax = 0.5;
            if (h > hmax)
                h = hmax;
            if (h < hmin)
                h = hmin;

            int n = y0.GetSize();
            VectorR k1 = new VectorR(n);
            VectorR k2 = new VectorR(n);
            VectorR k3 = new VectorR(n);
            VectorR k4 = new VectorR(n);
            VectorR k5 = new VectorR(n);
            VectorR k6 = new VectorR(n);

            while (x0 < x)
            {

                k1 = h * f(x0, y0);
                k2 = h * f(x0 + 0.25 * h, y0 + 0.25 * k1);
                k3 = h * f(x0 + 3 * h / 8, y0 + 3 * k1 / 32 + 9 * k2 / 32);
                k4 = h * f(x0 + 12 * h / 13, y0 + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197);
                k5 = h * f(x0 + h, y0 + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104);
                k6 = h * f(x0 + 0.5 * h, y0 - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40);
                
                VectorR e1 = (k1 / 360 - 128 * k3 / 4275 - 2197 * k4 / 75240 + k5 / 50 + 2 * k6 / 55) / h;
                double error = Math.Sqrt(Math.Abs(VectorR.DotProduct(e1, e1)) / n);
                double s = Math.Pow(0.5 * tolerance / error, 0.25);

                if (error < tolerance)
                {
                    y0 = y0 + 25 * k1 / 216 + 1408 * k3 / 2565 + 2197 * k4 / 4104 - 0.2 * k5;
                    x0 += h;
                }

                if (s < 0.1)
                    s = 0.1;
                if (s > 4)
                    s = 4;
                h *= s;
                if (h > hmax)
                    h = hmax;
                if (h < hmin)
                    h = hmin;
                if (h > x - x0)
                    h = x - x0;
            }
            return y0;
        }

    }
}
