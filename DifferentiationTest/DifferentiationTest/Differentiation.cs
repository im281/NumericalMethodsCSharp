using System;

namespace XuMath
{
    public class Differentiation
    {
        public delegate double Function(double x);
        private const double badResult = double.NaN;

        public static double Forward1(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (-3 * f(x) + 4 * f(x + h) - f(x + 2 * h)) / 2 / h;
        }

        public static double Forward1(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 3 || yindex < 0 || yindex > yarray.Length - 3 || h == 0)
                return badResult;
            return (-3 * yarray[yindex] + 4 * yarray[yindex + 1] - yarray[yindex + 2]) / 2 / h;
        }

        public static double Forward2(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (2 * f(x) - 5 * f(x + h) + 4 * f(x + 2 * h) - f(x + 3 * h)) / h / h;
        }

        public static double Forward2(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 4 || yindex < 0 || yindex > yarray.Length - 4 || h == 0)
                return badResult;
            return (2 * yarray[yindex] - 5 * yarray[yindex + 1] + 4 * yarray[yindex + 2] - yarray[yindex + 3]) / h / h;
        }

        public static double Forward3(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (-5 * f(x) + 18 * f(x + h) - 24 * f(x + 2 * h) + 14 * f(x + 3 * h) - 3 * f(x + 4 * h)) / 2 / h / h / h;
        }

        public static double Forward3(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 5 || yindex < 0 || yindex > yarray.Length - 5 || h == 0)
                return badResult;
            return (-5 * yarray[yindex] + 18 * yarray[yindex + 1] - 24 * yarray[yindex + 2] + 14 * yarray[yindex + 3] -
                     3 * yarray[yindex + 4]) / 2 / h / h / h;
        }

        public static double Forward4(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (3 * f(x) - 14 * f(x + h) + 26 * f(x + 2 * h) - 24 * f(x + 3 * h) + 11 * f(x + 4 * h) - 2 * f(x + 5 * h)) / h / h / h / h;
        }

        public static double Forward4(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 6 || yindex < 0 || yindex > yarray.Length - 6 || h == 0)
                return badResult;
            return (3 * yarray[yindex] - 14 * yarray[yindex + 1] + 26 * yarray[yindex + 2] - 24 * yarray[yindex + 3] +
                    11 * yarray[yindex + 4] - 2 * yarray[yindex + 5]) / h / h / h / h;
        }

        public static double Backward1(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (3 * f(x) - 4 * f(x - h) + f(x - 2 * h)) / 2 / h;
        }

        public static double Backward1(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 3 || yindex < 0 || yindex < 3 || h == 0)
                return badResult;
            return (3 * yarray[yindex] - 4 * yarray[yindex - 1] + yarray[yindex - 2]) / 2 / h;
        }

        public static double Backward2(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (2 * f(x) - 5 * f(x - h) + 4 * f(x - 2 * h) - f(x - 3 * h)) / h / h;
        }

        public static double Backward2(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 4 || yindex < 0 || yindex < 4 || h == 0)
                return badResult;
            return (2 * yarray[yindex] - 5 * yarray[yindex - 1] + 4 * yarray[yindex - 2] - yarray[yindex - 3]) / h / h;
        }

        public static double Backward3(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (5 * f(x) - 18 * f(x - h) + 24 * f(x - 2 * h) - 14 * f(x - 3 * h) + 3 * f(x - 4 * h)) / 2 / h / h / h;
        }

        public static double Backward3(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 5 || yindex < 0 || yindex < 5 || h == 0)
                return badResult;
            return (5 * yarray[yindex] - 18 * yarray[yindex - 1] + 24 * yarray[yindex - 2] - 14 * yarray[yindex - 3] +
                     3 * yarray[yindex - 4]) / 2 / h / h / h;
        }

        public static double Backward4(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (3 * f(x) - 14 * f(x - h) + 26 * f(x - 2 * h) - 24 * f(x - 3 * h) + 11 * f(x - 4 * h) - 2 * f(x - 5 * h)) / h / h / h / h;
        }

        public static double Backward4(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 6 || yindex < 0 || yindex < 6 || h == 0)
                return badResult;
            return (3 * yarray[yindex] - 14 * yarray[yindex - 1] + 26 * yarray[yindex - 2] - 24 * yarray[yindex - 3] +
                    11 * yarray[yindex - 4] - 2 * yarray[yindex - 5]) / h / h / h / h;
        }

        public static double Central1(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (f(x + h) - f(x - h)) / 2 / h;
        }

        public static double Central1(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 3 || yindex < 1 || yindex > yarray.Length - 2 || h == 0)
                return badResult;
            return (yarray[yindex + 1] - yarray[yindex - 1]) / 2 / h;
        }

        public static double Central2(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (f(x - h) - 2 * f(x) + f(x + h)) / h / h;
        }

        public static double Central2(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 3 || yindex < 1 || yindex > yarray.Length - 2 || h == 0)
                return badResult;
            return (yarray[yindex - 1] - 2 * yarray[yindex] + yarray[yindex + 1]) / h / h;
        }

        public static double Central3(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (-f(x - 2 * h) + 2 * f(x - h) - 2 * f(x + h) + f(x + 2 * h)) / 2 / h / h / h;
        }

        public static double Central3(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 5 || yindex < 2 || yindex > yarray.Length - 3 || h == 0)
                return badResult;
            return (-yarray[yindex - 2] + 2 * yarray[yindex - 1] - 2 * yarray[yindex + 1] + yarray[yindex + 2]) / 2 / h / h / h;
        }

        public static double Central4(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (f(x - 2 * h) - 4 * f(x - h) + 6 * f(x) - 4 * f(x + h) + f(x + 2 * h)) / h / h / h / h;
        }

        public static double Central4(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 5 || yindex < 2 || yindex > yarray.Length - 3 || h == 0)
                return badResult;
            return (yarray[yindex - 2] - 4 * yarray[yindex - 1] + 6 * yarray[yindex] -
                    4 * yarray[yindex + 1] + yarray[yindex + 2]) / h / h / h / h;
        }

        public static double ExtendedCentral1(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (f(x - 2 * h) - 8 * f(x - h) + 8 * f(x + h) - f(x + 2 * h)) / 12 / h;
        }

        public static double ExtendedCentral1(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 5 || yindex < 2 || yindex > yarray.Length - 3 || h == 0)
                return badResult;
            return (yarray[yindex - 2] - 8 * yarray[yindex - 1] + 8 * yarray[yindex + 1] - yarray[yindex + 2]) / 12 / h;
        }

        public static double ExtendedCentral2(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (-f(x - 2 * h) + 16 * f(x - h) - 30 * f(x) + 16 * f(x + h) - f(x + 2 * h)) / 12 / h / h;
        }

        public static double ExtendedCentral2(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 5 || yindex < 2 || yindex > yarray.Length - 3 || h == 0)
                return badResult;
            return (-yarray[yindex - 2] + 16 * yarray[yindex - 1] - 30 * yarray[yindex] + 16 * yarray[yindex + 1] - 
                     yarray[yindex + 2]) / 12 / h / h;
        }

        public static double ExtendedCentral3(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (f(x - 3 * h) - 8 * f(x - 2 * h) + 13 * f(x - h) - 13 * f(x + h) + 
                    8 * f(x + 2 * h) - f(x + 3 * h)) / 8 / h / h / h;
        }

        public static double ExtendedCentral3(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 7 || yindex < 3 || yindex > yarray.Length - 4 || h == 0)
                return badResult;
            return (yarray[yindex - 3] - 8 * yarray[yindex - 2] + 13 * yarray[yindex - 1] - 13 * yarray[yindex + 1] + 
                    8 * yarray[yindex + 2] - yarray[yindex + 3]) / 8 / h / h / h;
        }

        public static double ExtendedCentral4(Function f, double x, double h)
        {
            h = (h == 0) ? 0.01 : h;
            return (-f(x - 3 * h) + 12 * f(x - 2 * h) - 39 * f(x - h) + 56 * f(x) - 
                    39 * f(x + h) + 12 * f(x + 2 * h) - f(x + 3 * h)) / 6 / h / h / h / h;
        }

        public static double ExtendedCentral4(double[] yarray, int yindex, double h)
        {
            if (yarray == null || yarray.Length < 7 || yindex < 3 || yindex > yarray.Length - 4 || h == 0)
                return badResult;
            return (-yarray[yindex - 3] + 12 * yarray[yindex - 2] - 39 * yarray[yindex - 1] + 56 * yarray[yindex] - 
                    39 * yarray[yindex + 1] + 12 * yarray[yindex + 2] - yarray[yindex + 3]) / 6 / h / h / h / h;
        }

        public static double Richardson1(Function f, double x, double h, string flag)
        {
            double result = badResult;
            if (flag == "Backward")
            {
                result = (4 * Backward1(f, x, h / 2) - Backward1(f, x, h)) / 3;
            }
            else if (flag == "Forward")
            {
                result = (4 * Forward1(f, x, h / 2) - Forward1(f, x, h)) / 3;
            }
            else if (flag == "Central")
            {
                result = (4 * Central1(f, x, h / 2) - Central1(f, x, h)) / 3;
            }
            return result;
        }

        public static double Richardson1(double[] yarray, int yindex, double h, string flag)
        {
            double result = badResult;
            if (flag == "Backward")
            {
                result = (4 * Backward1(yarray, yindex, h / 2) - Backward1(yarray, yindex, h)) / 3;
            }
            else if (flag == "Forward")
            {
                result = (4 * Forward1(yarray, yindex, h / 2) - Forward1(yarray, yindex, h)) / 3;
            }
            else if (flag == "Central")
            {
                result = (4 * Central1(yarray, yindex, h / 2) - Central1(yarray, yindex, h)) / 3;
            }
            return result;
        }

        public static double Richardson2(Function f, double x, double h, string flag)
        {
            double result = badResult;
            if (flag == "Backward")
            {
                result = (4 * Backward2(f, x, h / 2) - Backward2(f, x, h)) / 3;
            }
            else if (flag == "Forward")
            {
                result = (4 * Forward2(f, x, h / 2) - Forward2(f, x, h)) / 3;
            }
            else if (flag == "Central")
            {
                result = (4 * Central2(f, x, h / 2) - Central2(f, x, h)) / 3;
            }
            return result;
        }

        public static double Richardson2(double[] yarray, int yindex, double h, string flag)
        {
            double result = badResult;
            if (flag == "Backward")
            {
                result = (4 * Backward2(yarray, yindex, h / 2) - Backward2(yarray, yindex, h)) / 3;
            }
            else if (flag == "Forward")
            {
                result = (4 * Forward2(yarray, yindex, h / 2) - Forward2(yarray, yindex, h)) / 3;
            }
            else if (flag == "Central")
            {
                result = (4 * Central2(yarray, yindex, h / 2) - Central2(yarray, yindex, h)) / 3;
            }
            return result;
        }

        public static double Richardson3(Function f, double x, double h, string flag)
        {
            double result = badResult;
            if (flag == "Backward")
            {
                result = (4 * Backward3(f, x, h / 2) - Backward3(f, x, h)) / 3;
            }
            else if (flag == "Forward")
            {
                result = (4 * Forward3(f, x, h / 2) - Forward3(f, x, h)) / 3;
            }
            else if (flag == "Central")
            {
                result = (4 * Central3(f, x, h / 2) - Central3(f, x, h)) / 3;
            }
            return result;
        }

        public static double Richardson3(double[] yarray, int yindex, double h, string flag)
        {
            double result = badResult;
            if (flag == "Backward")
            {
                result = (4 * Backward3(yarray, yindex, h / 2) - Backward3(yarray, yindex, h)) / 3;
            }
            else if (flag == "Forward")
            {
                result = (4 * Forward3(yarray, yindex, h / 2) - Forward3(yarray, yindex, h)) / 3;
            }
            else if (flag == "Central")
            {
                result = (4 * Central3(yarray, yindex, h / 2) - Central3(yarray, yindex, h)) / 3;
            }
            return result;
        }

        public static double Richardson4(Function f, double x, double h, string flag)
        {
            double result = badResult;
            if (flag == "Backward")
            {
                result = (4 * Backward4(f, x, h / 2) - Backward4(f, x, h)) / 3;
            }
            else if (flag == "Forward")
            {
                result = (4 * Forward4(f, x, h / 2) - Forward4(f, x, h)) / 3;
            }
            else if (flag == "Central")
            {
                result = (4 * Central4(f, x, h / 2) - Central4(f, x, h)) / 3;
            }
            return result;
        }

        public static double Richardson4(double[] yarray, int yindex, double h, string flag)
        {
            double result = badResult;
            if (flag == "Backward")
            {
                result = (4 * Backward4(yarray, yindex, h / 2) - Backward4(yarray, yindex, h)) / 3;
            }
            else if (flag == "Forward")
            {
                result = (4 * Forward4(yarray, yindex, h / 2) - Forward4(yarray, yindex, h)) / 3;
            }
            else if (flag == "Central")
            {
                result = (4 * Central4(yarray, yindex, h / 2) - Central4(yarray, yindex, h)) / 3;
            }
            return result;
        }

        public static double Interpolation0(double[] xarray, double[] yarray, double x)
        {
            double result = badResult;
            int n = yarray.Length;
            for (int i = 1; i < n - 1; i++)
            {
                if (x > xarray[i - 1] && x < xarray[i + 1])
                {
                    result = yarray[i - 1] * (x - xarray[i]) * (x - xarray[i + 1]) / ((xarray[i - 1] - xarray[i]) * (xarray[i - 1] - xarray[i + 1])) +
                             yarray[i] * (x - xarray[i - 1]) * (x - xarray[i + 1]) / ((xarray[i] - xarray[i - 1]) * (xarray[i] - xarray[i + 1])) +
                             yarray[i + 1] * (x - xarray[i - 1]) * (x - xarray[i]) / ((xarray[i + 1] - xarray[i - 1]) * (xarray[i + 1] - xarray[i]));
                }
            }
            return result;
        }

        public static double Interpolation1(double[] xarray, double[] yarray, double x)
        {
            double result = badResult;
            int n = yarray.Length;
            for (int i = 1; i < n - 1; i++)
            {
                if (x > xarray[i - 1] && x < xarray[i + 1])
                {
                    result = yarray[i - 1] * (2 * x - xarray[i] - xarray[i + 1]) / ((xarray[i - 1] - xarray[i]) * (xarray[i - 1] - xarray[i + 1])) +
                             yarray[i] * (2 * x - xarray[i - 1] - xarray[i + 1]) / ((xarray[i] - xarray[i - 1]) * (xarray[i] - xarray[i + 1])) +
                             yarray[i + 1] * (2 * x - xarray[i - 1] - xarray[i]) / ((xarray[i + 1] - xarray[i - 1]) * (xarray[i + 1] - xarray[i]));
                }
            }
            return result;
        }

        public static double Interpolation2(double[] xarray, double[] yarray, double x)
        {
            double result = badResult;
            int n = yarray.Length;
            for (int i = 1; i < n - 1; i++)
            {
                if (x > xarray[i - 1] && x < xarray[i + 1])
                {
                    result = yarray[i - 1] * 2 / ((xarray[i - 1] - xarray[i]) * (xarray[i - 1] - xarray[i + 1])) +
                             yarray[i] * 2 / ((xarray[i] - xarray[i - 1]) * (xarray[i] - xarray[i + 1])) +
                             yarray[i + 1] * 2 / ((xarray[i + 1] - xarray[i - 1]) * (xarray[i + 1] - xarray[i]));
                }
            }
            return result;
        }
    }
}
