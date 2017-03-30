using System;
using System.Collections;

namespace XuMath
{
    public static class RandomGenerators
    {
        private static Random rand = new Random();

        public static Random Rand
        {
            get { return rand; }
        }

        public static double NextBinomial(int n, double p)
        {
            double result = 0.0;
            for (int i = 0; i < n; i++)
            {
                if (rand.NextDouble() < p)
                {
                    result++;
                }
            }
            return result;
        }

        public static double[] NextBinomial(int n, double p, int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextBinomial(n, p);
            }
            return array;
        }

        public static double[] NextBeta(int a, int b, int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextBeta(a, b);
            }
            return array;
        }

        public static double NextBeta(int a, int b)
        {
            double gamma1 = NextGamma(a, 1);
            double gamma2 = NextGamma(b, 1);
            return gamma1 / (gamma1 + gamma2);
        }

        public static double[] NextGamma(int r, double alpha, int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextGamma(r,alpha);
            }
            return array;
        }

        public static double NextGamma(int r, double alpha)
        {
            if (r <= 0)
            {
                throw new ArgumentOutOfRangeException(
                 "r", r, "r must be > zero!");
            } 
            if (alpha <= 0.0)
            {
                throw new ArgumentOutOfRangeException(
                 "alpha", alpha, "alpha must be > zero!");
            }

            double result = 0.0;
            for (int i = 0; i < r; i++)
            {
                result += -Math.Log(rand.NextDouble()) / alpha;
            }
            return result;
        }

        public static double[] NextExponential(double alpha, int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextExponential(alpha);
            }
            return array;
        }

        public static double NextExponential(double alpha)
        {
            if (alpha <= 0.0)
            {
                throw new ArgumentOutOfRangeException(
                 "alpha", alpha, "alpha must be > zero!");
            }
            return -Math.Log(rand.NextDouble()) / alpha;
        }

        public static double NextCauchy()
        {
            return Math.Tan(Math.PI * rand.NextDouble() - 0.5);
        }

        public static double[] NextCauchy(int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextCauchy();
            }
            return array;
        }

        public static double NextCauchy(double a, double b)
        {
            if (b <= 0.0)
            {
                throw new ArgumentOutOfRangeException(
                 "b", b, "b must be positive!");
            }
            return a + b * (Math.Tan(Math.PI * rand.NextDouble() - 0.5));
        }

        public static double[] NextCauchy(double a, double b, int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextCauchy(a, b);
            }
            return array;
        }

        public static double NextStudentT(int n)
        {
            if (n <= 0)
            {
                throw new ArgumentOutOfRangeException(
                 "n", n, "n must be positive!");
            } 
            return NextNormal(0, 1) / NextChi(n);
        }

        public static double[] NextStudentT(int n, int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextStudentT(n);
            }
            return array;
        }

        public static double NextChi(int n, double sigma)
        {
            return Math.Sqrt(NextChiSquare(n, sigma) / n);
        }

        public static double[] NextChi(int n, double sigma, int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextChi(n, sigma);
            }
            return array;
        }

        public static double NextChi(int n)
        {
            return Math.Sqrt(NextChiSquare(n) / n);
        }

        public static double[] NextChi(int n, int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextChi(n);
            }
            return array;
        }

        public static double NextChiSquare(int n, double sigma)
        {
            if (n <= 0)
            {
                throw new ArgumentOutOfRangeException(
                 "n", n, "n must be positive!");
            } 
            if (sigma <= 0.0)
            {
                throw new ArgumentOutOfRangeException(
                 "sigma", sigma, "sigma must be positive!");
            }

            double result = 0.0;
            for (int i = 0; i < n; i++)
            {
                result += Math.Pow(NextNormal(0, sigma * sigma), 2);
            }
            return result;
        }

        public static double[] NextChiSquare(int n, double sigma, int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextChiSquare(n, sigma);
            }
            return array;
        }

        public static double NextChiSquare(int n)
        {
            if (n <= 0)
            {
                throw new ArgumentOutOfRangeException(
                 "n", n, "n must be positive!");
            }

            double result = 0.0;
            for (int i = 0; i < n; i++)
            {
                result += Math.Pow(NextNormal(0, 1), 2);
            }
            return result;
        }

        public static double[] NextChiSquare(int n, int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextChiSquare(n);
            }
            return array;
        }

        public static double NextNormal(double mu, double sigma)
        {
            double v1 = 0.0, v2 = 0.0, v12 = 0.0, y = 0.0;
            while (v12 >= 1.0 || v12 == 0.0)
            {
                v1 = 2.0 * rand.NextDouble() - 1.0;
                v2 = 2.0 * rand.NextDouble() - 1.0;
                v12 = v1 * v1 + v2 * v2;
            }
            y = Math.Sqrt(-2.0 * Math.Log(v12) / v12);
            
            return v1 * y * sigma + mu;
        }

        public static double[] NextNormal(double mu, double sigma, int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextNormal(mu, sigma);
            }
            return array;
        }

        public static double NextPoisson(double lambda)
        {
            if (lambda < 0.0)
            {
                throw new ArgumentOutOfRangeException("lambda", lambda, "lambda must be positive!");
            }
            int i = 0;
            double f = Math.Exp(-lambda);
            double p = f;
            double u = rand.NextDouble();
            while (f <= u)
            {
                p *= (lambda / (i + 1.0));
                f += p;
                i++;
            }
            return i;
        }

        public static double[] NextPoisson(double lambda, int nLength)
        {
            double[] array = new double[nLength];
            for (int i = 0; i < nLength; i++)
            {
                array[i] = NextPoisson(lambda);
            }
            return array;
        }

        public static ArrayList HistogramData(double[] data, double min, double max, int nBins)
        {
            ArrayList aList = new ArrayList();
            double dataSpacing = (max - min) / nBins;
            for (int i = 1; i < nBins + 1; i++)
            {
                int nCounts = 0;
                for (int j = 0; j < data.Length; j++)
                {
                    if (data[j] >= min + (i - 1) * dataSpacing && data[j] < min + i * dataSpacing)
                    {
                        nCounts++;
                    }
                }
                aList.Add((double)nCounts);
            }
            return aList;
        }

        public static ArrayList HistogramData(double[] data, int nBins)
        {
            ArrayList aList = new ArrayList();
            for (int i = 0; i < nBins + 1; i++)
            {
                int nCounts = 0;
                for (int j = 0; j < data.Length; j++)
                {
                    if (data[j] == i)
                    {
                        nCounts++;
                    }
                }
                aList.Add((double)nCounts);
            }
            return aList;
        }

        public static double ArrayMax(double[] array)
        {
            double max = array[0];
            for (int i = 1; i < array.Length; i++)
            {
                max = Math.Max(max, array[i]);
            }
            return max;
        }

        public static double ArrayMin(double[] array)
        {
            double min = array[0];
            for (int i = 1; i < array.Length; i++)
            {
                min = Math.Min(min, array[i]);
            }
            return min;
        }

        public static int[] RandomPermutation(int n)
        {
            ArrayList numbers = new ArrayList();
            int[] permutation = new int[n];

            // create a list that holds the numbser 0, 1, 2 ... nDimension
            for (int i = 0; i < n; i++)
            {
                numbers.Add(i);
            }

            // for each entry in the permutation list,
            // grab the number from a random position in the number list
            for (int i = 0; i < n; i++)
            {
                int n1 = rand.Next(numbers.Count);
                permutation[i] = (int)numbers[n1];
                numbers.RemoveAt(n1);
            }

            return permutation;
        }
    }
}
