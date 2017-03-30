using System;
using System.Collections;
using XuMath;

namespace DistributionFunctionsTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestNormal();
            //TestExponential();
            //TestChi();
            //TestCauchy();
            //TestStudentT();
            //TestGamma();
            //TestBeta();
            //TestPoisson();
            TestBinomial();
            Console.ReadLine();
        }

        static void TestNormal()
        {            
            int nBins = 20;
            int nPoints = 1000;
            double xmin = -1;
            double xmax = 5;

            double[] rand = RandomGenerators.NextNormal(2.0, 1.0, nPoints);
            ArrayList aList = RandomGenerators.HistogramData(rand, xmin, xmax, nBins);
            double[] xdata = new double[nBins];
            double[] ydata = new double[nBins];

            double[] ydistribution = new double[nBins];
            for (int i = 0; i < nBins; i++)
            {
                xdata[i] = xmin + (i + 0.5) * (xmax - xmin) / nBins;
                ydata[i] = (double)aList[i];
                ydistribution[i] = DistributionFunctions.Normal(xdata[i], 2.0, 1.0);
            }
            double normalizeFactor = RandomGenerators.ArrayMax(ydata) / 
                                     RandomGenerators.ArrayMax(ydistribution);
            Console.WriteLine("");
            for (int i = 0; i < nBins; i++)
            {
                Console.WriteLine(" x = {0,4:n1}, Normal random data = {1,3:n0}, Normal distribution = {2,3:n0}",
                        xdata[i], ydata[i], Math.Round(ydistribution[i] * normalizeFactor, 0));
            }
        }

        static void TestExponential()
        {
            int nBins = 20;
            int nPoints = 2000;
            double xmin = 0;
            double xmax = 5;

            double[] rand = RandomGenerators.NextExponential(1.5, nPoints);
            ArrayList aList = RandomGenerators.HistogramData(rand, xmin, xmax, nBins);
            double[] xdata = new double[nBins];
            double[] ydata = new double[nBins];

            double[] ydistribution = new double[nBins];
            for (int i = 0; i < nBins; i++)
            {
                xdata[i] = xmin + (i + 0.5) * (xmax - xmin) / nBins;
                ydata[i] = (double)aList[i];
                ydistribution[i] = DistributionFunctions.Exponential(xdata[i], 1.5);
            }
            double normalizeFactor = RandomGenerators.ArrayMax(ydata) / RandomGenerators.ArrayMax(ydistribution);
            Console.WriteLine("");
            for (int i = 0; i < nBins; i++)
            {
                Console.WriteLine(" x = {0,4:n2}, Random data = {1,3:n0},  density distribution = {2,3:n0}",
                        xdata[i], ydata[i], Math.Round(ydistribution[i] * normalizeFactor, 0));
            }
        }

        static void TestChi()
        {
            int nBins = 20;
            int nPoints = 2000;
            double xmin = 0;
            double xmax = 4;

            double[] rand1 = RandomGenerators.NextChi(2, 1.0, nPoints);
            double[] rand2 = RandomGenerators.NextChiSquare(2, 1.0, nPoints);
            ArrayList aList1 = RandomGenerators.HistogramData(rand1, xmin, xmax, nBins);
            ArrayList aList2 = RandomGenerators.HistogramData(rand2, xmin, xmax, nBins);
            double[] xdata = new double[nBins];
            double[] ydata1 = new double[nBins];
            double[] ydata2 = new double[nBins];

            double[] ychi = new double[nBins];
            double[] ychisquare = new double[nBins];
            for (int i = 0; i < nBins; i++)
            {
                xdata[i] = xmin + (i + 0.5) * (xmax - xmin) / nBins;
                ydata1[i] = (double)aList1[i];
                ydata2[i] = (double)aList2[i];
                ychi[i] = DistributionFunctions.Chi(xdata[i], 2, 1.0);
                ychisquare[i] = DistributionFunctions.ChiSquare(xdata[i], 2, 1.0);
            }
            double normalizeFactor1 = RandomGenerators.ArrayMax(ydata1) / RandomGenerators.ArrayMax(ychi);
            double normalizeFactor2 = RandomGenerators.ArrayMax(ydata2) / RandomGenerators.ArrayMax(ychisquare);
            Console.WriteLine("\n Chi distribution");
            for (int i = 0; i < nBins; i++)
            {
                Console.WriteLine(" x = {0:n1}, Random data = {1,3:n0}, Density distribution = {2,3:n0}",
                        xdata[i], ydata1[i], Math.Round(ychi[i] * normalizeFactor1, 0));

            }

            Console.WriteLine("\n Chi-square distribution");
            for (int i = 0; i < nBins; i++)
            {
                Console.WriteLine(" x = {0,3:n1}, Random data = {1,3:n0}, Density distribution = {2,3:n0}",
                        xdata[i], ydata2[i], Math.Round(ychisquare[i] * normalizeFactor2, 0));
            }
        }

        static void TestCauchy()
        {
            int nBins = 20;
            int nPoints = 2000;
            double xmin = -4;
            double xmax = 4;

            double[] rand = RandomGenerators.NextCauchy(0.0, 0.5, nPoints);
            ArrayList aList = RandomGenerators.HistogramData(rand, xmin, xmax, nBins);
            double[] xdata = new double[nBins];
            double[] ydata = new double[nBins];

            double[] ydistribution = new double[nBins];
            for (int i = 0; i < nBins; i++)
            {
                xdata[i] = xmin + (i + 0.5) * (xmax - xmin) / nBins;
                ydata[i] = (double)aList[i];
                ydistribution[i] = DistributionFunctions.Cauchy(xdata[i], 0.0, 0.5);
            }
            double normalizeFactor = RandomGenerators.ArrayMax(ydata) / RandomGenerators.ArrayMax(ydistribution);
            Console.WriteLine("");
            for (int i = 0; i < nBins; i++)
            {
                Console.WriteLine(" x = {0,4:n1}, Random data = {1,3:n0},  density distribution = {2,3:n0}",
                        xdata[i], ydata[i], Math.Round(ydistribution[i] * normalizeFactor, 0));
            }
        }

        static void TestStudentT()
        {
            int nBins = 20;
            int nPoints = 2000;
            double xmin = -5;
            double xmax = 5;

            double[] rand = RandomGenerators.NextStudentT(5, nPoints);
            ArrayList aList = RandomGenerators.HistogramData(rand, xmin, xmax, nBins);
            double[] xdata = new double[nBins];
            double[] ydata = new double[nBins];

            double[] ydistribution = new double[nBins];
            for (int i = 0; i < nBins; i++)
            {
                xdata[i] = xmin + (i + 0.5) * (xmax - xmin) / nBins;
                ydata[i] = (double)aList[i];
                ydistribution[i] = DistributionFunctions.StudentT(xdata[i], 5);
            }
            double normalizeFactor = RandomGenerators.ArrayMax(ydata) / RandomGenerators.ArrayMax(ydistribution);
            Console.WriteLine("");
            for (int i = 0; i < nBins; i++)
            {
                Console.WriteLine(" x = {0,4:n1}, Random data = {1,3:n0},  density distribution = {2,3:n0}",
                        xdata[i], ydata[i], Math.Round(ydistribution[i] * normalizeFactor, 0));
            }
        }

        static void TestGamma()
        {
            int nBins = 20;
            int nPoints = 2000;
            double xmin = 0;
            double xmax = 15;

            double[] rand = RandomGenerators.NextGamma(2, 0.5, nPoints);
            ArrayList aList = RandomGenerators.HistogramData(rand, xmin, xmax, nBins);
            double[] xdata = new double[nBins];
            double[] ydata = new double[nBins];

            double[] ydistribution = new double[nBins];
            for (int i = 0; i < nBins; i++)
            {
                xdata[i] = xmin + (i + 0.5) * (xmax - xmin) / nBins;
                ydata[i] = (double)aList[i];
                ydistribution[i] = DistributionFunctions.Gamma(xdata[i], 2, 0.5);
            }
            double normalizeFactor = RandomGenerators.ArrayMax(ydata) / RandomGenerators.ArrayMax(ydistribution);
            Console.WriteLine("");
            for (int i = 0; i < nBins; i++)
            {
                Console.WriteLine(" x = {0,5:n2}, Random data = {1,3:n0},  density distribution = {2,3:n0}",
                        xdata[i], ydata[i], Math.Round(ydistribution[i] * normalizeFactor, 0));
            }
        }

        static void TestBeta()
        {
            int nBins = 10;
            int nPoints = 2000;
            double xmin = 0;
            double xmax = 1;

            double[] rand = RandomGenerators.NextBeta(2, 5, nPoints);
            ArrayList aList = RandomGenerators.HistogramData(rand, xmin, xmax, nBins);
            double[] xdata = new double[nBins];
            double[] ydata = new double[nBins];

            double[] ydistribution = new double[nBins];
            for (int i = 0; i < nBins; i++)
            {
                xdata[i] = xmin + (i + 0.5) * (xmax - xmin) / nBins;
                ydata[i] = (double)aList[i];
                ydistribution[i] = DistributionFunctions.Beta(xdata[i], 2, 5);
            }
            double normalizeFactor = RandomGenerators.ArrayMax(ydata) / RandomGenerators.ArrayMax(ydistribution);
            Console.WriteLine("");
            for (int i = 0; i < nBins; i++)
            {
                Console.WriteLine(" x = {0,4:n2}, Random data = {1,3:n0},  density distribution = {2,3:n0}",
                        xdata[i], ydata[i], Math.Round(ydistribution[i] * normalizeFactor, 0));
            }
        }

        static void TestPoisson()
        {
            int nBins = 15;
            int nPoints = 2000;

            double[] rand = RandomGenerators.NextPoisson(4.0, nPoints);
            ArrayList aList = RandomGenerators.HistogramData(rand, nBins);
            double[] xdata = new double[nBins];
            double[] ydata = new double[nBins];

            double[] ydistribution = new double[nBins];
            for (int i = 0; i < nBins; i++)
            {
                xdata[i] = i;
                ydata[i] = (double)aList[i];
                ydistribution[i] = DistributionFunctions.Poisson(i, 4.0);
            }
            double normalizeFactor = RandomGenerators.ArrayMax(ydata) / RandomGenerators.ArrayMax(ydistribution);
            Console.WriteLine("");
            for (int i = 0; i < nBins; i++)
            {
                Console.WriteLine(" x = {0,2:n0}, Random data = {1,3:n0},  density distribution = {2,3:n0}",
                        xdata[i], ydata[i], Math.Round(ydistribution[i] * normalizeFactor, 0));
            }
        }

        static void TestBinomial()
        {
            int nBins = 21;
            int nPoints = 2000;

            double[] rand = RandomGenerators.NextBinomial(20, 0.5, nPoints);
            ArrayList aList = RandomGenerators.HistogramData(rand, nBins);
            double[] xdata = new double[nBins];
            double[] ydata = new double[nBins];

            double[] ydistribution = new double[nBins];
            for (int i = 0; i < nBins; i++)
            {
                xdata[i] = i;
                ydata[i] = (double)aList[i];
                ydistribution[i] = DistributionFunctions.Binomial(i, 20, 0.5);
            }
            double normalizeFactor = RandomGenerators.ArrayMax(ydata) / RandomGenerators.ArrayMax(ydistribution);
            Console.WriteLine("");
            for (int i = 0; i < nBins; i++)
            {
                Console.WriteLine(" x = {0,2:n0}, Random data = {1,3:n0},  density distribution = {2,3:n0}",
                        xdata[i], ydata[i], Math.Round(ydistribution[i] * normalizeFactor, 0));
            }
        }
    }
}
