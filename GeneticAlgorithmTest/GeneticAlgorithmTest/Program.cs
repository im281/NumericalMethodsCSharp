using System;
using XuMath;

namespace GeneticAlgorithmTest
{
    class Program
    {
        static void Main(string[] args)
        {
            TestDifferentialEvolution();
            //TestBinaryGeneticAlgorithm();

            Console.ReadLine();
        }

        static void TestBinaryGeneticAlgorithm()
        {
            //GeneticAlgorithm.BinaryGeneticAlgorithm(f);
            //MatrixR v = new MatrixR(new double[,] { { 1, 2, 3, 4 }, { 5, 6, 7, 8 }, { 9, 10, 11, 12 } });
            //MatrixR v1 = GeneticAlgorithm.MatrixReshape(v, 2, 6);
            //Console.WriteLine(v1.ToString());

            /*VectorR v = new VectorR(new double[] { 0.3, 0.1, 0.4, 0.2 });
            int[] ind = GeneticAlgorithm.VectorMeanIndex(v);
            foreach (int i in ind)
                Console.WriteLine(i.ToString());*/
            //GeneticAlgorithm.BinaryGeneticAlgorithm(f1);


            
        }

        static VectorR f1(MatrixR m)
        {
            int rows = m.GetRows();
            int cols = m.GetCols();
            MatrixR b = 10 * (m - 0.5);
            MatrixR b1 = new MatrixR(rows, cols);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    b1[i, j] = b[i, j] * b[i, j] - 10 * Math.Cos(2.0 * Math.PI * b[i, j]);
                }
            }
            VectorR v1 = new VectorR(cols);
            for (int i = 0; i < cols; i++)
            {
                v1[i] = 1.0;
            }
            return 10 * cols + MatrixR.Transform(b1, v1);
        }

        static void TestDifferentialEvolution()
        {
            VectorR bestmem;
            double bestval;
            GeneticAlgorithm.MaxIterations = 50;
            GeneticAlgorithm.MinCost = -50;
            GeneticAlgorithm.Refresh = 1;

          
            //GeneticAlgorithm.Xmax = new VectorR(new double[] { 0, -2 });
            //GeneticAlgorithm.Xmax = new VectorR(new double[] { 0.5, -1 });
            GeneticAlgorithm.DifferentialEvolution(Peaks, out bestmem, out bestval);
            Console.WriteLine("best memeber = \n{0}, \n best value = {1}", bestmem, bestval);

        }

        static double f(VectorR x)
        {
            return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
        }

        static double Peaks(VectorR x)
        {
            double z = 3 * (1 - x[0]) * (1 - x[0]) * Math.Exp(-x[0] * x[0] -
                       (x[1] + 1) * (x[1] + 1)) - 10 * (x[0] / 5 - Math.Pow(x[0], 3) -
                       Math.Pow(x[1], 5)) * Math.Exp(-x[0] * x[0] - x[1] * x[1])
                       - 1 / 3 * Math.Exp(-(x[0] + 1) * (x[0] + 1) - x[1] * x[1]);
            return z;
        }

    }
}
