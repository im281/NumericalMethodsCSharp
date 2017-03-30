using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using XuMath;

namespace LinearSystemTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestGaussJordan();
            TestLU();
            //TestIterations();
            Console.ReadLine();
        }

        static void TestGaussJordan()
        {
            LinearSystem ls = new LinearSystem();
            MatrixR A = new MatrixR(new double[3, 3] { { 2, 1, -1 }, { -3, -1, 2 }, { -2, 1, 2 } });
            VectorR b = new VectorR(new double[3] { 8, -11, -3 });
            VectorR x = ls.GaussJordan(A, b);
            Console.WriteLine("Solution x = {0}", x);
        }

        static void TestLU()
        {
            LinearSystem ls = new LinearSystem();
            MatrixR A = new MatrixR(new double[3, 3] { { 2, 1, -1 }, { -3, -1, 2 }, { -2, 1, 2 } });
            VectorR b = new VectorR(new double[3] { 8, -11, -3 });
            MatrixR AA = A.Clone();
            MatrixR BB = A.Clone();
            double d = ls.LUCrout(A, b);
            MatrixR inv = ls.LUInverse(AA);
            Console.WriteLine("\n Inverse of A = \n {0}", (inv));
            Console.WriteLine("\n Solution of the equations = {0}", b);
            Console.WriteLine("\n Determinant of A = {0}", d);
            Console.WriteLine("\n Test Inverse: BB*Inverse = \n {0}", BB * inv);
        }

        static void TestIterations()
        {
            LinearSystem ls = new LinearSystem();
            MatrixR A = new MatrixR(new double[3, 3] { { 5, 1, 2 }, { 1, 4, 1 }, { 2, 1, 3 } });
            VectorR b = new VectorR(new double[3] { 8, 6, 6 });
            MatrixR A1 = A.Clone();
            VectorR b1 = b.Clone();

            VectorR x = ls.GaussJacobi(A, b, 10, 1.0e-4);
            Console.WriteLine("\n Solusion from the Gauss-Jacobi iteration:");
            Console.WriteLine(" x[0] = {0}", x[0]);
            Console.WriteLine(" x[1] = {0}", x[1]);
            Console.WriteLine(" x[2] = {0}", x[2]);

            VectorR x1 = ls.GaussSeidel(A1, b1, 10, 1.0e-4);
            Console.WriteLine("\n Solusion from the Gauss-Seidel iteration:");
            Console.WriteLine(" x1[0] = {0}", x1[0]);
            Console.WriteLine(" x1[1] = {0}", x1[1]);
            Console.WriteLine(" x1[2] = {0}", x1[2]);
        }
    }
}
