using System;
using XuMath;

namespace RealMatrixTest
{
    class Program
    {
        static void Main(string[] args)
        {
            TestMatrixR();
            Console.ReadLine();
        }

        static void TestMatrixR()
        {
            // Create a matrix using a 2D double array:
            MatrixR m1 = new MatrixR(new double[3, 3] { {1, 1, 1}, 
                                                        {1, 2, 3}, 
                                                        {1, 3, 6}});
            // Create a matrix by directly defining its elements:
            MatrixR m2 = new MatrixR(3, 3);
            m2[0, 0] = 8; m2[0, 1] = 1; m2[0, 2] = 6;
            m2[1, 0] = 3; m2[1, 1] = 5; m2[1, 2] = 7;
            m2[2, 0] = 4; m2[2, 1] = 9; m2[2, 2] = 2;
            VectorR v = new VectorR(new double[] { 2, 0, -1 });
            Console.WriteLine("\n Original matrix: m1 = \n{0}", m1);
            Console.WriteLine("\n Original matrix: m2 = \n{0}", m2);
            Console.WriteLine("\n v = {0}", (v));
            Console.WriteLine("\n m1 + m2 = \n{0}", (m1 + m2));
            Console.WriteLine("\n m1 - m2 = \n{0}", (m1 - m2));
            Console.WriteLine("\n m1 * m2 = \n{0}", (m1 * m2));
            Console.WriteLine("\n m2 * m1 = \n{0}", (m2 * m1));
            Console.WriteLine("\n m1 * v = {0}", MatrixR.Transform(m1, v));
            Console.WriteLine("\n v * m1 = {0}", MatrixR.Transform(v, m2));
            Console.WriteLine("\n Inverse of m1 = \n{0}", MatrixR.Inverse(m1));
            Console.WriteLine("\n Inverse of m2 = \n{0}", MatrixR.Inverse(m2));
            Console.WriteLine("\n Determinant of m1 = {0}", MatrixR.Determinant(m1));
            Console.WriteLine("\n Determinant of m2 = {0}", MatrixR.Determinant(m2));
        }
    }
}
