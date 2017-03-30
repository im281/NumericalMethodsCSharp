using System;
using XuMath;

namespace ComplexMatrixTest
{
    class Program
    {
        static void Main(string[] args)
        {
            TestMatrixC();
            Console.ReadLine();
        }

        static void TestMatrixC()
        {
            // Create a complex matrix using a 2D complex array:
            MatrixC m1 = new MatrixC(new Complex[,]{{new Complex(1,1),new Complex(1,2), new Complex(1,3)},
                                                    {new Complex(2,1),new Complex(2,2), new Complex(2,3)},
                                                    {new Complex(3,1),new Complex(3,2), new Complex(3,3)}});

            //Create a complex matrix by directly defining its elements:
            MatrixC m2 = new MatrixC(3, 3);
            m2[0, 0] = new Complex(1, 2);
            m2[0, 1] = new Complex(7, -3);
            m2[0, 2] = new Complex(3, 4);
            m2[1, 0] = new Complex(6, -2);
            m2[1, 1] = new Complex(0, 9);
            m2[1, 2] = new Complex(4, 7);
            m2[2, 0] = new Complex(2, 1);
            m2[2, 1] = new Complex(3, -1);
            m2[2, 2] = new Complex(3, 0);
            VectorC v = new VectorC(new Complex[] { new Complex(1, 1), 
                                                    new Complex(1, 2), 
                                                    new Complex(1, 3) });
            Console.WriteLine("\n Original matrix: m1 = \n{0}", m1);
            Console.WriteLine("\n Original matrix: m2 = \n{0}", m2);
            Console.WriteLine("\n v = {0}", (v));

            Console.WriteLine("\n m1 + m2 = \n{0}", (m1 + m2));
            Console.WriteLine("\n m1 - m2 = \n{0}", (m1 - m2));
            Console.WriteLine("\n m1 * m2 = \n{0}", (m1 * m2));
            Console.WriteLine("\n m2 * m1 = \n{0}", (m2 * m1));
            Console.WriteLine("\n m1 * v = {0}", MatrixC.Transform(m1, v));
            Console.WriteLine("\n v * m1 = {0}", MatrixC.Transform(v, m2));
            Console.WriteLine("\n Inverse of m2 = \n{0}", MatrixC.Inverse(m2));
            Console.WriteLine("\n Determinant of m1 = {0}", MatrixC.Determinant(m1));
            Console.WriteLine("\n Determinant of m2 = {0}", MatrixC.Determinant(m2));
        }

    }
}
