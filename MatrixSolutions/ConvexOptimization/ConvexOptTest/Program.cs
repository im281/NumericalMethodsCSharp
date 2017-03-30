using ConvexOptimization;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConvexOptTest
{
    class Program
    {
        static void Main(string[] args)
        {

            #region matrix operations
            // Create a matrix using a 2D double array:
            //MatrixC m1 = new MatrixC(new double[3, 3] { {1, 1, 1}, 
            //                                            {1, 2, 3}, 
            //                                            {1, 3, 6}});
            //// Create a matrix by directly defining its elements:
            //MatrixC m2 = new MatrixC(3, 3);
            //m2[0, 0] = 8; m2[0, 1] = 1; m2[0, 2] = 6;
            //m2[1, 0] = 3; m2[1, 1] = 5; m2[1, 2] = 7;
            //m2[2, 0] = 4; m2[2, 1] = 9; m2[2, 2] = 2;
            //VectorC v = new VectorC(new double[] { 2, 0, -1 });
            //Console.WriteLine("\n Original matrix: m1 = \n{0}", m1);
            //Console.WriteLine("\n Original matrix: m2 = \n{0}", m2);
            //Console.WriteLine("\n v = {0}", (v));
            //Console.WriteLine("\n m1 + m2 = \n{0}", (m1 + m2));
            //Console.WriteLine("\n m1 - m2 = \n{0}", (m1 - m2));
            //Console.WriteLine("\n m1 * m2 = \n{0}", (m1 * m2));
            //Console.WriteLine("\n m2 * m1 = \n{0}", (m2 * m1));
            //Console.WriteLine("\n m1 * v = {0}", MatrixC.Transform(m1, v));
            //Console.WriteLine("\n v * m1 = {0}", MatrixC.Transform(v, m2));
            //Console.WriteLine("\n Inverse of m1 = \n{0}", MatrixC.Inverse(m1));
            //Console.WriteLine("\n Inverse of m2 = \n{0}", MatrixC.Inverse(m2));
            //Console.WriteLine("\n Determinant of m1 = {0}", MatrixC.Determinant(m1));
            //Console.WriteLine("\n Determinant of m2 = {0}", MatrixC.Determinant(m2));
            #endregion

            //Example
            MatrixC m1 = new MatrixC(new double[5, 3] { {1, 1, 1},
                                                        {1, 2, 3},
                                                        {1, 3, 6},
                                                        {1, 5, 14},
                                                        {1, 17, 2}});

            VectorC v1 = new VectorC(new double[] { 2, 3, 4, 6, 18 });

            ConvexOptimization.ConvexOptimization test = new ConvexOptimization.ConvexOptimization();

            //ADMM (NNLS)
            var r = test.MatchLibraryMatrix(m1, v1);

            //LS
            var ls = test.MatchLibraryMatrixLS(m1, v1);
        }
    }
}
