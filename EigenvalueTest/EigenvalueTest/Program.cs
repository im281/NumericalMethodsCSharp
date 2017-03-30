using System;
using XuMath;

namespace EigenvalueTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestJacobi();
            //TestPower();
            //TestInverse();
            //Testrayleigh();
            //TestrayleighQuotient();
            //TestTridiagonalize();
            //TestTridiagonalEigenvalues();
            TestLargeMatrix();
            Console.ReadLine();
        }

        static void TestJacobi()
        {
            MatrixR A = new MatrixR(new double[,] { { 4,3,6}, {3,7,-1}, {6,-1,9} });
            MatrixR x;
            VectorR lambda;
            Eigenvalue.Jacobi(A, 1e-10, out x, out lambda);
            Console.WriteLine("\n x = \n {0}", x);
            Console.WriteLine("\n lambda = \n {0}",lambda);
        }

        static void TestPower()
        {
            MatrixR A = new MatrixR(new double[,] { { 4, 3, 6 }, { 3, 7, -1 }, { 6, -1, 9 } });
            VectorR x;
            double lambda;
            Eigenvalue.Power(A, 1e-5, out x, out lambda);
            Console.WriteLine("\n lambda = {0} \n x= {1}", lambda, x);
        }

        static void TestInverse()
        {
            MatrixR A = new MatrixR(new double[,] { { 4, 3, 6 }, { 3, 7, -1 }, { 6, -1, 9 } });
            VectorR x;
            double lambda;
            Eigenvalue.Inverse(A, 7, 1e-5, out x, out lambda);
            Console.WriteLine("\n lambda = {0} \n x= {1}", lambda, x);
        }

        static void Testrayleigh()
        {
            MatrixR A = new MatrixR(new double[,] { { 4, 3, 6 }, { 3, 7, -1 }, { 6, -1, 9 } });
            VectorR x;
            double lambda;
            Eigenvalue.Rayleigh(A, 1e-8, out x, out lambda);
            Console.WriteLine("\n lambda = {0} \n x= {1}", lambda, x);
        }

        static void TestrayleighQuotient()
        {
            MatrixR A = new MatrixR(new double[,] { { 4, 3, 6 }, { 3, 7, -1 }, { 6, -1, 9 } });
            VectorR x;
            double lambda;
            Eigenvalue.RayleighQuotient(A, 1e-8, 1, out x, out lambda);
            Console.WriteLine("\n Results for an initial vector filled with random numbers:");
            Console.WriteLine(" lambda = {0} \n x= {1}", lambda, x);

            x = new VectorR(3);
            lambda = 0.0;
            Eigenvalue.RayleighQuotient(A, 1e-8, 2, out x, out lambda);
            Console.WriteLine("\n\n Results for an  initial vector generated from the Rayleigh method:");
            Console.WriteLine(" lambda = {0} \n x= {1}", lambda, x);
        }

        static void TestTridiagonalize()
        {
            MatrixR A = new MatrixR(new double[,]{{ 5, 1, 2, 2, 4 },
                                                  { 1, 1, 2, 1, 0},
                                                  { 2, 2, 0, 2, 1},
                                                  { 2, 1, 2, 1, 2},
                                                  { 4, 0, 1, 2, 4}});

            MatrixR V = Eigenvalue.Tridiagonalize(A);
            MatrixR T = Eigenvalue.SetTridiagonalMatrix();
            Console.WriteLine("\n Tridiagonal matrix T = \n{0}", T);
            Console.WriteLine("\n Transform matrix V = \n{0}", V);
        }

        static void TestTridiagonalEigenvalues()
        {
            MatrixR A = new MatrixR(new double[,]{{ 5, 1, 2, 2, 4 },
                                                  { 1, 1, 2, 1, 0},
                                                  { 2, 2, 0, 2, 1},
                                                  { 2, 1, 2, 1, 2},
                                                  { 4, 0, 1, 2, 4}});
            int nn = 5;
            MatrixR xx = new MatrixR(A.GetCols(), nn);
            MatrixR V = Eigenvalue.Tridiagonalize(A);
            double[] lambda = Eigenvalue.TridiagonalEigenvalues(nn);
            for (int i = 0; i < nn; i++)
            {
                double s = lambda[i] * 1.001;
                double lam;
                VectorR x = Eigenvalue.TridiagonalEigenvector(s, 1e-8, out lam);
                for (int j = 0; j < A.GetCols(); j++)
                    xx[j, i] = x[j];
            }
            xx = V * xx;

            Console.WriteLine("\n Results from the tridiagonalization method:");
            Console.WriteLine("\n Eigenvalues: \n ({0,10:n6}  {1,10:n6}  {2,10:n6}  {3,10:n6}  {4,10:n6})", lambda[0],lambda[1],lambda[2],lambda[3],lambda[4]);
            Console.WriteLine("\n Eigenvectors:");
            for (int i = 0; i < 5; i++)
            {
                Console.WriteLine(" ({0,10:n6}  {1,10:n6}  {2,10:n6}  {3,10:n6}  {4,10:n6})", xx[i,0],xx[i,1],xx[i,2],xx[i,3],xx[i,4]);
            }



            A = new MatrixR(new double[,]{{ 5, 1, 2, 2, 4 },
                                          { 1, 1, 2, 1, 0},
                                          { 2, 2, 0, 2, 1},
                                          { 2, 1, 2, 1, 2},
                                          { 4, 0, 1, 2, 4}});

            MatrixR xm;
            VectorR lamb;
            Eigenvalue.Jacobi(A, 1e-8, out xm, out lamb);

            Console.WriteLine("\n\n Results from the Jacobi method:");
            Console.WriteLine("\n Eigenvalues: \n ({0,10:n6}  {1,10:n6}  {2,10:n6}  {3,10:n6}  {4,10:n6})", lamb[4], lamb[3], lamb[2], lamb[1], lamb[0]);
            Console.WriteLine("\n Eigenvectors:");
            for (int i = 0; i < 5; i++)
            {
                Console.WriteLine(" ({0,10:n6}  {1,10:n6}  {2,10:n6}  {3,10:n6}  {4,10:n6})", xm[i, 4], xm[i, 3], xm[i, 2], xm[i, 1], xm[i, 0]);
            }
        }

        static void TestLargeMatrix()
        {
            int n = 200;
            Eigenvalue.Alpha = new double[n];
            Eigenvalue.Beta = new double[n-1];
            Eigenvalue.Alpha[0] = 3.0;
            for (int i = 0; i < n-1; i++)
            {
                Eigenvalue.Alpha[i+1] = 3.0;
                Eigenvalue.Beta[i] = -1.5;
            }

            double[] lambda = Eigenvalue.TridiagonalEigenvalues(20);
            Console.WriteLine("\n 20 smallest eigenvalues of a 200 x 200 tridiagonal matrix:\n");
            Console.WriteLine(" {0,8:n6}  {1,8:n6}  {2,8:n6}  {3,8:n6}  {4,8:n6}", lambda[0], lambda[1], lambda[2], lambda[3], lambda[4]);
            Console.WriteLine(" {0,8:n6}  {1,8:n6}  {2,8:n6}  {3,8:n6}  {4,8:n6}", lambda[5], lambda[6], lambda[7], lambda[8], lambda[9]);
            Console.WriteLine(" {0,8:n6}  {1,8:n6}  {2,8:n6}  {3,8:n6}  {4,8:n6}", lambda[10], lambda[11], lambda[12], lambda[13], lambda[14]);
            Console.WriteLine(" {0,8:n6}  {1,8:n6}  {2,8:n6}  {3,8:n6}  {4,8:n6}", lambda[15], lambda[16], lambda[17], lambda[18], lambda[19]);
        }
    }
}
