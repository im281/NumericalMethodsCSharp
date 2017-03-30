using System;
using System.Collections.Generic;
using XuMath;

namespace RealVectorTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestEquals();
            //TestNorm();
            //TestMathOperators();
            TestMultiplications();
            Console.ReadLine();
        }

        static void TestEquals()
        {
            double[] a1 = new double[3] { 1.0, 2.0, 3.0 };
            double[] a2 = new double[3] { 4.0, 5.0, 6.0 };
            VectorR v1 = new VectorR(a1);
            VectorR v2 = new VectorR(a2);
            bool b = v1 == v2;
            Console.WriteLine("\n b = {0}", b);            
        }

        static void TestNorm()
        {
            VectorR v = new VectorR(10);
            for (int i = 0; i < 10; i++)
            {
                v[i] = 0.5 * i;
            }
            double result = v.GetNorm();
            Console.WriteLine("\n Norm of the Vector = {0}", result);
        }

        static void TestMathOperators()
        {
            VectorR v1 = new VectorR(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 });
            VectorR v2 = new VectorR(new double[] { 6.0, 7.0, 8.0, 9.0, 10.0 });
            double d = 20;
            Console.WriteLine("\n v1 = {0}", v1);
            Console.WriteLine(" v2 = {0}", v2);
            Console.WriteLine(" d = {0}", d);
            Console.WriteLine(" v2 + v1 = {0}", (v2 + v1));
            Console.WriteLine(" v2 - v1 = {0}", (v2 - v1));
            Console.WriteLine(" v1 * d = {0}", (v1 * d));
            Console.WriteLine(" v1 / d = {0}", (v1 / d));
        }

        static void TestMultiplications()
        {
            VectorR v1 = new VectorR(new double[] { 1.0, 2.0, -1.0 });
            VectorR v2 = new VectorR(new double[] { 0.0, 1.0, 1.0 });
            VectorR v3 = new VectorR(new double[] { 1.0, -1.0, 0.0 });
            Console.WriteLine("\n v1 = {0}", v1);
            Console.WriteLine(" v2 = {0}", v2);
            Console.WriteLine(" v3 = {0}", v3);
            Console.WriteLine(" Dot product of v1 and v2 = {0}", VectorR.DotProduct(v1, v2));
            Console.WriteLine(" Cross product of v1 and v2 = {0}", VectorR.CrossProduct(v1, v2).ToString());
            Console.WriteLine(" Triple scalar product of v1, v2, and v3 = {0}", VectorR.TriScalarProduct(v1, v2, v3));
            Console.WriteLine(" Triple vector product of v1, v2, and v3 = {0}", VectorR.TriVectorProduct(v1, v2, v3));
        }
    }
}
