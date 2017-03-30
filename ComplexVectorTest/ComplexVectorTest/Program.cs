using System;
using XuMath;

namespace ComplexVectorTest
{
    class Program
    {
        static void Main(string[] args)
        {
            TestVectorC();
            Console.ReadLine();
        }

        static void TestVectorC()
        {
            VectorC v1 = new VectorC(new Complex[] {new Complex(1, 2), 
                                                    new Complex(2, 3), 
                                                    new Complex(3, 4),
                                                    new Complex(4, 5)});
            VectorC v2 = new VectorC(new Complex[] {new Complex(2, 1), 
                                                    new Complex(3, 2), 
                                                    new Complex(4, 3),
                                                    new Complex(5, 4)});
            Complex c = new Complex(5, 10);

            Console.WriteLine("\n v1 = {0}", v1);
            Console.WriteLine(" v2 = {0}", v2);
            Console.WriteLine(" c = {0}", c);
            Console.WriteLine(" v2 + v1 = {0}", (v2 + v1));
            Console.WriteLine(" v2 - v1 = {0}", (v2 - v1));
            Console.WriteLine(" v2 * c = {0}", (v2 * c));
            Console.WriteLine(" v2 / c = {0}", (v2 / c));
            Console.WriteLine(" Product of v1 and v2 = {0}", VectorC.Product(v1, v2));
            Console.WriteLine(" Dot product of v1 and v2 = {0}", VectorC.DotProduct(v1, v2));
        } 
    }
}
