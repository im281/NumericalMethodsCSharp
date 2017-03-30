using System;

namespace BasicMathOperations
{
    class Program
    {
        static void Main(string[] args)
        {
            double x;
            x = 5 + 3;
            Console.WriteLine("\n Addition: \n 4 + 5 = {0}\n", x);
            x = 5 - 3;
            Console.WriteLine(" Subtraction: \n 5 - 3 = {0}\n", x);
            x = 5 * 3; ;
            Console.WriteLine(" Multiplication: \n 5 * 3 = {0}\n", x);
            x = 5 / 3.0;
            Console.WriteLine(" Division: \n 5 / 3 = {0}", x);
            Console.ReadLine();
        }
    }
}
