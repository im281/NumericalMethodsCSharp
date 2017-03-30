using System;

namespace MethodDemo
{
    class Program
    {
        static void Main(string[] args)
        {
            Program p = new Program();
            p.DisplaySquareRoot(3);

            Console.ReadLine();
        }

        void DisplaySquareRoot(double x)
        {
            if (x < 0)
                return;
            Console.WriteLine("Square root of x = {0}", Math.Sqrt(Math.Abs(x)));
        }
    }
}
