using System;

namespace DoLoop
{
    class Program
    {
        static void Main(string[] args)
        {
            int n = 0;
            do
            {
                Console.WriteLine(n);
                n++;
            } while (n < 10);
            Console.ReadLine();
        }
    }
}
