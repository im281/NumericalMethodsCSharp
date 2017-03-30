using System;

namespace ForeachLoop
{
    class Program
    {
        static void Main(string[] args)
        {
            string[] students = { "Anna", "Betty", "Tyler" };
            foreach (string student in students)
            {
                Console.WriteLine(student);
            }
            Console.ReadLine();
        }
    }
}
