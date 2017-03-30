using System;

namespace IfStatement
{
    class Program
    {
        static void Main(string[] args)
        {
            string myString;
            int stringLength;

            Console.Write("Please enter a string: ");
            myString = Console.ReadLine();
            stringLength = myString.Length;

            // Single decision:
            if (stringLength < 1)
            {
                Console.WriteLine("The string length (= {0}) is less than one.", stringLength);
            }

            // Either-or decision:
            if (stringLength != 1)
            {
                Console.WriteLine("The string length (= {0}) is not equal to one.", stringLength);
            }
            else
            {
                Console.WriteLine("The string length (= {0}) is equal to one.", stringLength);
            }

            // Multiple-case decision:
            if (stringLength < 1 || stringLength == 1)
            {
                Console.WriteLine("The string length (= {0}) is less than or equal to one.", stringLength);
            }
            else if (stringLength > 1 && stringLength <= 5)
            {
                Console.WriteLine("The string length (= {0}) is in the range from 2 to 5.", stringLength);
            }
            else if (stringLength > 5 && stringLength <= 10)
            {
                Console.WriteLine("The string length (= {0}) is in the range from 5 to 10.", stringLength);
            }
            else if (stringLength > 10 && stringLength <= 15)
            {
                Console.WriteLine("The string length (= {0}) is in the range from 11 to 15.", stringLength);
            }
            else
            {
                Console.WriteLine("The string length (= {0}) is greater than 15.", stringLength);
            }

            Console.ReadLine();
        }
    }
}
