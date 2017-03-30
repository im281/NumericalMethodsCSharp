using System;

namespace SwitchStatement
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Which major do you like: 1 = Math, 2 = English, 3 = Physics");
            Console.Write("Please enter your selection: ");
            int selection = int.Parse(Console.ReadLine());
            string mySelection;
            switch (selection)
            {
                case 1:
                    mySelection = "I like math.";
                    break;
                case 2:
                    mySelection="I like English.";
                    break;
                case 3:
                    mySelection = "I like Physics.";
                    break;
                default:
                    mySelection = "I like none of the above.";
                    break;
            }
            Console.WriteLine(mySelection);
            Console.ReadLine();
        }
    }
}
