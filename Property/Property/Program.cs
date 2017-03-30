using System;

namespace Property
{
    class Program
    {
        
        static void Main(string[] args)
        {
            /*Square square1 = new Square();
            square1.Side = 10.5;


            Square square2 = new Square();
            square2.Side = 22.5;

            Console.WriteLine(square1.Side);        // Output = 10.5
            Console.WriteLine(square2.Side);        // Output = 22.5
            
            Square square3 = new Square();
            square3.Side = -5.5;                    // Throws an overflow exception.*/

            /*Square square = new Square();
            square.Side = 10;
            Console.WriteLine(square.Area);         // Output = 100*/

            /*Square s = new Square();
            s.Side = 10;
            s.Perimeter = 4 * s.Side;               // Output = 40*/

            Customer customer = new Customer { ID = 1, Name = "Tyler", Purchases = 10.50 };

            Console.WriteLine("Customer ID = {0}, Name = {1}, Puerchases = {2}", customer.ID, customer.Name, customer.Purchases);


            Console.ReadLine();
        }
    }

    class Square
    {
        private double side;

        public double Side
        {
            get { return side; }
            set
            {
                if (value < 0)
                    throw new OverflowException();
                side = value;
            }
        }

        public double Area
        {
            get { return side * side; }
        }

        private double perimeter;
        public double Perimeter
        {
            set
            {
                perimeter = value;
                Console.WriteLine("The perimeter of the square = {0}", perimeter);
            }
        }
    }

    class Customer
    {
        public int ID { get; set; }
        public string Name { get; set; }
        public double Purchases { get; set; }
    }
}
