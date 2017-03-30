using System;
using XuMath;

namespace ComplexNumberTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestProperties();
            //TestMethods();
            //TestEquals();
            //TestAddition();
            //TestOperators();
            //TestMathFunctions();
            //TestTrigonometricFunctions();
            //TestInverseTrigonometricFunctions();
            //TestHyperbolicTrigonometricFunctions();
            TestInverseHyperbolicTrigonometricFunctions();
            Console.ReadLine();
        }

        static void TestProperties()
        {
            // Test public properties:
            Complex c = new Complex(2, 3);
            Console.WriteLine("\n Complex nember = {0}\n", + c);
            Console.WriteLine(" Test public properties:");
            Console.WriteLine(" Real part = {0}", c.Real);
            Console.WriteLine(" Imaginary part = {0}", c.Imaginary);
            Console.WriteLine(" Conjugate = {0}", c.Conjugate);
            Console.WriteLine(" Modulus = {0}", c.Modulus);
            Console.WriteLine(" Angle = {0}\n", c.Angle);
        }

        static void TestMethods()
        {
            //Test public methods:
            Complex c = new Complex(2, 3);
            Console.WriteLine(" Test public methods:");
            Console.WriteLine(" Real part = {0}", Complex.Re(c));
            Console.WriteLine(" Imaginary part = {0}", Complex.Im(c));
            Console.WriteLine(" Conjugate = {0}", Complex.Conj(c));
            Console.WriteLine(" Modulus = {0}", Complex.Mod(c));
            Console.WriteLine(" Angle = {0}\n", Complex.Arg(c));
        }

        static void TestEquals()
        {
            //Test Equals:
            Complex c1 = new Complex(2, 3);
            Complex c2 = new Complex(3, 4);
            bool b = c1 == c2;
            Console.WriteLine(b.ToString());
        }

        static void TestAddition()
        {
            Complex c1 = new Complex(2, 3);
            Complex c2 = new Complex(3, 4);
            double d = 5.0;
            Complex c3 = c1 + c2;      // result: 5 + 7i
            Complex c4 = c1 + d;       // result: 7 + 3i
            Complex c5 = d + c2;       // result: 8 + 4i

            Console.WriteLine("\n c3 = {0}\n c4 = {1} \n c5 = {2}", c3, c4, c5);
        }

        static void TestOperators()
        {
            Complex c1 = new Complex(4, 5);
            Complex c2 = new Complex(2, 3);

            Console.WriteLine("\n c1 = {0}", c1);
            Console.WriteLine(" c2 = {0}", c2);
            Console.WriteLine(" c1 + c2 = {0}", c1 + c2);
            Console.WriteLine(" c1 - c2 = {0}", c1 - c2);
            Console.WriteLine(" c1 * c2 = {0}", c1 * c2);
            Console.WriteLine(" c1 / c2 = {0}", c1 / c2);
        }

        static void TestMathFunctions()
        {
            Complex c1 = new Complex(4, 5);
            Complex c2 = new Complex(2, 3);
            Console.WriteLine("\n Test Math Functions:");
            Console.WriteLine(" c1 = {0}", c1);
            Console.WriteLine(" c2 = {0}", c2);
            Console.WriteLine(" Sqrt(c1) = {0}", Complex.Sqrt(c1));
            Console.WriteLine(" Pow(c1, c2) = {0}", Complex.Pow(c1, c2));
            Console.WriteLine(" Exp(c1) = {0}", Complex.Exp(c1));
            Console.WriteLine(" Log(c1) = {0}", Complex.Log(c1));
        }

        static void TestTrigonometricFunctions()
        {
            Console.WriteLine("\n Test complex trigonometric functions:");
            Complex c = new Complex(2, 3);
            Console.WriteLine(" c = {0}", c);
            Console.WriteLine(" Sin(c) = {0}", Complex.Sin(c));
            Console.WriteLine(" Cos(c) = {0}", Complex.Cos(c));
            Console.WriteLine(" Tan(c) = {0}", Complex.Tan(c));
        }

        static void TestInverseTrigonometricFunctions()
        {
            Console.WriteLine("\n Test complex inverse trigonometric functions:");
            Complex c = new Complex(2, 3);
            Console.WriteLine(" c = {0}", c);
            Console.WriteLine(" Asin(c) = {0}", Complex.Asin(c));
            Console.WriteLine(" Acos(c) = {0}", Complex.Acos(c));
            Console.WriteLine(" Atan(c) = {0}", Complex.Atan(c));
        }

        static void TestHyperbolicTrigonometricFunctions()
        {
            Console.WriteLine("\n Test complex hyperbolic trigonometric functions:");
            Complex c = new Complex(2, 3);
            Console.WriteLine(" c = {0}", c);
            Console.WriteLine(" Sinh(c) = {0}", Complex.Sinh(c));
            Console.WriteLine(" Cosh(c) = {0}", Complex.Cosh(c));
            Console.WriteLine(" Tanh(c) = {0}", Complex.Tanh(c));
        }

        static void TestInverseHyperbolicTrigonometricFunctions()
        {
            Console.WriteLine("\n Test complex inverse hyperbolic trigonometric functions:");
            Complex c = new Complex(2, 3);
            Console.WriteLine(" c = {0}", c);
            Console.WriteLine(" Asinh(c) = {0}", Complex.Asinh(c));
            Console.WriteLine(" Acosh(c) = {0}", Complex.Acosh(c));
            Console.WriteLine(" Atanh(c) = {0}", Complex.Atanh(c));
        }
    }
}

