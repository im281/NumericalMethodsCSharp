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
            Console.WriteLine("Complex nember: " + c.ToString() + "\n");
            Console.WriteLine("Test public properties:");
            Console.WriteLine("Real part: " + c.Real.ToString());
            Console.WriteLine("Imaginary part: " + c.Imaginary.ToString());
            Console.WriteLine("Conjugate: " + c.Conjugate.ToString());
            Console.WriteLine("Modulus: " + c.Modulus.ToString());
            Console.WriteLine("Angle: " + c.Angle.ToString() + "\n");
        }

        static void TestMethods()
        {            
            //Test public methods:
            Complex c = new Complex(2, 3);
            Console.WriteLine("Test public methods:");
            Console.WriteLine("Real part: " + Complex.Re(c).ToString());
            Console.WriteLine("Imaginary part: " + Complex.Im(c).ToString());
            Console.WriteLine("Conjugate: " + Complex.Conj(c).ToString());
            Console.WriteLine("Modulus: " + Complex.Mod(c).ToString());
            Console.WriteLine("Angle: " + Complex.Arg(c) + "\n");
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

            Console.WriteLine("c3 = " + c3.ToString());
            Console.WriteLine("c4 = " + c4.ToString());
            Console.WriteLine("c5 = " + c5.ToString());
        }

        static void TestOperators()
        {
            Complex c1 = new Complex(4, 5);
            Complex c2 = new Complex(2, 3);

            Console.WriteLine("c1 = " + c1.ToString());
            Console.WriteLine("c2 = " + c2.ToString());
            Console.WriteLine("c1 + c2 = " + (c1 + c2).ToString());
            Console.WriteLine("c1 - c2 = " + (c1 - c2).ToString());
            Console.WriteLine("c1 * c2 = " + (c1 * c2).ToString());
            Console.WriteLine("c1 / c2 = " + (c1 / c2).ToString());
        }

        static void TestMathFunctions()
        {
            Complex c1 = new Complex(4, 5);
            Complex c2 = new Complex(2, 3);
            Console.WriteLine("Test Math Functions:");
            Console.WriteLine("c1 = " + c1.ToString());
            Console.WriteLine("c2 = " + c2.ToString());
            Console.WriteLine("Sqrt(c1) = " + Complex.Sqrt(c1).ToString());
            Console.WriteLine("Pow(c1, c2) = " + Complex.Pow(c1, c2).ToString());
            Console.WriteLine("Exp(c1) = " + Complex.Exp(c1).ToString());
            Console.WriteLine("Log(c1) = " + Complex.Log(c1).ToString());
        }

        static void TestTrigonometricFunctions()
        {
            Console.WriteLine("Test complex trigonometric functions:");
            Complex c = new Complex(2, 3);
            Console.WriteLine("c = " + c.ToString());
            Console.WriteLine("Sin(c) = " + Complex.Sin(c).ToString());
            Console.WriteLine("Cos(c) = " + Complex.Cos(c).ToString());
            Console.WriteLine("Tan(c) = " + Complex.Tan(c).ToString());
        }

        static void TestInverseTrigonometricFunctions()
        {
            Console.WriteLine("Test complex inverse trigonometric functions:");
            Complex c = new Complex(2, 3);
            Console.WriteLine("c = " + c.ToString());
            Console.WriteLine("Asin(c) = " + Complex.Asin(c).ToString());
            Console.WriteLine("Acos(c) = " + Complex.Acos(c).ToString());
            Console.WriteLine("Atan(c) = " + Complex.Atan(c).ToString());
        }

        static void TestHyperbolicTrigonometricFunctions()
        {
            Console.WriteLine("Test complex hyperbolic trigonometric functions:");
            Complex c = new Complex(2, 3);
            Console.WriteLine("c = " + c.ToString());
            Console.WriteLine("Sinh(c) = " + Complex.Sinh(c).ToString());
            Console.WriteLine("Cosh(c) = " + Complex.Cosh(c).ToString());
            Console.WriteLine("Tanh(c) = " + Complex.Tanh(c).ToString());
        }

        static void TestInverseHyperbolicTrigonometricFunctions()
        {
            Console.WriteLine("Test complex inverse hyperbolic trigonometric functions:");
            Complex c = new Complex(2, 3);
            Console.WriteLine("c = " + c.ToString());
            Console.WriteLine("Asinh(c) = " + Complex.Asinh(c).ToString());
            Console.WriteLine("Acosh(c) = " + Complex.Acosh(c).ToString());
            Console.WriteLine("Atanh(c) = " + Complex.Atanh(c).ToString());
        }
    }
}
