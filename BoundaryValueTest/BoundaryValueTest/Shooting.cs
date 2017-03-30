using System;

namespace XuMath
{
    public class Shooting
    {
        private double xa = 0.0;
        private double xb = 2.0;
        private double u1 = 1.0;
        private double u2 = 2.0;
        private double h = 0.1;

        public void Shooting2()
        {
            double u = NonlinearSystem.FalsePosition(r, u1, u2, 1e-3);
            VectorR y0 = new VectorR(new double[] { 0, u });
            Console.WriteLine(u.ToString());

            for (int i = 0; i < 11; i++)
            {
                double x = 0.2 * i;
                VectorR y = ODE.MultiRungeKutta4(f, xa, y0, h, x);
                Console.WriteLine("x = {0:n2}, y[0] = {1:n6}, y[1] = {2:n6}", x, y[0], y[1]);
            }
        }

        private double r(double u)
        {
            VectorR y0 = new VectorR(new double[] { 0, u });
            VectorR y = ODE.MultiRungeKutta4(f, xa, y0, h, xb);
            return y[0] - 1.0;
        }

        private VectorR f(double x, VectorR y)
        {
            VectorR result = new VectorR(2);
            result[0] = y[1];
            result[1] = -3.0 * y[0] * y[1];
            return result;
        }


    }
}
