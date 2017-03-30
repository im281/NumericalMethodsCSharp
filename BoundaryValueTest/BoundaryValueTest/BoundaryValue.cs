using System;

namespace XuMath
{
    public class BoundaryValue
    {
        // Shooting method for the second-order differential equation: y" = f(x,y,y')
        public double xa { get; set; }
        public double xb { get; set; }
        public double ya { get; set; }
        public double yb { get; set; }
        public double u1 { get; set; }
        public double u2 { get; set; }
        public double StepSize { get; set; }
        public double xOut { get; set; }

        public ODE.MultiFunction F1 { get; set; }

        public VectorR Shooting2()
        {
            double u = NonlinearSystem.FalsePosition(ResidualFunction, u1, u2, 1e-8);
            VectorR y0 = new VectorR(new double[] { ya, u });
            return ODE.MultiRungeKutta4(F1, xa, y0, StepSize, xOut);
       }

        private double ResidualFunction(double u)
        {
            VectorR y0 = new VectorR(new double[] { ya, u });
            VectorR y = ODE.MultiRungeKutta4(F1, xa, y0, StepSize, xb);
            return y[0] - yb;
        }






        // Finite Difference method for the second-order linear differential equation: P(x)y" + Q(x)y' + R(x)y + S(x)=0
        // Note xa, xb, ya, and yb are the same as in the shooting method.

        public double va { get; set; }            // va = y'(a)
        public double vb { get; set; }            // vb = y'(b)
        public int n { get; set; }
        
        private int[] boundaryFlag = new int[] { 0, 0 };
        public int[] BoundaryFlag
        {
            get { return boundaryFlag; }
            set { boundaryFlag = value; }
        }


        public delegate VectorR MultiFunction(double x);

        public VectorR FiniteDifferenceLinear2(MultiFunction f, out double[] x)
        {
            double h = (xb - xa) / n;
            x = new double[n + 1];
            for (int i = 0; i < n + 1; i++)
                x[i] = xa + i * h;

            double[] a = new double[n + 1];
            double[] b = new double[n + 1];
            double[] c = new double[n + 1];
            double[] d = new double[n + 1];

            if (boundaryFlag[0] != 1)
            {
                a[0] = 0.0;
                b[0] = 1.0;
                c[0] = 0.0;
                d[0] = ya;
            }
            else
            {
                a[0] = 0.0;
                b[0] = f(x[0])[2] - 2 * f(x[0])[0] / h / h;
                c[0] = 2 * f(x[0])[0] / h / h;
                d[0] = 2 * h * va * (f(x[0])[0] / h / h - f(x[0])[1] / 2 / h) - f(x[0])[3];
            }
            if (boundaryFlag[1] != 1)
            {
                a[n] = 0.0;
                b[n] = 1.0;
                c[n] = 0.0;
                d[n] = yb;
            }
            else
            {
                a[n] = 2 * f(x[n])[0] / h / h;
                b[n] = f(x[n])[2] - 2 * f(x[n])[0] / h / h;
                c[n] = 0.0;
                d[n] = -2 * h * vb * (f(x[n])[0] / h / h + f(x[n])[1] / 2 / h) - f(x[n])[3];
            }

            for (int i = 1; i < n; i++)
            {
                a[i] = f(x[i])[0] / h / h - f(x[i])[1] / 2 / h;
                b[i] = f(x[i])[2] - 2 * f(x[i])[0] / h / h;
                c[i] = f(x[i])[0] / h / h + f(x[i])[1] / 2 / h;
                d[i] = -f(x[i])[3];
            }

            MatrixR A1 = new MatrixR(n + 1, n + 1);
            A1[0, 0] = b[0];
            for (int i = 1; i <= n; i++)
            {
                A1[i, i] = b[i];
                A1[i, i - 1] = a[i];
                A1[i - 1, i] = c[i - 1];
            }
            VectorR b1 = new VectorR(d);

            LinearSystem ls = new LinearSystem();
            double d1 = ls.LUCrout(A1, b1);
            
            return b1;
        }


        // Finite difference for second-order nonlinear differential equation: y" = f(x,y,y')
        public delegate double FDFunction(double x, double y, double yprime);
        public FDFunction fd { get; set; }

        public VectorR FiniteDifferenceNonlinear2(out double[] x)
        {
            double h = (xb - xa) / n;
            x = new double[n + 1];
            VectorR y = new VectorR(n + 1);
            for (int i = 0; i < n + 1; i++)
            {
                x[i] = xa + i * h;
                y[i] = 0.5 * x[i];
            }
            y = NonlinearSystem.NewtonMultiEquations(VF, y, 1e-8);
            return y;
        }

        private VectorR VF(VectorR y)
        {
            double h = (xb - xa) / n;
            double[] x = new double[n + 1];
            for (int i = 0; i < n + 1; i++)
                x[i] = xa + i * h;

            VectorR result = new VectorR(n + 1);

            if (boundaryFlag[0] != 1)
            {
                result[0] = y[0] - ya;
            }
            else
            {
                result[0] = -2 * y[0] + 2 * y[1] - h * h * fd(x[0], y[0], va) - 2 * h * va;
            }
            if (boundaryFlag[1] != 1)
            {
                result[n] = y[n] - yb;
            }
            else
            {
                result[n] = 2 * y[n - 1] - 2 * y[n] - h * h * fd(x[n], y[n], vb) + 2 * h * yb;
            }
            for (int i = 1; i < n; i++)
            {
                result[i] = y[i - 1] - 2 * y[i] + y[i + 1] - h * h * fd(x[i], y[i], (y[i + 1] - y[i - 1]) / 2 / h);
            }
            return result;
        }



        // Finite difference method for fourth-order linear differential equation: y"" = P(x)*y + Q(x)
        public double ga { get; set; }
        public double gb { get; set; }

        public VectorR FiniteDifferenceLinear4(MultiFunction f, out double[] x)
        {
            double h = (xb - xa) / n;
            double h3 = h * h * h;
            double h4 = h * h * h * h;
            x = new double[n + 1];
            for (int i = 0; i < n + 1; i++)
                x[i] = xa + i * h;

            double[] a = new double[n + 1];
            double[] b = new double[n + 1];
            double[] c = new double[n + 1];
            double[] d = new double[n + 1];
            double[] e = new double[n + 1];
            double[] t = new double[n + 1];

            if (boundaryFlag[0] != 1)
            {
                c[0] = 1.0;
                d[0] = 0.0;
                e[0] = 0.0;
                t[0] = ya;
            }
            else
            {
                c[0] = 6 - h4 * f(x[0])[0];
                d[0] = -8.0;
                e[0] = 2.0;
                t[0] = h4 * f(x[0])[1] + 2 * h3 * ga - 4 * h * va;
            }

            b[1] = -4.0;
            c[1] = 7 - h4 * f(x[1])[0];
            d[1] = -4.0;
            e[1] = 1.0;
            t[1] = h4 * f(x[1])[1] + 2 * h * va;

            for (int i = 2; i < n - 1; i++)
            {
                a[i] = 1.0;
                b[i] = -4.0;
                c[i] = 6 - h4 * f(x[i])[0];
                d[i] = -4.0;
                e[i] = 1.0;
                t[i] = h4 * f(x[i])[1];
            }

            a[n - 1] = 1.0;
            b[n - 1] = -4.0;
            c[n - 1] = 7 - h4 * f(x[n - 1])[0];
            d[n - 1] = -4.0;
            t[n - 1] = h4 * f(x[n - 1])[1] - 2 * h * vb;

            if (boundaryFlag[1] != 1)
            {
                a[n] = 0.0;
                b[n] = 0.0;
                c[n] = 1.0;
                t[n] = yb;
            }
            else
            {
                a[n] = 2.0;
                b[n] = -8.0;
                c[n] = 6 - h4 * f(x[n])[0];
                t[n] = h4 * f(x[n])[1] - 2 * h3 * gb + 4 * h * vb;
            }

            MatrixR A1 = new MatrixR(n + 1, n + 1);
            A1[0, 0] = c[0];
            A1[0, 1] = d[0];
            A1[1, 0] = b[1];
            A1[1, 1] = c[1];
            for (int i = 2; i <= n; i++)
            {
                A1[i, i - 2] = a[i];
                A1[i, i - 1] = b[i];
                A1[i, i] = c[i];
                A1[i - 1, i] = d[i - 1];
                A1[i - 2, i] = e[i - 2];
            }

            VectorR b1 = new VectorR(t);
            LinearSystem ls = new LinearSystem();
            double d1 = ls.LUCrout(A1, b1);

            return b1;
        }
    }
}
