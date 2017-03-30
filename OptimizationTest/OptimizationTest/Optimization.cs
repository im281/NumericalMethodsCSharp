using System;
using System.Collections;

namespace XuMath
{
    public class Optimization
    {
        public delegate double Function(double x);
        public static double Bisection(Function f, double xa, double xb, double tolerance)
        {
            double xm, fa, fb, fm;
            fa = f(xa);
            fb = f(xb);
            do
            {
                xm = 0.5 * (xa + xb);
                fm = f(xm);
                
                if (Derivative(f, xm) * Derivative(f, xa) > 0)
                {
                    xa = xm;
                    fm = f(xa);
                }
                else
                {
                    xb = xm;
                    fm = f(xb);
                }

            }
            while (Math.Abs(xb - xa) > tolerance);
            return xm;
        }

        private static double Derivative(Function f, double x)
        {
            double dx = (Math.Abs(x) > 1) ? 0.01 * x : 0.01;
            return (f(x + dx) - f(x - dx)) / (2.0 * dx);
        }

        public static double GoldenSearch(Function f, double xa, double xb, double tolerance)
        {
            double x1, x2, f1, f2;
            double g = 1.0 - (Math.Sqrt(5.0) - 1.0) / 2.0;
            x1 = xa + g * (xb - xa);
            x2 = xb - g * (xb - xa);
            f1 = f(x1);
            f2 = f(x2);
            do
            {
                if (f1 < f2)
                {
                    xb = x2;
                    x2 = x1;
                    x1 = xa + g * (xb - xa);
                    f2 = f1;
                    f1 = f(x1);
                }
                else
                {
                    xa = x1;
                    x1 = x2;
                    x2 = xb - g * (xb - xa);
                    f1 = f2;
                    f2 = f(x2);
                }
            }
            while (Math.Abs(xb - xa) > tolerance);
            return 0.5 * (xa + xb);
        }

        public static double Newton(Function f, double x, double tolerance)
        {
            double dx;
            double fm, f0, fp;
            double d1, d2;
            do
            {
                dx = (Math.Abs(x) > 1) ? 0.01 * x : 0.01;
                f0 = f(x);
                fm = f(x - dx);
                fp = f(x + dx);
                d1 = (fp - fm) / (2.0 * dx);
                d2 = (fp - 2.0 * f0 + fm) / dx / dx;
                x -= d1 / d2;
            }
            while (Math.Abs(d1 / d2) > tolerance);
            return x;
        }

        public static double Brent(Function f, double xa, double xb, double tolerance)
        {
            double x1 = 0;
            double x2 = 0;
            double bx = 0;
            double xd = 0;
            double xe = 0;
            double xtemp = 0;
            double fu = 0;
            double fv = 0;
            double fw = 0;
            double fx = 0;
            double p = 0;
            double q = 0;
            double r = 0;
            double xu = 0;
            double xv = 0;
            double xw = 0;
            double x = 0;
            double xm = 0;
            double tao = 0.5 * (3.0 - Math.Sqrt(5));

            bx = 0.5 * (xa + xb);
            if (xa < xb)
            {
                x1 = xa;
            }
            else
            {
                x1 = xb;
            }
            if (xa > xb)
            {
                x2 = xa;
            }
            else
            {
                x2 = xb;
            }
            xv = bx;
            xw = xv;
            x = xv;
            xe = 0.0;
            fx = f(x);
            fv = fx;
            fw = fx;

            do
            {
                xm = 0.5 * (x1 + x2);
                if (Math.Abs(xe) > tolerance)
                {
                    r = (x - xw) * (fx - fv);
                    q = (x - xv) * (fx - fw);
                    p = (x - xv) * q - (x - xw) * r;
                    q = 2 * (q - r);
                    if (q > 0)
                    {
                        p = -p;
                    }
                    q = Math.Abs(q);
                    xtemp = xe;
                    xe = xd;
                    if (!(Math.Abs(p) >= Math.Abs(0.5 * q * xtemp) | p <= q * (x1 - x) | p >= q * (x2 - x)))
                    {
                        xd = p / q;
                        xu = x + xd;
                        if (xu - x1 < tolerance * 2 | x2 - xu < tolerance * 2)
                        {
                            xd = Math.Sign(xm - x) * tolerance;
                        }
                    }
                    else
                    {
                        if (x >= xm)
                        {
                            xe = x1 - x;
                        }
                        else
                        {
                            xe = x2 - x;
                        }
                        xd = tao * xe;
                    }
                }
                else
                {
                    if (x >= xm)
                    {
                        xe = x1 - x;
                    }
                    else
                    {
                        xe = x2 - x;
                    }
                    xd = tao * xe;
                }
                if (Math.Abs(xd) >= tolerance)
                {
                    xu = x + xd;
                }
                else
                {
                    xu = x + Math.Sign(xd) * tolerance;
                }
                fu = f(xu);
                if (fu <= fx)
                {
                    if (xu >= x)
                    {
                        x1 = x;
                    }
                    else
                    {
                        x2 = x;
                    }
                    xv = xw;
                    fv = fw;
                    xw = x;
                    fw = fx;
                    x = xu;
                    fx = fu;
                }
                else
                {
                    if (xu < x)
                    {
                        x1 = xu;
                    }
                    else
                    {
                        x2 = xu;
                    }
                    if (fu <= fw | xw == x)
                    {
                        xv = xw;
                        fv = fw;
                        xw = xu;
                        fw = fu;
                    }
                    else
                    {
                        if (fu <= fv | xv == x | xv == 2)
                        {
                            xv = xu;
                            fv = fu;
                        }
                    }
                }
            }
            while (Math.Abs(x - xm) > tolerance * 2 - 0.5 * (x2 - x1));
            return x;
        }

        public delegate double MultiFunction(VectorR x);
        public static VectorR multiNewton(MultiFunction f, double[] xarray, double tolerance)
        {
            for (int i = 0; i < xarray.Length; i++)
            {
                double dx, fm, f0, fp;
                double d1, d2;
                double x = xarray[i];
                do
                {
                    dx = (Math.Abs(x) > 1) ? 0.01 * x : 0.01;
                    xarray[i] = x - dx;
                    fm = f(new VectorR(xarray));
                    xarray[i] = x + dx;
                    fp = f(new VectorR(xarray));
                    xarray[i] = x;
                    f0 = f(new VectorR(xarray));
                    d1 = (fp - fm) / (2.0 * dx);
                    d2 = (fp + fm - 2.0 * f0) / dx / dx;
                    x -= d1 / d2;
                    xarray[i] = x;
                }
                while (Math.Abs(d1 / d2) > tolerance);
            }
            return new VectorR(xarray);
        }

        public static VectorR Simplex(MultiFunction f, MatrixR x, int MaxIterations)
        {
            double reflect = 1.0;
            double expand = 2.0;
            double contract = 0.5;
            bool flag;

            VectorR Y = new VectorR(x.GetRows());
            int nv = x.GetRows() - 1;

            int iw, ib;
            double y1, y2, x0;
            VectorR x1 = new VectorR(nv);
            VectorR x2 = new VectorR(nv);
            VectorR centroid = new VectorR(nv);

            // calculate Y using x1:
            for (int i = 0; i <= nv; i++)
            {
                for (int j = 0; j < nv; j++)
                {
                    x1[j] = x[i, j];
                }
                Y[i] = f(x1);
            }

            int iteration = 0;
            do
            {
                iteration++;

                // find worst and best points:
                iw = 0;
                ib = 0;

                for (int i = 1; i <= nv; i++)
                {
                    if (Y[i] < Y[ib])
                        ib = i;
                    else if (Y[i] > Y[iw])
                        iw = i;
                }

                // calculate centriod:
                for (int i = 0; i < nv; i++)
                {
                    centroid[i] = 0;
                    for (int j = 0; j <= nv; j++)
                    {
                        if (j != iw)
                            centroid[i] += x[j, i];
                    }
                    centroid[i] /= nv;
                }

                // calculate reflected points:
                for (int i = 0; i < nv; i++)
                {
                    x1[i] = (1 + reflect) * centroid[i] - reflect * x[iw, i];
                }
                y1 = f(x1);

                if (y1 < Y[ib])
                {
                    // calculate expended points:
                    for (int i = 0; i < nv; i++)
                    {
                        x2[i] = (1 + expand) * x1[i] - expand * centroid[i];
                    }
                    y2 = f(x2);
                    if (y2 < Y[ib])
                    {
                        for (int i = 0; i < nv; i++)
                        {
                            x[iw, i] = x2[i];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < nv; i++)
                        {
                            x[iw, i] = x1[i];
                        }
                    }
                }
                else
                {
                    flag = true;
                    for (int i = 0; i <= nv; i++)
                    {
                        if (i != iw && y1 <= Y[i])
                        {
                            flag = false;
                            break;
                        }
                    }
                    if (flag)
                    {
                        if (y1 < Y[iw])
                        {
                            for (int i = 0; i < nv; i++)
                            {
                                x[iw, i] = x1[i];
                            }
                            Y[iw] = y1;
                        }

                        // calculate contracted points:
                        for (int i = 0; i < nv; i++)
                        {
                            x2[i] = contract * x[iw, i] + (1 - contract) * centroid[i];
                        }
                        y2 = f(x2);
                        if (y2 > Y[iw])
                        {
                            for (int i = 0; i < nv; i++)
                            {
                                x2[i] = x[ib, i];
                            }
                            for (int j = 0; j <= nv; j++)
                            {
                                for (int i = 0; i < nv; i++)
                                {
                                    x[j, i] = 0.5 * (x2[i] + x[j, i]);
                                }
                            }
                        }
                        else
                        {
                            for (int i = 0; i < nv; i++)
                            {
                                x[iw, i] = x2[i];
                            }
                        }
                    }
                    else
                    {
                        for (int i = 0; i < nv; i++)
                        {
                            x[iw, i] = x1[i];
                        }
                    }

                    // calculate Y[] using x1:
                    for (int i = 0; i <= nv; i++)
                    {
                        for (int j = 0; j < nv; j++)
                        {
                            x1[j] = x[i, j];
                        }
                        Y[i] = f(x1);
                    }
                }

                // find the best points:
                ib = 0;
                for (int i = 1; i <= nv; i++)
                {
                    if (Y[i] < Y[ib])
                    {
                        ib = i;
                    }
                }
                if (ib != 0)
                {
                    for (int i = 0; i < nv; i++)
                    {
                        x0 = x[0, i];
                        x[0, i] = x[ib, i];
                        x[ib, i] = x0;
                    }
                    y1 = Y[0];
                    Y[0] = Y[ib];
                    Y[ib] = y1;
                }

            }
            while (iteration < MaxIterations);
            return x.GetRowVector(0);
        }
               
        public static VectorR Anneal(MultiFunction f, double[] xarray, double Tmin, int nEquilibrium)
        {
            double T0 = 1.0;
            double T = T0;

            double e0 = f(new VectorR(xarray));
            double e1;
            double de;
            double[] xcurrent;
            Random rand = new Random();

            int j = 0;
            do
            {
                j++;

                int i = 0;
                while (i < nEquilibrium)
                {
                    i++;
                    
                    xcurrent = RandomPerturbation((double[])xarray.Clone());
                    e1 = f(new VectorR(xcurrent));
                    de = e1 - e0;

                    if (de < 0)
                    {
                        xarray = xcurrent;
                        e0 = e1;
                    }
                    else
                    {
                        if (rand.NextDouble() < Math.Exp(-de / T))
                        {
                            xarray = xcurrent;
                            e0 = e1;
                        }
                    }
                }

                T *= Math.Pow(0.9, j);
                //T *= 0.9;
            }
            while (T > Tmin);
            return new VectorR(xarray);
        }
            
        public static double[] RandomPerturbation(double[] xarray)
        {
            int[] r = RandomGenerators.RandomPermutation(xarray.Length);
            for (int i = 0; i < xarray.Length; i++)
            {
                if (r[i] == xarray.Length - 1)
                {
                    xarray[i] += RandomGenerators.NextNormal(0,1) / 100.0;
                }
            }
            return xarray;
        }




        /**********************************************************************************************
         * 
         *  Differential Evolution
         * 
         **********************************************************************************************/

        private static double minCost = -100;
        private static int nVaraibles = 2;
        private static int nPopulations = 10;
        private static int maxIterations = 200;
        private static double stepSize = 0.8;
        private static double crossover = 0.5;
        private static VectorR xmin = new VectorR(new double[] { -2, -2 });
        private static VectorR xmax = new VectorR(new double[] { 2, 2 });
        private static int strategy = 2;
        private static int refresh = 10;

        public static double MinCost
        {
            get { return minCost; }
            set { minCost = value; }
        }
        public static int NVariables
        {
            get { return nVaraibles; }
            set { nVaraibles = value; }
        }
        public static int NPopulations
        {
            get { return nPopulations; }
            set { nPopulations = value; }
        }
        public static int MaxIterations
        {
            get { return maxIterations; }
            set { maxIterations = value; }
        }
        public static double StepSize
        {
            get { return stepSize; }
            set { stepSize = value; }
        }
        public static double CrossOver
        {
            get { return crossover; }
            set { crossover = value; }
        }

        public static VectorR Xmin
        {
            get { return xmin; }
            set { xmin = value; }
        }
        public static VectorR Xmax
        {
            get { return xmax; }
            set { xmax = value; }
        }

        public static int Strategy
        {
            get { return strategy; }
            set { strategy = value; }
        }
        public static int Refresh
        {
            get { return refresh; }
            set { refresh = value; }
        }

        public static void DifferentialEvolution(MultiFunction f, out VectorR bestMember, out double bestValue)
        {
            int n = NVariables;
            int np = NPopulations * n;
            bestMember = new VectorR(n);
            bestValue = 0.0;

            if (np < 5)
                np = 5;
            if (CrossOver < 0 || CrossOver > 1)
                CrossOver = 0.5;
            if (MaxIterations <= 0)
                MaxIterations = 200;

            MatrixR population = new MatrixR(np, n);

            Random rand = new Random();
            for (int i = 0; i < np; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    population[i, j] = Xmin[j] + rand.NextDouble() * (Xmax[j] - Xmin[j]);
                }
            }

            MatrixR population0 = new MatrixR(np, n);
            VectorR values = new VectorR(np);
            VectorR bestMemit = new VectorR(n);

            int ibest = 0;
            values[0] = f(population.GetRowVector(ibest));
            bestValue = values[0];
            for (int i = 1; i < np; i++)
            {
                values[i] = f(population.GetRowVector(i));
                if (values[i] < bestValue)
                {
                    ibest = i;
                    bestValue = values[i];
                }
            }

            bestMemit = population.GetRowVector(ibest);
            double bestvalueit = bestValue;
            bestMember = bestMemit;

            MatrixR pm1 = new MatrixR(np, n);
            MatrixR pm2 = new MatrixR(np, n);
            MatrixR pm3 = new MatrixR(np, n);
            MatrixR pm4 = new MatrixR(np, n);
            MatrixR pm5 = new MatrixR(np, n);
            MatrixR bm = new MatrixR(np, n);
            MatrixR ui = new MatrixR(np, n);
            MatrixR mui = new MatrixR(np, n);
            MatrixR mpo = new MatrixR(np, n);
            int[] rot = new int[np];
            int[] rotn = new int[n];
            for (int i = 0; i < np; i++)
                rot[i] = i;
            for (int i = 0; i < n; i++)
                rotn[i] = i;
            int[] rt = new int[np];
            int[] rtn = new int[n];
            int[] a1 = new int[np];
            int[] a2 = new int[np];
            int[] a3 = new int[np];
            int[] a4 = new int[np];
            int[] a5 = new int[np];
            int[] ind = new int[4];

            int iterations = 1;
            do
            {
                population0 = population.Clone();
                ind = RandomGenerators.RandomPermutation(4);
                a1 = RandomGenerators.RandomPermutation(np);

                for (int i = 0; i < np; i++)
                {
                    rt[i] = (rot[i] + ind[0]) % np;
                    a2[i] = a1[rt[i]];
                }
                for (int i = 0; i < np; i++)
                {
                    rt[i] = (rot[i] + ind[1]) % np;
                    a3[i] = a2[rt[i]];
                }
                for (int i = 0; i < np; i++)
                {
                    rt[i] = (rot[i] + ind[2]) % np;
                    a4[i] = a3[rt[i]];
                }
                for (int i = 0; i < np; i++)
                {
                    rt[i] = (rot[i] + ind[3]) % np;
                    a5[i] = a4[rt[i]];
                }

                double[,] randArray = new double[np, n];
                for (int i = 0; i < np; i++)
                {
                    for (int j = 0; j < n; j++)
                        randArray[i, j] = rand.NextDouble();
                }

                for (int i = 0; i < np; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        pm1[i, j] = population0[a1[i], j];
                        pm2[i, j] = population0[a2[i], j];
                        pm3[i, j] = population0[a3[i], j];
                        pm4[i, j] = population0[a4[i], j];
                        pm5[i, j] = population0[a5[i], j];

                        bm[i, j] = bestMemit[j];
                        mui[i, j] = 0;

                        if (randArray[i, j] < CrossOver)
                            mui[i, j] = 1;
                    }
                }

                int st = Strategy;
                if (Strategy > 5)
                    st = Strategy - 5;
                else
                {
                    st = Strategy;
                    mui.Transpose();
                    mui = MatrixSort(mui);
                    for (int i = 0; i < np; i++)
                    {
                        int nn = (int)Math.Floor((decimal)rand.NextDouble() * n);
                        if (nn > 0)
                        {
                            for (int j = 0; j < n; j++)
                            {
                                rtn[j] = (rotn[j] + n) % n;
                                mui[j, i] = mui[rtn[j], i];
                            }
                        }
                    }
                    mui.Transpose();
                }

                for (int i = 0; i < np; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (mui[i, j] < 0.5)
                            mpo[i, j] = 1;
                    }
                }

                for (int i = 0; i < np; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (st == 1)
                        {

                            ui[i, j] = bm[i, j] + StepSize * (pm1[i, j] - pm2[i, j]);
                            ui[i, j] = population0[i, j] * mpo[i, j] + ui[i, j] * mui[i, j];
                        }
                        else if (st == 2)
                        {
                            ui[i, j] = pm3[i, j] + StepSize * (pm1[i, j] - pm2[i, j]);
                            ui[i, j] = population0[i, j] * mpo[i, j] + ui[i, j] * mui[i, j];
                        }
                        else if (st == 3)
                        {
                            ui[i, j] = population0[i, j] + StepSize * (bm[i, j] - population0[i, j]) + StepSize * (pm1[i, j] - pm2[i, j]);
                            ui[i, j] = population0[i, j] * mpo[i, j] + ui[i, j] * mui[i, j];
                        }
                        else if (st == 4)
                        {
                            ui[i, j] = bm[i, j] + StepSize * (pm1[i, j] - pm2[i, j] + pm3[i, j] - pm4[i, j]);
                            ui[i, j] = population0[i, j] * mpo[i, j] + ui[i, j] * mui[i, j];
                        }
                        else if (st == 5)
                        {
                            ui[i, j] = pm5[i, j] + StepSize * (pm1[i, j] - pm2[i, j] + pm3[i, j] - pm4[i, j]);
                            ui[i, j] = population0[i, j] * mpo[i, j] + ui[i, j] * mui[i, j];
                        }
                    }
                }
                // Select which vectors are allowed to enter the next population:
                for (int i = 1; i < np; i++)
                {
                    double temp = f(ui.GetRowVector(i));
                    if (temp <= values[i])
                    {
                        for (int j = 0; j < n; j++)
                        {
                            population[i, j] = ui[i, j];
                        }
                        values[i] = temp;
                        if (temp < bestValue)
                        {
                            bestValue = temp;
                            bestMember = ui.GetRowVector(i);
                        }
                    }
                }
                bestMemit = bestMember;

                if (Refresh > 0)
                {
                    if (iterations % Refresh == 0)
                    {
                        Console.WriteLine("\n Iterations = {0}, Best value = {1}", iterations, bestValue);
                        Console.WriteLine(" Best member = {0}", bestMember);
                    }
                }
                iterations++;
            }
            while (iterations < MaxIterations && bestValue > MinCost);
        }


        public static MatrixR MatrixSort(MatrixR m)
        {
            int rows = m.GetRows();
            int cols = m.GetCols();

            VectorR v = new VectorR(cols);
            MatrixR m1 = new MatrixR(rows, cols);
            for (int i = 0; i < cols; i++)
            {
                v = VectorSort(m.GetColVector(i));

                for (int j = 0; j < rows; j++)
                {
                    m1[j, i] = v[j];
                }
            }
            return m1;
        }

        public static VectorR VectorSort(VectorR v)
        {
            int n = v.GetSize();
            double[] d = new double[n];
            for (int i = 0; i < n; i++)
                d[i] = v[i];
            Array.Sort(d);
            return new VectorR(d);
        }

        public static MatrixR RandomPopulation(int nChromosomes, int nVariables)
        {
            Random random = new Random();
            MatrixR pop = new MatrixR(nChromosomes, nVariables);
            for (int i = 0; i < nChromosomes; i++)
            {
                for (int j = 0; j < nVariables; j++)
                {
                    pop[i, j] = random.NextDouble();
                }
            }
            return pop;
        }

        public static MatrixR BinaryPopulation(int nChromosomes, int nStates)
        {
            Random random = new Random();
            MatrixR pop = new MatrixR(nChromosomes, nStates);
            for (int i = 0; i < nChromosomes; i++)
            {
                for (int j = 0; j < nStates; j++)
                {
                    pop[i, j] = Math.Round(random.NextDouble());
                }
            }
            return pop;
        }
    }
}
