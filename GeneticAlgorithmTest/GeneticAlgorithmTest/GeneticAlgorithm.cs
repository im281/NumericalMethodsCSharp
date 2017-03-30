using System;
using System.Collections;

namespace XuMath
{
    public class GeneticAlgorithm
    {
        public delegate double Function(VectorR x);

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

        public static void DifferentialEvolution(Function f, out VectorR bestMember, out double bestValue)
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

                Console.WriteLine("Strategy = {0}", st);

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
                        Console.WriteLine(" Best member = \n{0}", bestMember);
                    }
                }

                iterations++;

            }
            while (iterations < MaxIterations && bestValue > MinCost);          
        }



        //====================================================================================
        //
        //  Evolutionary Algorithm
        //
        //====================================================================================

        public void EvolutionAlgorithm(Function f, out VectorR bestMember, out double bestValue)
        {
            /*HashMain hm = new HashMain();
            ArrayList aTargetParameters = hm.HashMainSearch("frmTargetParameters");
            ArrayList aDeviceParameters = hm.HashMainSearch("DeviceParameters");
            ArrayList aEAParameters = hm.HashMainSearch("EADEParameters");

            if (aTargetParameters == null || aEAParameters == null || aDeviceParameters == null)
            {
                MessageBox.Show("Please set device parameters and specify target function.");
                return;
            }

            ArrayList aListDeviceProperties = (ArrayList)aDeviceParameters[2];

            string[] sArrayStructures = (string[])aDeviceParameters[4];
            int nSubstrate = 0;
            for (int i = 1; i < sArrayStructures.Length; i++)
            {
                if (sArrayStructures[i] == "Substrate")
                    nSubstrate = i - 1;
            }

            // Target parameters:
            string[] sThickness = (string[])aTargetParameters[1];
            string[] sOptim = (string[])aTargetParameters[3];
            float fTheta0 = (float)aTargetParameters[4];
            int nFlag = (int)aTargetParameters[5];

            ArrayList aListProperties = (ArrayList)aTargetParameters[2];
            ArrayList aListTarget = (ArrayList)aTargetParameters[6];
            double[] dLambda = (double[])aListTarget[0];

            float[] fLambda = uMath.DoubleArrayToFloatArray(dLambda);

            float[] fThickness = new float[sThickness.Length - 2];
            string[] sOptimization = new string[sThickness.Length - 2];
            for (int i = 1; i < sThickness.Length - 1; i++)
            {
                fThickness[i - 1] = Convert.ToSingle(sThickness[i]);
                sOptimization[i - 1] = sOptim[i];
            }

            float[] fThickness0 = fThickness;

            ComplexFloatMatrix mNK = new ComplexFloatMatrix(fThickness.Length, fLambda.Length);
            ComplexFloatVector vNKIncidentMedium = new ComplexFloatVector(fLambda.Length);
            ComplexFloatVector vNKExitMedium = new ComplexFloatVector(fLambda.Length);

            ArrayList aListFirstLayer = (ArrayList)aListProperties[0];
            ArrayList aListLastLayer = (ArrayList)aListProperties[aListProperties.Count - 1];
            double[] fNNIncidentMedium = (double[])aListFirstLayer[1];
            double[] fKKIncidentMedium = (double[])aListFirstLayer[2];
            double[] fNNExitMedium = (double[])aListLastLayer[1];
            double[] fKKExitMedium = (double[])aListLastLayer[2];

            for (int i = 0; i < fLambda.Length; i++)
            {
                vNKIncidentMedium[i] = new ComplexFloat((float)fNNIncidentMedium[i], -(float)fKKIncidentMedium[i]);
                vNKExitMedium[i] = new ComplexFloat((float)fNNExitMedium[i], -(float)fKKExitMedium[i]);
            }

            ArrayList aNK = new ArrayList();
            ArrayList aN0 = new ArrayList();
            double dLambda0 = 550;
            float[] fN0 = new float[fThickness.Length];
            for (int i = 1; i < sThickness.Length - 1; i++)
            {
                aNK = (ArrayList)aListProperties[i];
                aN0 = (ArrayList)aListDeviceProperties[i];

                double[] fNN = (double[])aNK[1];
                double[] fKK = (double[])aNK[2];
                dLambda0 = (double)aN0[0];
                double dN0 = (double)aN0[1];
                fN0[i - 1] = (float)dN0;
                for (int j = 0; j < fLambda.Length; j++)
                {
                    mNK[i - 1, j] = new ComplexFloat((float)fNN[j], -(float)fKK[j]);
                }
            }*/

            // Set EA parameters:
            /*float fThicknessFrom = (float)aEAParameters[1];
            float fThicknessTo = (float)aEAParameters[2];
            int nOptimLayers = (int)aEAParameters[3];
            int nPopulation = (int)aEAParameters[4];
            int nGenerations = (int)aEAParameters[5];
            float fStepSize = (float)aEAParameters[6];
            int nTotalLayers = fThickness.Length;*/

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

            float[] fMeritValue = new float[nPopulation];
            for (int i = 0; i < nPopulation; i++)
            {
                fMeritValue[i] = 0.0f;
            }
            float[] fBestPopulationMember = new float[nTotalLayers];
            float[] fBestPopulationMemberInIterations = new float[nTotalLayers];
            for (int i = 0; i < nOptimLayers; i++)
            {
                fBestPopulationMember[i] = 0.0f;
                fBestPopulationMemberInIterations[i] = 0.0f;
            }

            // Evaluate the best member after initialization starting with first population member

            int nBest = 1;
            for (int j = 0; j < nTotalLayers; j++)
            {
                fThickness[j] = fPopulation[nBest, j];
            }


            ArrayList aListResult = MeritFunctionResults(nFlag, fLambda, fTheta0, fThickness,
                                   aListTarget, mNK, vNKIncidentMedium, vNKExitMedium);


            fMeritValue[0] = (float)aListResult[0];

            ArrayList aListBestResult = aListResult;
            float fBestValue = fMeritValue[0];

            for (int i = 1; i < nPopulation; i++)
            {
                for (int j = 0; j < nTotalLayers; j++)
                {
                    fThickness[j] = fPopulation[i, j];
                }
                aListResult = MeritFunctionResults(nFlag, fLambda, fTheta0, fThickness,
                                       aListTarget, mNK, vNKIncidentMedium, vNKExitMedium);

                fMeritValue[i] = (float)aListResult[0];
                if (fMeritValue[i] < fBestValue)
                {
                    nBest = i;
                    fBestValue = fMeritValue[i];
                    aListBestResult = aListResult;
                }
            }


            for (int j = 0; j < nTotalLayers; j++)
            {
                fBestPopulationMemberInIterations[j] = fPopulation[nBest, j];
            }
            fBestPopulationMember = fBestPopulationMemberInIterations;

            //****************************************************************************
            //
            //  EA minimization
            //
            //****************************************************************************

            // Set EA Step Size:
            float[,] fSigma = new float[nPopulation, nTotalLayers];
            float[,] fSigmaOld = new float[nPopulation, nTotalLayers];
            float[,] fSigmaC = new float[nPopulation, nTotalLayers];
            float[,] fPC = new float[nPopulation, nTotalLayers];
            float[,] fZi = new float[nPopulation, nTotalLayers];
            float[,] fZ0i = new float[nPopulation, nTotalLayers];
            float[,] fPR = new float[nPopulation, nTotalLayers];
            float[,] fPT = new float[nPopulation, nTotalLayers];

            for (int i = 0; i < nPopulation; i++)
            {
                for (int j = 0; j < nTotalLayers; j++)
                {
                    fSigma[i, j] = fStepSize;
                    fPC[i, j] = 0;
                }
            }

            float fTau0 = 1.0f / (float)Math.Sqrt(2 * nTotalLayers);
            float fTaoi = 1.0f / (float)Math.Sqrt(2 * Math.Sqrt(nTotalLayers));

            ArrayList aListIterations = new ArrayList();
            ArrayList aListMeritValue = new ArrayList();

            int nIterations = 1;
            while (nIterations < nGenerations)
            {
                fPopulationOld = fPopulation;
                fSigmaOld = fSigma;

                // Choose Parents:
                int[] nDR1 = new int[nPopulation];
                int[] nDR2 = new int[nPopulation];
                for (int i = 0; i < nPopulation; i++)
                {
                    nDR1[i] = (int)Math.Floor(nPopulation * rand.NextDouble());
                    nDR2[i] = (int)Math.Floor(nPopulation * rand.NextDouble());
                }

                // Discrete Recombinations:
                for (int i = 0; i < nPopulation; i++)
                {
                    for (int j = 0; j < nTotalLayers; j++)
                    {
                        fPC[i, j] = fPopulationOld[nDR1[i], j];
                    }
                }

                int[] nDummy = new int[nTotalLayers];
                for (int i = 0; i < nTotalLayers; i++)
                {
                    nDummy[i] = (int)Math.Round(rand.NextDouble());
                }

                for (int i = 0; i < nPopulation; i++)
                {
                    for (int j = 0; j < nTotalLayers; j++)
                    {
                        if (nDummy[j] == 0)
                        {
                            fPC[i, j] = fPopulationOld[nDR2[i], j];
                        }
                    }
                }

                // Intermediate Recombinations:
                for (int i = 0; i < nPopulation; i++)
                {
                    for (int j = 0; j < nTotalLayers; j++)
                    {
                        fSigmaC[i, j] = 0.5f * (fSigma[nDR1[i], j] + fSigma[nDR2[i], j]);
                    }
                }

                // Mutations:
                float fZ0 = fTau0 * (float)RandomProvider.NextNormal();
                for (int i = 0; i < nPopulation; i++)
                {
                    for (int j = 0; j < nTotalLayers; j++)
                    {
                        fZi[i, j] = fTaoi * (float)RandomProvider.NextNormal();
                        fZ0i[i, j] = (float)Math.Exp(fZ0 + fZi[i, j]);
                    }
                }
                fPR = uMath.DoubleArrayToFloatArray2D(uMath.RandomNormal(nPopulation, nTotalLayers));
                fPT = uMath.DoubleArrayToFloatArray2D(uMath.RandomStudentT2D(nPopulation, nTotalLayers));

                for (int i = 0; i < nPopulation; i = i + 2)
                {
                    for (int j = 0; j < nTotalLayers; j++)
                    {
                        fSigmaC[i, j] = fSigmaC[i, j] * fZ0i[i, j];
                        fPC[i, j] = fPC[i, j] + fPR[i, j] * fSigmaC[i, j];
                    }
                }
                for (int i = 1; i < nPopulation; i = i + 2)
                {
                    for (int j = 0; j < nTotalLayers; j++)
                    {
                        fSigmaC[i, j] = fSigmaC[i, j] * fZ0i[i, j];
                        fPC[i, j] = fPC[i, j] + fPT[i, j] * fSigmaC[i, j];
                    }
                }

                // Constraints:
                for (int i = 0; i < nPopulation; i++)
                {
                    for (int j = 0; j < nTotalLayers; j++)
                    {
                        if (fPC[i, j] < fThicknessFrom)
                        {
                            fPC[i, j] = fThicknessFrom;
                        }
                        if (fPC[i, j] > fThicknessTo)
                        {
                            fPC[i, j] = fThicknessTo;
                        }
                    }
                }

                // Select which vectors are allowed to enter the new population:
                for (int i = 0; i < nPopulation; i++)
                {
                    for (int j = 0; j < nTotalLayers; j++)
                    {
                        if (sOptimization[j] == "No")
                        {
                            fPC[i, j] = fThickness[j];
                        }
                    }
                }

                for (int i = 1; i < nPopulation; i++)
                {
                    string[] sThick = new string[fThickness.Length];
                    for (int j = 0; j < nTotalLayers; j++)
                    {
                        fThickness[j] = fPC[i, j];
                    }

                    aListResult = MeritFunctionResults(nFlag, fLambda, fTheta0, fThickness,
                                           aListTarget, mNK, vNKIncidentMedium, vNKExitMedium);
                    float fTempMeritValue = (float)aListResult[0];
                    if (fTempMeritValue <= fMeritValue[i])
                    {
                        for (int j = 0; j < nTotalLayers; j++)
                        {
                            fPopulation[i, j] = fPC[i, j];
                        }
                        fMeritValue[i] = fTempMeritValue;

                        // Update best merit function only in case of success to save time

                        if (fTempMeritValue < fBestValue)
                        {
                            fBestValue = fTempMeritValue;
                            for (int j = 0; j < nTotalLayers; j++)
                            {
                                fBestPopulationMember[j] = fPC[i, j];
                                if (sOptimization[j] == "Yes")
                                {
                                    fThickness[j] = fPC[i, j];
                                    fThickness0[j] = fPC[i, j];
                                }
                            }
                            aListBestResult = aListResult;
                        }
                    }
                }

                fBestPopulationMemberInIterations = fBestPopulationMember;

                // Pass intermediate results to frmOptimization to plot:
                aListIterations.Add(nIterations);
                aListMeritValue.Add(fBestValue);
                int[] nArrayIterations = (int[])aListIterations.ToArray(typeof(int));
                float[] fArrayMeritValue = (float[])aListMeritValue.ToArray(typeof(float));

                ArrayList aListResults = new ArrayList();
                aListResults.Add(nIterations);
                aListResults.Add(dLambda);
                aListResults.Add(aListTarget);
                aListResults.Add(aListBestResult);
                aListResults.Add(nArrayIterations);
                aListResults.Add(fArrayMeritValue);
                aListResults.Add(nSubstrate);
                aListResults.Add(fThickness0);
                aListResults.Add(fN0);
                aListResults.Add(dLambda0);
                aListResults.Add(nGenerations);
                aListResults.Add(nFlag);

                Thread.Sleep(1);

                frmOptim.Invoke(frmOptim.m_DelegateAddResults, new Object[] { aListResults });

                // check if thread is cancelled
                if (m_EventStop.WaitOne(0, true))
                {
                    m_EventStopped.Set();
                    return;
                }

                // Make asynchronous call to frmOptimization to inform it that thread finished
                //frmOptim.Invoke(frmOptim.m_DelegateThreadFinished, null);            

                nIterations = nIterations + 1;
            } // While loop

            ArrayList aListThick = new ArrayList();
            aListThick.Add(fThickness0);
            hm.HashMainOverWrite("OptimizedThickness", aListThick);

            m_EventStopped.Set();
        }




        public static MatrixR MatrixSort(MatrixR m)
        {
            int rows = m.GetRows();
            int cols = m.GetCols();

            VectorR v = new VectorR(cols);
            MatrixR m1 = new MatrixR(rows,cols);
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
