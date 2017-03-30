using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConvexOptimization
{
    public class IsotopeClusterDeconvolution
    {
        ConvexOptimization solver = new ConvexOptimization();

        static readonly double[,] intensities = new[,] {
	{0000, 1.00, 0.00, 0.00, 0.00, 0.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00, 0.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{1000, 1.00, 0.53, 0.17, 0.04, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00, 0.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{1500, 1.00, 0.82, 0.42, 0.16, 0.05, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00, 0.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{2000, 0.92, 1.00, 0.63, 0.29, 0.11, 0.03, 0.01, 0.0, 0.0, 0.0, 0.0, 0.00, 0.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{2500, 0.74, 1.00, 0.76, 0.41, 0.18, 0.07, 0.02, 0.01, 0.0, 0.0, 0.0, 0.00, 0.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{3000, 0.62, 1.00, 0.88, 0.56, 0.28, 0.12, 0.04, 0.01, 0.0, 0.0, 0.0, 0.00, 0.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{4000, 0.39, 0.85, 1.00, 0.83, 0.55, 0.30, 0.14, 0.06, 0.02, 0.01, 0.0, 0.00, 0.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{5000, 0.26, 0.70, 1.00, 1.00, 0.78, 0.51, 0.28, 0.26, 0.14, 0.06, 0.02, 0.00, 0.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{6000, 0.15, 0.50, 0.85, 1.00, 0.92, 0.69, 0.44, 0.25, 0.27, 0.13, 0.06, 0.03, 0.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{7000, 0.09, 0.35, 0.69, 0.94, 1.00, 0.87, 0.66, 0.43, 0.25, 0.13, 0.07, 0.03, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{8000, 0.06, 0.24, 0.54, 0.84, 1.00, 0.98, 0.82, 0.60, 0.39, 0.23, 0.13, 0.06, 0.03, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0},
	{9000,0.033, 0.16, 0.41, 0.70, 0.93, 1.00, 0.92, 0.74, 0.53, 0.35, 0.21, 0.11, 0.06, 0.03, 0.01, 0.0, 0.0, 0.0, 0.0},
	{10000,.019, 0.10, 0.29, 0.55, 0.81, 0.98, 1.00, 0.89, 0.71, 0.51, 0.34, 0.20, 0.12, 0.06, 0.03, 0.0, 0.0, 0.0, 0.0},
	{15000,.010, 0.04, 0.12, 0.26, 0.45, 0.66, 0.85, 0.97, 1.00, 0.94, 0.81, 0.65, 0.48, 0.34, 0.22, 0.14, 0.08, 0.05, 0.03},
    {50000,.0001,0.0002, 0.0004, 0.008, 0.016, 0.032, 0.06, 0.12, 0.24, 0.48, 1.00, 0.47, 0.30, 0.19, 0.10, 0.05, 0.03, 0.0, 0.0}};

        /// <summary>
        /// GetModelIntensities - for a given mass, return the isotope pattern as a list of 19 doubles,
        /// representing relative intensities. This is done via a look-up table and linear interpolation,
        /// so the values are only approximate, but useful to find the A0 m/z in noisy data
        /// </summary>
        /// <param name="mass">The deconvolved mass</param>
        /// <returns></returns>
        private static List<double> GetModelIntensities(double mass)
        {
            List<double> result = new List<double>();

            // we don't handle huge masses
            if (mass > 0 && mass < 50000)
            {
                int index = 0;
                while (intensities[index, 0] < mass)
                {
                    index++;
                }
                // since there was an entry for mass 0, at this point index is the higher mass
                // and index-1 is the lower mass
                double offsetIntoRange = mass - intensities[index - 1, 0];
                double positionInRange = offsetIntoRange / (intensities[index, 0] - intensities[index - 1, 0]);

                // create the linear interpolation for each value to be returned
                int width = intensities.GetLength(1);
                for (int j = 1; j < width; j++)
                {
                    double lowValue = intensities[index - 1, j];
                    double highValue = intensities[index, j];
                    double deltaValue = highValue - lowValue;
                    double fractionOfDelta = deltaValue * positionInRange;
                    double actualValue = lowValue + fractionOfDelta;
                    result.Add(actualValue);
                }
            }
            return result;
        }

        private List<double> NormalizeIntensities(List<double> rawIntensities)
        {
            var ratio = 1.0 / rawIntensities.Max();
            List<double> normalizedIntensities = rawIntensities.Select(i => i * ratio).ToList();
            return normalizedIntensities;
        }

        public List<double> NonNegLeastSquares(double[,] matrix, List<double> values)
        {

            List<double> normalizedValues = this.NormalizeIntensities(values);

            //to handle < 20 intensiites for isotopes
            //otherwise vector-matrix length mismatch and crash
            List<double> tempValues = new List<double>();
            for (int i = 0; i < 20; i++)
            {
                if (i < normalizedValues.Count)
                {
                    tempValues.Add(normalizedValues[i]);
                }
                else
                {
                    tempValues.Add(0);
                }
            }

            List<double> results = new List<double>();
            double[,] tempMatrix = new double[tempValues.Count, 1];

            //hardcoded to handle only 20 isotope intensities
            for (int i = 0; i < 20; i++)
            {
                if (i < normalizedValues.Count)
                {
                    tempMatrix[i, 0] = normalizedValues[i];
                }
                else
                {
                    tempMatrix[i, 0] = 0;
                }
            }

            MatrixC theoretical = new MatrixC(tempMatrix);

            double[] experimental = new double[values.Count];

            for (int i = 0; i < values.Count; i++)
            {
                experimental[i] = values[i];
            }

            VectorC v = new VectorC(experimental);

            var result = solver.MatchLibraryMatrix(theoretical, v);


            for (int i = 0; i < result.GetSize(); i++)
            {
                results.Add(result[i]);
            }

            //calcualte %of each component
            double sum = 0;
            for (int j = 0; j < result.GetSize(); j++)
            {
                sum += results[j];
            }

            for (int k = 0; k < result.GetSize(); k++)
            {
                results[k] = results[k] / sum;
            }
            return results;

        }
    }
}
