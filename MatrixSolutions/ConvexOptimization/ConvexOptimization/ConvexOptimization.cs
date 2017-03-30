using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace ConvexOptimization
{
    public class ConvexOptimization
    {

        /// <summary>
        /// Implementation Convex optimization using ADMM( Stephen Boyd, Stanford) It is essentially least squares fitting with non-negativity constraints
        /// An ANOVA computation is added below for significance testing although this implementation is NOT VALIDATED
        /// The matrix is any reference matrix (database) and y is the experimental data. The result is a vector of coefficients for the model
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public VectorC MatchLibraryMatrix(MatrixC matrix, VectorC y)
        {
       
            #region ADMM (Alternating direction method of multipliers)

            //regularized matrix version of least squares
            //x = (A'*A)*A'*b

            //rho, the wiggle factor
            double rho = .001;

            //Ax = b where x is the prediction, b is the experimental and A is the reference database matrix  

            // A = x^t*x
            var A = matrix.Transpose() * matrix;
        
            //b = xt*y
            MatrixC y1 = new MatrixC(matrix.GetRows(), 1);
            y1.ReplaceCol(y, 0);
            
            //var b = MatrixC.Transform(matrix.Transpose(), y); 
            var temp1 = matrix.Transpose() * y1;
            VectorC b = new VectorC(temp1.GetRows());
            for(int i = 0; i < temp1.GetRows();i++)
            {
                b[i] = temp1[i,0];
            }
            //create identity matrix
            var I = A.Identity();

            //A' = (A + rho*I)
            var Aprime = MatrixC.Inverse(A + rho * I);

            double[] components = new double[b.GetSize()];
            //add zeros to vectors
            for (int i = 0; i < b.GetSize(); i++)
            {
                components[i] = 0;
            }

            VectorC u = new VectorC(components); //accululated error in the negative direction
            VectorC z = new VectorC(components); //non-negative version of x
            VectorC x = new VectorC(components); //regularized least squares
           
            u = u.Clone();
            z = z.Clone();
            x = x.Clone();

            var test = b + rho * (z - u);
            double sumSquares = 100;

            for (int i = 0; i < 1500; i++)
            {
                x = MatrixC.Transform(Aprime, (b + rho * (z - u)));
         
                //for each element in x + u
                //z = (x + u);
                for (int j = 0; j < x.GetSize(); j++)
                {
                    var temp = x[j] + u[j];
                    if (temp > 0)
                    {
                        z[j] = temp;
                    }
                    else
                    {
                        z[j] = 0;
                    }

                }
                //for each element in x - z
                // u = u + (x - z);
                for (int k = 0; k < u.GetSize(); k++)
                {
                    u[k] = u[k] + (x[k] - z[k]);
                }

                //calculate RMS error
                for (int ii = 0; ii < x.GetSize(); ii++)
                {
                    sumSquares += (x[ii] - z[ii])*(x[ii] - z[ii]);
                }

                var RMS = Math.Sqrt(sumSquares / x.GetSize());
                if (RMS < .001)
                {
                    break;
                }


            }

            #endregion

            #region ANOVA

            //Matrix yMatrix = Matrix.Create(y);
            //Matrix database = (m.Transpose()*m).GetInverse()*m.Transpose();

            ////ADMM
            //var predicted = m*z;

            ////Hat matrix. The hat matrix is used to reduce overfitting by finding outliers 
            ////inv(x'*x)*x'*x'
            //var hatMatrix = (m*(m.Transpose()*m).GetInverse())*m.Transpose();

            ////J matrix square matrix of ones
            //Matrix JMatrix = Matrix.Create(hatMatrix.ColumnCount, hatMatrix.RowCount);
            //for (int i = 0; i < JMatrix.ColumnCount; i++)
            //{
            //    for (int j = 0; j < JMatrix.RowCount; j++)
            //    {
            //        JMatrix[i, j] = 1;
            //    }
            //}

            ////SSr  (2000 is the observations or m/z in our case)
            ////y'*[H - (1/n)J]*y
            //var SSr = Y.Transpose()*((hatMatrix - (1/2000)*JMatrix))*Y;

            ////Mean Square Regression
            //var MSR = SSr/matrix.GetLongLength(1);

            ////Identity matrix
            //// The following constructs a nxn identity matrix:
            //DenseMatrix identityMatrix = DenseMatrix.GetIdentity(2001);

            ////SSE
            ////y'*[I-H]*y                          
            //var SSe = Y.Transpose()*(identityMatrix - hatMatrix)*Y;

            //double MSquareRegression = MSR[0, 0];

            ////Mean Square Error
            //var MSE = SSe/(2000 - (matrix.GetLongLength(1) + 1));

            //double MSquareError = MSE[0, 0];

            ////F-statistic
            ////The fstat conducts a hypothesis test. It simply says, is there at least one coefficient that 
            ////id different than zero?
            //var f = MSquareRegression/MSquareError;



            ////get covariance matrix
            //var C = MSquareError*(m.Transpose()*m).GetInverse();

            ////get diangonal
            //var d = C.GetDiagonal();

            ////conduct t-test for all coefficients (vectorResult for OLS)
            //for (int i = 0; i < z.Length; i++)
            //{
            //    double tTest = z[i]/Math.Sqrt(d[i]);

            //}

            ////populate correlation scores (vectorResult for OLS)
            //for (int i = 0; i < z.Length; i++)
            //{

            //}

            

            #endregion

            return z;

        }

        public VectorC MatchLibraryMatrixLS(MatrixC A, VectorC b)
        {    
            //regularized matrix version of least squares Ax = b solution 
            //x = Inverse((A'*A))*A'*b
            return MatrixC.Transform(MatrixC.Inverse((A.Transpose() * A)) * ( A.Transpose()), b);    
       
        }
        public void ANOVA()
        {
            #region ANOVA

            //Matrix yMatrix = Matrix.Create(y);
            //Matrix database = (m.Transpose()*m).GetInverse()*m.Transpose();

            ////ADMM
            //var predicted = m*z;

            ////Hat matrix. The hat matrix is used to reduce overfitting by finding outliers 
            ////inv(x'*x)*x'*x'
            //var hatMatrix = (m*(m.Transpose()*m).GetInverse())*m.Transpose();

            ////J matrix square matrix of ones
            //Matrix JMatrix = Matrix.Create(hatMatrix.ColumnCount, hatMatrix.RowCount);
            //for (int i = 0; i < JMatrix.ColumnCount; i++)
            //{
            //    for (int j = 0; j < JMatrix.RowCount; j++)
            //    {
            //        JMatrix[i, j] = 1;
            //    }
            //}

            ////SSr  (2000 is the observations or m/z in our case)
            ////y'*[H - (1/n)J]*y
            //var SSr = Y.Transpose()*((hatMatrix - (1/2000)*JMatrix))*Y;

            ////Mean Square Regression
            //var MSR = SSr/matrix.GetLongLength(1);

            ////Identity matrix
            //// The following constructs a nxn identity matrix:
            //DenseMatrix identityMatrix = DenseMatrix.GetIdentity(2001);

            ////SSE
            ////y'*[I-H]*y                          
            //var SSe = Y.Transpose()*(identityMatrix - hatMatrix)*Y;

            //double MSquareRegression = MSR[0, 0];

            ////Mean Square Error
            //var MSE = SSe/(2000 - (matrix.GetLongLength(1) + 1));

            //double MSquareError = MSE[0, 0];

            ////F-statistic
            ////The fstat conducts a hypothesis test. It simply says, is there at least one coefficient that 
            ////id different than zero?
            //var f = MSquareRegression/MSquareError;



            ////get covariance matrix
            //var C = MSquareError*(m.Transpose()*m).GetInverse();

            ////get diangonal
            //var d = C.GetDiagonal();

            ////conduct t-test for all coefficients (vectorResult for OLS)
            //for (int i = 0; i < z.Length; i++)
            //{
            //    double tTest = z[i]/Math.Sqrt(d[i]);

            //}

            ////populate correlation scores (vectorResult for OLS)
            //for (int i = 0; i < z.Length; i++)
            //{

            //}



            #endregion
        }
    }
}

