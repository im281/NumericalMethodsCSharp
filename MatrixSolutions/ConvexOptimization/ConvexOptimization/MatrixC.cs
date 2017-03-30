using System;

namespace ConvexOptimization
{
    public struct MatrixC:ICloneable
    {
        private int Rows;
        private int Cols;
        private double[,] matrix;

        public MatrixC(int Rows, int Cols)
        {
            this.Rows = Rows;
            this.Cols = Cols;
            this.matrix = new double[Rows, Cols];
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Cols; j++)
                {
                    matrix[i, j] = 0.0;
                }
            }
        }

        public MatrixC(double[,] matrix)
        {
            this.Rows = matrix.GetLength(0);
            this.Cols = matrix.GetLength(1);
            this.matrix = matrix;
        }

        public MatrixC(MatrixC m)
        {
            Rows = m.GetRows();
            Cols = m.GetCols();
            matrix = m.matrix;
        }

        #region Make a deep copy
        public MatrixC Clone()
        {
            // returns a deep copy of the matrix
            MatrixC m = new MatrixC(matrix);
            m.matrix = (double[,])matrix.Clone();
            return m;
        }

        object ICloneable.Clone()
        {
            return Clone();
        }     
        #endregion

        #region Indexing:
        public double this[int m, int n]
        {
            get
            {
                if (m < 0 || m > Rows)
                {
                    throw new ArgumentOutOfRangeException(
                     "m", m, "m is out of range!");
                }
                if (n < 0 || n > Cols)
                {
                    throw new ArgumentOutOfRangeException(
                     "n", n, "n is out of range!");
                }
                return matrix[m, n];
            }
            set { matrix[m, n] = value; }
        }
        #endregion

        #region Identity matrix:
        public MatrixC Identity()
        {
            MatrixC m = new MatrixC(Rows, Cols);
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Cols; j++)
                {
                    if (i == j)
                    {
                        m[i, j] = 1;
                    }
                }
            }
            return m;
        }
        #endregion

        #region Equals and Hashing:
        public override bool Equals(object obj)
        {
            return (obj is MatrixC) && this.Equals((MatrixC)obj);
        }

        public bool Equals(MatrixC m)
        {
            return matrix == m.matrix;
        }

        public override int GetHashCode()
        {
            return matrix.GetHashCode();
        }

        public static bool operator ==(MatrixC m1, MatrixC m2)
        {
            return m1.Equals(m2);
        }

        public static bool operator !=(MatrixC m1, MatrixC m2)
        {
            return !m1.Equals(m2);
        }
        #endregion

        #region Matrix dimension
        public int GetRows()
        {
            return Rows;
        }

        public int GetCols()
        {
            return Cols;
        }

        public bool IsSquared()
        {
            if (Rows == Cols)
                return true;
            else
                return false;
        }

        public static bool CompareDimension(MatrixC m1, MatrixC m2)
        {
            if (m1.GetRows() == m2.GetRows() && m1.GetCols() == m2.GetCols())
                return true;
            else
                return false;
        }
        #endregion

        #region Get row or column vector from a matrix:
        public VectorC GetRowVector(int n)
        {
            if (n < 0 || n > Rows)
            {
                throw new ArgumentOutOfRangeException(
                 "n", n, "n is out of range!");
            }
            VectorC v = new VectorC(Cols);
            for (int i = 0; i < Cols; i++)
            {
                v[i] = matrix[n, i];
            }
            return v;
        }

        public VectorC GetColVector(int n)
        {
            if (n < 0 || n > Cols)
            {
                throw new ArgumentOutOfRangeException(
                 "n", n, "n is out of range!");
            }
            VectorC v = new VectorC(Rows);
            for (int i = 0; i < Rows; i++)
            {
                v[i] = matrix[i, n];
            }
            return v;
        }
        #endregion

        #region Replace row or column with vector
        public MatrixC ReplaceRow(VectorC v, int n)
        {
            if (n < 0 || n > Rows)
            {
                throw new ArgumentOutOfRangeException(
                 "n", n, "n is out of range!");
            }
            if (v.GetSize() != Cols)
            {
                throw new ArgumentOutOfRangeException(
                    "Vector size", v.GetSize(), "vector size is out of range!");
            }
            for (int i = 0; i < Cols; i++)
            {
                matrix[n, i] = v[i];
            }
            return new MatrixC(matrix);
        }

        public MatrixC ReplaceCol(VectorC v, int n)
        {
            if (n < 0 || n > Cols)
            {
                throw new ArgumentOutOfRangeException(
                 "n", n, "n is out of range!");
            }
            if (v.GetSize() != Rows)
            {
                throw new ArgumentOutOfRangeException(
                    "Vector size", v.GetSize(), "vector size is out of range!");
            }
            for (int i = 0; i < Rows; i++)
            {
                matrix[i, n] = v[i];
            }
            return new MatrixC(matrix);
        }
        #endregion

        #region Matrix transpose and trace:
        public MatrixC GetTranspose()
        {
            MatrixC m = this;
            m.Transpose();
            return m;
        }

        public MatrixC Transpose()
        {
            MatrixC m = new MatrixC(Cols, Rows);
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Cols; j++)
                {
                    m[j, i] = matrix[i, j];
                }
            }
            return m;
        }

        public double GetTrace()
        {
            double d = 0.0;
            for (int i = 0; i < Rows; i++)
            {
                if (i < Cols)
                    d += matrix[i, i];
            }
            return d;
        }
        #endregion

        #region Matrix transformation:
        public static VectorC Transform(MatrixC m, VectorC v)
        {
            //temporary_hack
            VectorC result = new VectorC(m.GetRows());
            //VectorC result = new VectorC(v.GetSize());
            //temporary_hack
            //if (!m.IsSquared())
            //{
            //    throw new ArgumentOutOfRangeException(
            //     "Dimension", m.GetRows(), "The matrix must be squared!");
            //}
            //if (m.GetCols() != v.GetSize())
            //{
            //    throw new ArgumentOutOfRangeException(
            //     "Size", v.GetSize(), "The size of the vector must be equal"
            //     + "to the number of rows of the matrix!");
            //}
            for (int i = 0; i < m.GetRows(); i++)
            {
                result[i] = 0.0;
                for (int j = 0; j < m.GetCols(); j++)
                {
                    result[i] += m[i, j] * v[j];
                }
            }
            return result;
        }

        public static VectorC Transform(VectorC v, MatrixC m)
        {
            VectorC result = new VectorC(v.GetSize());
            //if (!m.IsSquared())
            //{
            //    throw new ArgumentOutOfRangeException(
            //     "Dimension", m.GetRows(), "The matrix must be squared!");
            //}
            //if (m.GetRows() != v.GetSize())
            //{
            //    throw new ArgumentOutOfRangeException(
            //     "Size", v.GetSize(), "The size of the vector must be equal"
            //     + "to the number of rows of the matrix!");
            //}
            for (int i = 0; i < m.GetRows(); i++)
            {
                result[i] = 0.0;
                for (int j = 0; j < m.GetCols(); j++)
                {
                    result[i] += v[j] * m[j, i];
                }
            }
            return result;
        }

        public static MatrixC Transform(VectorC v1, VectorC v2)
        {
            if (v1.GetSize() != v2.GetSize())
            {
                throw new ArgumentOutOfRangeException(
                 "v1", v1.GetSize(), "The vectors must have the same size!");
            }
            MatrixC result = new MatrixC(v1.GetSize(), v1.GetSize());
            for (int i = 0; i < v1.GetSize(); i++)
            {
                for (int j = 0; j < v1.GetSize(); j++)
                {
                    result[j, i] = v1[i] * v2[j];
                }
            }
            return result;
        }
        #endregion

        #region Interchange rows or columns
        public MatrixC GetRowSwap(int m, int n)
        {
            double temp = 0.0;
            for (int i = 0; i < Cols; i++)
            {
                temp = matrix[m, i];
                matrix[m, i] = matrix[n, i];
                matrix[n, i] = temp;
            }
            return new MatrixC(matrix);
        }

        public MatrixC GetColSwap(int m, int n)
        {
            double temp = 0.0;
            for (int i = 0; i < Rows; i++)
            {
                temp = matrix[i, m];
                matrix[i, m] = matrix[i, n];
                matrix[i, n] = temp;
            }
            return new MatrixC(matrix);
        }
        #endregion

        #region Mathematical operators
        public static MatrixC operator +(MatrixC m)
        {
            return m;
        }

        public static MatrixC operator -(MatrixC m)
        {
            for (int i = 0; i < m.GetRows(); i++)
            {
                for (int j = 0; j < m.GetCols(); j++)
                {
                    m[i, j] = -m[i, j];
                }
            }
            return m;
        }

        public static MatrixC operator +(MatrixC m1, MatrixC m2)
        {
            if (!MatrixC.CompareDimension(m1, m2))
            {
                throw new ArgumentOutOfRangeException(
                 "Dimension", m1, "The dimensions of two matrices must be the same!");
            }
            MatrixC result = new MatrixC(m1.GetRows(), m1.GetCols());
            for (int i = 0; i < m1.GetRows(); i++)
            {
                for (int j = 0; j < m1.GetCols(); j++)
                {
                    result[i, j] = m1[i, j] + m2[i, j];
                }
            }
            return result;
        }

        public static MatrixC operator +(MatrixC m, double d)
        {
            MatrixC result = new MatrixC(m.GetRows(), m.GetCols());
            for (int i = 0; i < m.GetRows(); i++)
            {
                for (int j = 0; j < m.GetCols(); j++)
                {
                    result[i, j] = m[i, j] + d;
                }
            }
            return result;
        }

        public static MatrixC operator +(double d, MatrixC m)
        {
            MatrixC result = new MatrixC(m.GetRows(), m.GetCols());
            for (int i = 0; i < m.GetRows(); i++)
            {
                for (int j = 0; j < m.GetCols(); j++)
                {
                    result[i, j] = m[i, j] + d;
                }
            }
            return result;
        }

        public static MatrixC operator -(MatrixC m1, MatrixC m2)
        {
            if (!MatrixC.CompareDimension(m1, m2))
            {
                throw new ArgumentOutOfRangeException(
                 "Dimension", m1, "The dimensions of two matrices must be the same!");
            }
            MatrixC result = new MatrixC(m1.GetRows(), m1.GetCols());
            for (int i = 0; i < m1.GetRows(); i++)
            {
                for (int j = 0; j < m1.GetCols(); j++)
                {
                    result[i, j] = m1[i, j] - m2[i, j];
                }
            }
            return result;
        }

        public static MatrixC operator -(MatrixC m, double d)
        {
            MatrixC result = new MatrixC(m.GetRows(), m.GetCols());
            for (int i = 0; i < m.GetRows(); i++)
            {
                for (int j = 0; j < m.GetCols(); j++)
                {
                    result[i, j] = m[i, j] - d;
                }
            }
            return result;
        }

        public static MatrixC operator -(double d, MatrixC m)
        {
            MatrixC result = new MatrixC(m.GetRows(), m.GetCols());
            for (int i = 0; i < m.GetRows(); i++)
            {
                for (int j = 0; j < m.GetCols(); j++)
                {
                    result[i, j] = d - m[i, j];
                }
            }
            return result;
        }

        public static MatrixC operator *(MatrixC m1, MatrixC m2)
        {
            if (m1.GetCols()!=m2.GetRows())
            {
                throw new ArgumentOutOfRangeException(
                 "Columns", m1, "The numbers of columns of the first matrix must be equal to" +
                 " the number of rows of the second matrix!");
            }
            MatrixC result = new MatrixC(m1.GetRows(), m2.GetCols());
            VectorC v1 = new VectorC(m1.GetCols());
            VectorC v2 = new VectorC(m2.GetRows());
            for (int i = 0; i < m1.GetRows(); i++)
            {
                v1 = m1.GetRowVector(i);
                for (int j = 0; j < m2.GetCols(); j++)
                {
                    v2 = m2.GetColVector(j);
                    result[i, j] = VectorC.DotProduct(v1, v2);
                }
            }
            return result;
        }

        public static MatrixC operator *(MatrixC m, double d)
        {
            MatrixC result = new MatrixC(m.GetRows(), m.GetCols());
            for (int i = 0; i < m.GetRows(); i++)
            {
                for (int j = 0; j < m.GetCols(); j++)
                {
                    result[i, j] = m[i, j] * d;
                }
            }
            return result;
        }

        public static MatrixC operator *(double d, MatrixC m)
        {
            MatrixC result = new MatrixC(m.GetRows(), m.GetCols());
            for (int i = 0; i < m.GetRows(); i++)
            {
                for (int j = 0; j < m.GetCols(); j++)
                {
                    result[i, j] = m[i, j] * d;
                }
            }
            return result;
        }

        public static MatrixC operator /(MatrixC m, double d)
        {
            MatrixC result = new MatrixC(m.GetRows(), m.GetCols());
            for (int i = 0; i < m.GetRows(); i++)
            {
                for (int j = 0; j < m.GetCols(); j++)
                {
                    result[i, j] = m[i, j] / d;
                }
            }
            return result;
        }

        public static MatrixC operator /(double d, MatrixC m)
        {
            MatrixC result = new MatrixC(m.GetRows(), m.GetCols());
            for (int i = 0; i < m.GetRows(); i++)
            {
                for (int j = 0; j < m.GetCols(); j++)
                {
                    result[i, j] = d/ m[i, j];
                }
            }
            return result;
        }
        #endregion

        #region Determinant of a matrix
        public static double Determinant(MatrixC m)
        {
            double result = 0.0;
            if (!m.IsSquared())
            {
                throw new ArgumentOutOfRangeException(
                 "Dimension", m.GetRows(), "The matrix must be squared!");
            }
            if (m.GetRows() == 1)
                result = m[0, 0];
            else
            {
                for (int i = 0; i < m.GetRows(); i++)
                {
                    result += Math.Pow(-1, i)*m[0, i] * Determinant(MatrixC.Minor(m, 0, i));
                }
            }
            return result;
        }

        public static MatrixC Minor(MatrixC m, int row, int col)
        {
            MatrixC mm = new MatrixC(m.GetRows() - 1, m.GetCols() - 1);
            int ii = 0, jj = 0;
            for (int i = 0; i < m.GetRows(); i++)
            {
                if (i == row)
                    continue;
                jj = 0;
                for (int j = 0; j < m.GetCols(); j++)
                {
                    if (j == col)
                        continue;
                    mm[ii, jj] = m[i, j];
                    jj++;
                }
                ii++;
            }
            return mm;
        }

        public static MatrixC Inverse(MatrixC m)
        {
            if (Determinant(m) == 0)
            {
                throw new DivideByZeroException(
                    "Cannot inverse a matrix with a zero determinant!");
            }
            return (Adjoint(m) / Determinant(m));
        }

        public static MatrixC Adjoint(MatrixC m)
        {
            if (!m.IsSquared())
            {
                throw new ArgumentOutOfRangeException(
                 "Dimension", m.GetRows(), "The matrix must be squared!");
            }
            MatrixC ma = new MatrixC(m.GetRows(), m.GetCols());
            for (int i = 0; i < m.GetRows(); i++)
            {
                for (int j = 0; j < m.GetCols(); j++)
                {
                    ma[i, j] = Math.Pow(-1, i + j) * (Determinant(Minor(m, i, j)));
                }
            }
            return ma.GetTranspose();
        }

        #endregion

        #region Overrise ToString method
        public override string ToString()
        {
            string ss = "(";
            for (int i = 0; i < Rows; i++)
            {
                string s = "";
                for (int j = 0; j < Cols - 1; j++)
                {
                    s += matrix[i, j].ToString() + ", ";
                }
                s += matrix[i, Cols - 1].ToString();
                if (i != Rows - 1 && i == 0)
                    ss += s + "\n";
                else if (i != Rows - 1 && i != 0)
                    ss += " " + s + "\n";
                else
                    ss += " " + s + ")";
            }
            return ss;
        }
        #endregion
    }
}
