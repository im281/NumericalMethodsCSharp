using System;

namespace XuMath
{
    public struct VectorC:ICloneable
    {
        public int size;
        public Complex[] vector;

         public VectorC(int size)
        {
            this.size = size;
            this.vector = new Complex[size];
            for (int i = 0; i < size; i++)
            {
                vector[i] = Complex.Zero;
            }
        }

        public VectorC(Complex[] vector)
        {
            this.size = vector.Length;
            this.vector = vector;
        }

        #region Make a deep copy
        public VectorC Clone()
        {
            // returns a deep copy of the vector
            VectorC v = new VectorC(vector);
            v.vector = (Complex[])vector.Clone();
            return v;
        }

        object ICloneable.Clone()
        {
            return Clone();
        }
        #endregion

        #region Equals and Hashing:
        public override bool Equals(object obj)
        {
            return (obj is VectorC) && this.Equals((VectorC)obj);
        }

        public bool Equals(VectorC v)
        {
            return vector == v.vector;
        }

        public override int GetHashCode()
        {
            return vector.GetHashCode();
        }

        public static bool operator ==(VectorC v1, VectorC v2)
        {
            return v1.Equals(v2);
        }

        public static bool operator !=(VectorC v1, VectorC v2)
        {
            return !v1.Equals(v2);
        }
        #endregion

        #region Definition and basics:
        public Complex this[int n]
        {
            get
            {
                if (n < 0 || n > size)
                {
                    throw new ArgumentOutOfRangeException(
                     "n", n, "n is out of range!");
                }
                return vector[n];
            }
            set { vector[n] = value; }
        }

        public double GetNorm()
        {
            Complex result = Complex.Zero;
            for (int i = 0; i < size; i++)
            {
                result += vector[i].Conjugate * vector[i];
            }
            return Math.Sqrt(result.Real);
        }

        public double GetNormSquare()
        {
            Complex result = Complex.Zero;
            for (int i = 0; i < size; i++)
            {
                result += vector[i].Conjugate * vector[i];
            }
            return result.Real;
        }

        public void Normalize()
        {
            double norm = GetNorm();
            if (norm == 0)
            {
                throw new DivideByZeroException("Normalize a vector with norm of zero!");
            }
            for (int i = 0; i < size; i++)
            {
                vector[i] /= norm;
            }
        }

        public VectorC GetUnitVector()
        {
            VectorC result = new VectorC(vector);
            result.Normalize();
            return result;
        }

        public int GetSize()
        {
            return size;
        }

        public VectorC GetConjugate()
        {
            for (int i = 0; i < size; i++)
            {
                vector[i] = vector[i].Conjugate;
            }
            return new VectorC(vector);
        }

        public VectorC GetSwap(int m, int n)
        {
            Complex temp = vector[m];
            vector[m] = vector[n];
            vector[n] = temp;
            return new VectorC(vector);
        }
        #endregion

        #region Mathematical operations:

        public static VectorC operator +(VectorC v)
        {
            return v;
        }

        public static VectorC operator -(VectorC v)
        {
            Complex[] result = new Complex[v.GetSize()];
            for (int i = 0; i < v.GetSize(); i++)
            {
                result[i] = -v[i];
            }
            return new VectorC(result);
        }

        public static VectorC operator +(VectorC v1, VectorC v2)
        {
            VectorC result = new VectorC(v1.size);
            for (int i = 0; i < v1.size; i++)
            {
                result[i] = v1[i] + v2[i];
            }
            return result;
        }

        public static VectorC operator +(VectorC v, double d)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = v[i] + d;
            }
            return result;
        }

        public static VectorC operator +(double d, VectorC v)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = v[i] + d;
            }
            return result;
        }

        public static VectorC operator +(VectorC v, Complex c)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = v[i] + c;
            }
            return result;
        }

        public static VectorC operator +(Complex c, VectorC v)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = v[i] + c;
            }
            return result;
        }

        public static VectorC operator -(VectorC v1, VectorC v2)
        {
            VectorC result = new VectorC(v1.size);
            for (int i = 0; i < v1.size; i++)
            {
                result[i] = v1[i] - v2[i];
            }
            return result;
        }

        public static VectorC operator -(VectorC v, double d)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = v[i] - d;
            }
            return result;
        }

        public static VectorC operator -(double d, VectorC v)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = d - v[i];
            }
            return result;
        }

        public static VectorC operator -(VectorC v, Complex c)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = v[i] - c;
            }
            return result;
        }

        public static VectorC operator -(Complex c, VectorC v)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = c - v[i];
            }
            return result;
        }

        public static VectorC operator *(VectorC v, double d)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = v[i] * d;
            }
            return result;
        }

        public static VectorC operator *(double d, VectorC v)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = d * v[i];
            }
            return result;
        }

        public static VectorC operator *(VectorC v, Complex c)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = v[i] * c;
            }
            return result;
        }

        public static VectorC operator *(Complex c, VectorC v)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = c * v[i];
            }
            return result;
        }

        public static VectorC operator /(VectorC v, double d)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = v[i] / d;
            }
            return result;
        }

        public static VectorC operator /(double d, VectorC v)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = d / v[i];
            }
            return result;
        }

        public static VectorC operator /(VectorC v, Complex c)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = v[i] / c;
            }
            return result;
        }

        public static VectorC operator /(Complex c, VectorC v)
        {
            VectorC result = new VectorC(v.size);
            for (int i = 0; i < v.size; i++)
            {
                result[i] = c / v[i];
            }
            return result;
        }
        #endregion;

        #region Public methods:
        public static Complex DotProduct(VectorC v1, VectorC v2)
        {
            Complex result = Complex.Zero;
            for (int i = 0; i < v1.size; i++)
            {
                result += v1[i].Conjugate * v2[i];
            }
            return result;
        }

        public static VectorC Product(VectorC v1, VectorC v2)
        {
            VectorC result = new VectorC(v1.size);
            for (int i = 0; i < v1.size; i++)
            {
                result[i] = v1[i] * v2[i];
            }
            return result;
        }

        public static VectorC CrossProduct(VectorC v1, VectorC v2)
        {
            if (v1.size != 3)
            {
                throw new ArgumentOutOfRangeException(
                 "v1", v1, "Vector v1 must be 3 dimensional!");
            }
            VectorC result = new VectorC(3);
            result[0] = v1[1] * v2[2] - v1[2] * v2[1];
            result[1] = v1[2] * v2[0] - v1[0] * v2[2];
            result[2] = v1[0] * v2[1] - v1[1] * v2[0];
            return result;
        }
        #endregion;

        public override string ToString()
        {
            string s = "(";
            for (int i = 0; i < size -1; i++)
            {
                s += vector[i] + ", ";
            }
            s += vector[size - 1] + ")";
            return s;
        }

    }
}
