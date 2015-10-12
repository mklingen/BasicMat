// Copyright (c) 2015 Matt Klingensmith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef BASICMAT_H_
#define BASICMAT_H_

#include <string>
#include <sstream>

// An extremely simple, dumb single-header library for basic 2D floating point matricies of fixed size.
// Does no fancy range checking, asserts, or any other safety features.
namespace basic_mat
{
  // Matrix has N rows and M columns.
  template <size_t N, size_t M> class BasicMat
  {
        public:
            typedef BasicMat<N, M> Type;
            BasicMat()
            {
                   for(size_t i = 0; i < N * M; i++)
                   {
                       m[i] = 0.0f;
                   }
            }

            ~BasicMat()
            {

            }
            
            // Print the matrix to a string.
            inline std::string ToString()
            {
                std::stringstream ss;

                for (size_t r = 0; r < N; r++)
                {
                    for (size_t c = 0; c < M; c++)
                    {
                        ss << (*this)(r, c) << " ";
                    }

                    ss << "\n";
                }

                return ss.str();
            }

            inline size_t GetNumRows()
            {
                return N;
            }

            inline size_t GetNumCols()
            {
                return M;
            }
            // Get floating point value at linear index, row-major order.
            float& operator[](size_t idx)
            {
                return m[idx];
            }

            const float& operator[](size_t idx) const
            {
                return m[idx];
            }

            // Get the floating point value at a row and column.
            float& operator()(size_t r, size_t c = 0)
            {
                return (*this)[r * M + c];
            }

            const float& operator()(size_t r, size_t c  = 0) const
            {
                return (*this)[r * M + c];
            }

            void operator=(const Type& other)
            {
                for(size_t i = 0; i < N * M; i++)
                 {
                     m[i] = other[i];
                 }
            }

            // Matrix addition
            void operator+=(const Type& other)
            {
                for(size_t i = 0; i < N * M; i++)
                {
                    m[i] += other[i];
                }
            }

            friend Type operator+(Type lhs, const Type& rhs)
            {
                Type toReturn = lhs;
                toReturn += rhs;
                return toReturn;
            }

            // Scalar multiplication
            Type operator*=(const float scalar)
            {
                Type toReturn = *this;
                for(size_t i = 0; i < N * M; i++)
                {
                    toReturn[i] *= scalar;
                }
                return toReturn;
            }

            friend Type operator*(Type lhs, const float& rhs)
            {
                return lhs *= rhs;
            }
            
            friend Type operator*(const float& rhs, Type lhs)
            {
                return lhs *= rhs;
            }

            // Copies the matrix into a transposed version of itself.
            inline BasicMat<M, N> Transpose() const
            {
                BasicMat<M, N> toReturn;

                for(size_t r = 0; r < N; r++)
                {
                    for(size_t c = 0; c < M; c++)
                    {
                        toReturn(c, r) = (*this)(r, c);
                    }
                }
                return toReturn;
            }

            // Premultiply a matrix by another matrix.
            template <size_t N2> BasicMat<N2, M> PreMult(const BasicMat<N2, N>& lhs)
            {
                    BasicMat<N2, M> toReturn;

                    for (size_t r = 0; r < N2; r++)
                    {
                        for(size_t c = 0; c < M; c++)
                        {
                            for(size_t k = 0; k < N; k++)
                            {
                                toReturn(r, c) += lhs(r, k) * (*this)(k, r);
                            }
                        }
                    }
                    return toReturn;
            }

            // Post-multiply a matrix by another matrix.
            template <size_t M2> BasicMat<N, M2> PostMult(const BasicMat<M, M2>& rhs)
            {
                BasicMat<N, M2> toReturn;

                for (size_t r = 0; r < N; r++)
                {
                    for(size_t c = 0; c < M2; c++)
                    {
                        for(size_t k = 0; k < M; k++)
                        {
                            toReturn(r, c) += (*this)(r, k) * rhs(k, c);
                        }
                    }
                }
                return toReturn;
            }

            // Matrix multiplication overload.
            template <size_t M2> BasicMat<N, M2> operator*=(const BasicMat<M, M2>& rhs)
            {
                return PostMult(rhs);
            }

            // No need to implement the other way around. It will be taken care of in the other matrix
            // header.
            template <size_t M2> friend BasicMat<N, M2> operator*(Type lhs, const BasicMat<M, M2>& rhs)
            {
                return lhs *= rhs;
            }

            // The data is in a fixed-size array. Row-major order.
            float m[N * M];
    };
}
#endif // BASICMAT_H_ 

