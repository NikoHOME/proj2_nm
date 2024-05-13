

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstdint>
#include "types.hpp"
#include "solutions.hpp"

enum MatrixType { LOWER_TRIANGLE, UPPER_TRIANGLE, IDENTITY, DIAGONAL };


class MatrixReadOnly {
    public:
        uint32_t height;
        uint32_t width;
        float64_t *content;

        void print();
        void print(uint32_t x_limit, uint32_t y_limit);
        float64_t get(uint32_t x, uint32_t y);

        Matrix *solution(MatrixReadOnly *matrix);
        Matrix *solution_lu(MatrixReadOnly *matrix);
        Matrix *clone();
        float64_t norm();
        DirectSolution *solution_direct(MatrixReadOnly *matrix);
        IterativeSolution *solution_jacobi(MatrixReadOnly *input, uint32_t iter, float64_t err_norm);
        IterativeSolution *solution_gauss_seidel(MatrixReadOnly *input, uint32_t iter, float64_t err_norm);
        Matrix *gauss_matrix(MatrixReadOnly *matrix);
        void mutiply_into(MatrixReadOnly *matrix, Matrix *dest);
        
};


class MatrixShapeLookup : public MatrixReadOnly {
    protected:
        MatrixType type;
    public:
        MatrixShapeLookup(MatrixReadOnly *matrix, MatrixType type);
        float64_t get(uint32_t x, uint32_t y);
        void print();
        void print(uint32_t x_limit, uint32_t y_limit);
};


class Matrix : public MatrixReadOnly {
    public:
        Matrix(uint32_t width, uint32_t height);
        Matrix(uint32_t width, uint32_t height, float64_t *content);
        Matrix(uint32_t width, uint32_t height, float64_t val);
        Matrix(MatrixReadOnly *matrix);
        ~Matrix();

        void set(uint32_t x, uint32_t y, float64_t val);
        void add(MatrixReadOnly *matrix);
        void sub(MatrixReadOnly *matrix);
        void copy(MatrixReadOnly *matrix);
        void transform(MatrixType type);
        void mutiply(MatrixReadOnly *matrix);
        
        void dot(float64_t val);
        void swap_rows(uint32_t y_row1, uint32_t y_row2);
        
        void invert();
        void transpose();
        bool is_invertible();
        void divide(MatrixReadOnly *matrix);
        
};



#endif // MATRIX_HPP
