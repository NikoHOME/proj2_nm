#ifndef SOLUTIONS_HPP
#define SOLUTIONS_HPP

#include "types.hpp"

class Matrix;

class IterativeSolution {
    public:
        Matrix *matrix;
        float64_t *norm_history;
        uint32_t max_iterations;
        uint32_t iterations;
        uint32_t execution_time;

        IterativeSolution(Matrix *matrix, uint32_t max_iterations);
        ~IterativeSolution();
};

class DirectSolution {
    public:
        Matrix *matrix;
        uint32_t execution_time;
        DirectSolution(Matrix *matrix);
        ~DirectSolution();
};

#endif // SOLUTIONS_HPP