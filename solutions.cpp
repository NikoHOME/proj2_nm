#include "matrix.hpp"
#include "solutions.hpp"

IterativeSolution::IterativeSolution(Matrix *matrix, uint32_t max_iterations) {
    this->matrix = matrix;
    this->norm_history = new float64_t[max_iterations];
}

IterativeSolution::~IterativeSolution() {
    delete this->matrix;
    delete[] this->norm_history;
}


DirectSolution::DirectSolution(Matrix *matrix) {
    this->matrix = matrix;
}

DirectSolution::~DirectSolution() {
    delete this->matrix;
}