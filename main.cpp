
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>

#include "types.hpp"
#include "matrix.hpp"
#include "solutions.hpp"
#include <string> 


void make_three_arg_diagonal(Matrix *matrix, float64_t a1, float64_t a2, float64_t a3) {
    for (int32_t diag = 0; diag < matrix->width; diag += 1) {
        matrix->set(diag, diag, a1);

        if (diag - 1 >= 0) {
            matrix->set(diag - 1, diag, a2);
        }
        if (diag - 2 >= 0) {
            matrix->set(diag - 2, diag, a3);
        }
        if (diag + 1 < matrix->width) {
            matrix->set(diag + 1, diag, a2);
        }
        if (diag + 2 < matrix->width) {
            matrix->set(diag + 2, diag, a3);
        }
    }
}


#define MAX_ITERATIONS 1000
#define HIGH_ACCURACY 10e-12
#define LOW_ACCURACY 10e-12



class Solutions {
    public:
        DirectSolution **direct;
        IterativeSolution **jacobi;
        IterativeSolution **gauss_seidel;
        uint32_t size;
        Solutions(uint32_t size) {
            this->size = size;
            this->direct = new DirectSolution*[size];
            this->jacobi = new IterativeSolution*[size];
            this->gauss_seidel = new IterativeSolution*[size];
        }
        ~Solutions() {
            for (uint32_t i = 0; i < this->size; i += 1) {
                delete this->direct[i];
                delete this->jacobi[i];
                delete this->gauss_seidel[i];
            }
            delete[] this->direct;
            delete[] this->jacobi;
            delete[] this->gauss_seidel;
        }
        void solve(uint32_t index, MatrixReadOnly *matrix, MatrixReadOnly *vector) {
            direct[index] = matrix->solution_direct(vector);
            jacobi[index] = matrix->solution_jacobi(vector, MAX_ITERATIONS, HIGH_ACCURACY);
            gauss_seidel[index] = matrix->solution_gauss_seidel(vector, MAX_ITERATIONS, HIGH_ACCURACY);
        }

        void print(uint32_t i) {
            printf("SOLUTIONS: %d\n", i);
            if (this->direct[i] != NULL) {
                printf("Direct\n");
                printf("Time spent on execution = %lf seconds.\n", this->direct[i]->execution_time);
            }
            if (this->jacobi[i] != NULL) {
                printf("Jacobi\n");
                printf("Time spent on execution = %lf seconds.\n", this->jacobi[i]->execution_time);
                printf("Number of iterations = %u.\n", this->jacobi[i]->iterations);
            }
            if (this->gauss_seidel[i] != NULL) {
                printf("Gauss Seidel\n");
                printf("Time spent on execution = %lf seconds.\n", this->gauss_seidel[i]->execution_time);
                printf("Number of iterations = %u.\n", this->gauss_seidel[i]->iterations);
            }
            printf("\n");
        }

        void print() {
            for (uint32_t i = 0; i < this->size; i += 1) {
                this->print(i);
            }
        }

        void save() {
            for (uint32_t i = 0; i < this->size; i += 1) {
                
                if (this->direct[i] != NULL) {
                    std::string file_name = "save/data_direct" + std::to_string(i) + std::to_string(this->size) + ".txt";
                    std::ofstream outFile(file_name);
                    outFile << this->direct[i]->execution_time << "\n";
                }
                if (this->jacobi[i] != NULL) {
                    std::string file_name = "save/data_jaco" + std::to_string(i) + std::to_string(this->size) + ".txt";
                    std::ofstream outFile(file_name);
                    outFile << this->jacobi[i]->execution_time << "\n";
                    outFile << this->jacobi[i]->iterations << "\n";
                    for (uint32_t j= 0; j < this->jacobi[i]->iterations; j += 1) {
                        outFile << this->jacobi[i]->norm_history[j] << "\n";
                    }
                    
                }
                if (this->gauss_seidel[i] != NULL) {
                    std::string file_name = "save/data_gaus" + std::to_string(i) + std::to_string(this->size) + ".txt";
                    std::ofstream outFile(file_name);
                    outFile << this->gauss_seidel[i]->execution_time << "\n";
                    outFile << this->gauss_seidel[i]->iterations << "\n";
                    for (uint32_t j= 0; j < this->gauss_seidel[i]->iterations; j += 1) {
                        outFile << this->gauss_seidel[i]->norm_history[j] << "\n";
                    }
                }
                
            }
        }
        
};

int main() {
    uint32_t index = 193285;

    uint32_t c = index % 100 / 10;
    uint32_t d = index % 10;
    uint32_t e = index % 1000 / 100;
    uint32_t f = index % 10000 / 1000;
    
    uint32_t N = 900 + index % 100;

    float64_t a1 = 5 + e, a2 = -1, a3 = -1;

    Matrix a_matrix(N, N);
    Matrix b_matrix(1, N);
    Matrix c_matrix(N, N);

    for (uint32_t y = 0; y < b_matrix.height; y += 1) {
        b_matrix.set(0, y, sin(y * (f + 1)));
    }

    make_three_arg_diagonal(&a_matrix, a1, a2, a3);
    make_three_arg_diagonal(&c_matrix, 3, a2, a3);

    DirectSolution *solv = c_matrix.solution_direct(&b_matrix);
    Matrix *norm = c_matrix.clone();
    c_matrix.mutiply_into(solv->matrix, norm);
    norm->print(1,3);
    norm->sub(&b_matrix);
    norm->print(1,3);
    float64_t norm_val = norm->norm();

    printf("%lf", norm_val);

    //b_matrix.print();
    return 0;
    
    

    uint32_t e_matrices_count= 8;

    Matrix **e_matrices = new Matrix*[e_matrices_count];
    Matrix **e_b_matrices = new Matrix*[e_matrices_count];

    uint32_t sizes[] = {100, 500, 1000, 2000, 3000, 4000, 5000, 6000};
    
    for (uint32_t i = 0; i < e_matrices_count; i += 1) {
        e_matrices[i] = new Matrix(sizes[i], sizes[i]);
        e_b_matrices[i] = new Matrix(1, sizes[i]);
        make_three_arg_diagonal(e_matrices[i], a1, a2, a3);
        for (uint32_t y = 0; y < sizes[i]; y += 1) {
            e_b_matrices[i] -> set(0, y, sin(y * (f + 1)));
        }
    }

    Solutions solutions_a_c(2);
    solutions_a_c.solve(0, &a_matrix, &b_matrix);
    solutions_a_c.solve(1, &c_matrix, &b_matrix);
    solutions_a_c.print();
    solutions_a_c.save();
    Solutions solutions_e(e_matrices_count);
    for (uint32_t i = 0; i < e_matrices_count; i += 1) {
        solutions_e.solve(i, e_matrices[i], e_b_matrices[i]);
        solutions_e.print(i);
    }
    solutions_e.print();
    solutions_e.save();
    return 0;
}