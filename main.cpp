
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <cmath>
#include <chrono>

#include "types.hpp"
#include "matrix.hpp"
#include "solutions.hpp"


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

    make_three_arg_diagonal(&a_matrix, a1, a2, a3);
    for (uint32_t y = 0; y < b_matrix.height; y += 1) {
        b_matrix.set(0, y, sin(y * (f + 1)));
    }


    Matrix c_matrix(N, N);
    make_three_arg_diagonal(&c_matrix, 3, a2, a3);
    // float64_t test2[] {
    //     2, 1 
    //     -1, 4
    //      -2, 1
    //     -4, 1
    // };

    // float64_t test1[] {
    //     2, -1, -2, 
    //     -4, 2, 7, 2,
    //     8, 3 ,2, 0,
    //     1 , 2, 3 ,5
    // };

    // Matrix gauss_matrix_1(4, 4, test1);
    // Matrix gauss_matrix_2(4, 4, test1);
    // // make_three_arg_diagonal(&gauss_matrix_1, 3, a2, a3);
    // // make_three_arg_diagonal(&gauss_matrix_2, 6, 1, a3);
    // Matrix *gauss1 = gauss_matrix_1.solution_lu(&gauss_matrix_2);
    // Matrix *gauss2 = gauss_matrix_1.solution(&gauss_matrix_2);


    // Matrix *check1 = gauss_matrix_1.clone();
    // check1->mutiply(gauss1);
    // Matrix *check2 = gauss_matrix_1.clone();
    // check2->mutiply(gauss2);

    // gauss1->print(); printf("\n");
    // gauss2->print(); printf("\n");

    // check1->print(); printf("\n");
    // check2->print(); printf("\n");
    // return(0);


    // printf("\n");
    // DirectSolution *solution1 = a_matrix.solution_direct(&b_matrix);
    // solution1->matrix->print(1,5);
    // printf("Time spent on execution = %u seconds.\n", solution1->execution_time);
    // delete solution1;

    // printf("\n");
    // IterativeSolution *solution2 = a_matrix.solution_jacobi(&b_matrix, 100, 1.0e-10);
    // solution2->matrix->print(1,5);
    // printf("Time spent on prep = %u seconds.\n", solution2->prep_time);
    // printf("Time spent on execution = %u seconds.\n", solution2->execution_time);
    // printf("Number of iterations = %u.\n", solution2->iterations);
    // delete solution2;

    // printf("\n");
    // IterativeSolution *solution3 = a_matrix.solution_gauss_seidel(&b_matrix, 100, 1.0e-10);
    // solution3->matrix->print(1,5);
    // printf("Time spent on prep = %u seconds.\n", solution3->prep_time);
    // printf("Time spent on execution = %u seconds.\n", solution3->execution_time);
    // printf("Number of iterations = %u.\n", solution3->iterations);
    // delete solution3;


    printf("\n");
    DirectSolution *solution4 = c_matrix.solution_direct(&b_matrix);
    solution4->matrix->print(1,5);
    printf("Time spent on execution = %u seconds.\n", solution4->execution_time);
    delete solution4;

    printf("\n");
    IterativeSolution *solution5 = c_matrix.solution_jacobi(&b_matrix, 100, 1.0e-10);
    solution5->matrix->print(1,5);
    printf("Time spent on prep = %u seconds.\n", solution5->prep_time);
    printf("Time spent on execution = %u seconds.\n", solution5->execution_time);
    printf("Number of iterations = %u.\n", solution5->iterations);
    delete solution5;

    printf("\n");
    IterativeSolution *solution6 = c_matrix.solution_gauss_seidel(&b_matrix, 100, 1.0e-10);
    solution6->matrix->print(1,5);
    printf("Time spent on prep = %u seconds.\n", solution6->prep_time);
    printf("Time spent on execution = %u seconds.\n", solution6->execution_time);
    printf("Number of iterations = %u.\n", solution6->iterations);
    delete solution6;

    return 0;
}