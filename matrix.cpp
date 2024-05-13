#include "matrix.hpp"
#include <stdint.h>
#include "types.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <chrono>

typedef std::chrono::high_resolution_clock fast_clock_t;
typedef std::chrono::time_point<fast_clock_t> time_point_t;

time_point_t get_time() {
    return fast_clock_t::now();
}

uint32_t calc_time_diff_in_second(time_point_t start, time_point_t end) {
    return std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
}



Matrix::Matrix(uint32_t width, uint32_t height) {
    this->height = height;
    this->width = width;
    this->content = new float64_t[width * height]; 
    for (uint32_t i = 0; i < this->width * this->height; i += 1) {
        this->content[i] = 0;
    }
}

Matrix::Matrix(uint32_t width, uint32_t height, float64_t *content) {
    this->height = height;
    this->width = width;
    this->content = new float64_t[width * height]; 
    for (uint32_t x = 0; x < this->width; x += 1) {
        for (uint32_t y = 0; y < this->height; y += 1) {
            this->set(x, y, content[y * this->width + x]);
        }
    }
}

Matrix::Matrix(uint32_t width, uint32_t height, float64_t val) {
    this->height = height;
    this->width = width;
    this->content = new float64_t[width * height]; 
    for (uint32_t x = 0; x < this->width; x += 1) {
        for (uint32_t y = 0; y < this->height; y += 1) {
            this->set(x, y, val);
        }
    }
}

Matrix::Matrix(MatrixReadOnly *matrix) {
    this->height = matrix->height;
    this->width = matrix->width;
    this->content = new float64_t[width * height]; 
    for (uint32_t x = 0; x < this->width; x += 1) {
        for (uint32_t y = 0; y < this->height; y += 1) {
            this->set(x, y, matrix->get(x, y));
        }
    }
}


Matrix::~Matrix() {
    delete[] this->content;
}

void MatrixReadOnly::print() {
    for (uint32_t y = 0; y < this->height; y += 1) {
        for (uint32_t x = 0; x < this->width; x += 1) {
            printf("%lf,", this->get(x, y));
        }
        printf("\n");
    }
}

void MatrixReadOnly::print(uint32_t x_limit, uint32_t y_limit) {
    for (uint32_t y = 0; y < y_limit; y += 1) {
        for (uint32_t x = 0; x < x_limit; x += 1) {
            printf("%lf,", this->get(x, y));
        }
        printf("\n");
    }
}

float64_t MatrixReadOnly::get(uint32_t x, uint32_t y) {
    return this->content[y * this->width + x];
}

void Matrix::set(uint32_t x, uint32_t y, float64_t val) {
    this->content[y * this->width + x] = val;
}

void Matrix::add(MatrixReadOnly *matrix) {
    for (uint32_t x = 0; x < this->width; x += 1) {
        for (uint32_t y = 0; y < this->height; y += 1) {
            this->set(x, y, this->get(x, y) + matrix->get(x, y));
        }
    }
}

void Matrix::sub(MatrixReadOnly *matrix) {
    for (uint32_t x = 0; x < this->width; x += 1) {
        for (uint32_t y = 0; y < this->height; y += 1) {
            this->set(x, y, this->get(x,y) - matrix->get(x, y));
        }
    }
}

Matrix *MatrixReadOnly::clone() {
    return new Matrix(this->width, this->height, this->content);
}


void Matrix::copy(MatrixReadOnly *matrix) {
    this->height = matrix->height;
    this->width = matrix->width;
    this->content = (float64_t *)realloc(
        this->content,
        matrix->height * matrix->width * sizeof(float64_t)
    );
    for (uint32_t x = 0; x < matrix->width; x += 1) {
        for (uint32_t y = 0; y < matrix->height; y += 1) {
            this->set(x, y, matrix->get(x, y));
        }
    }
}

void Matrix::transform(MatrixType type) {
    switch(type) {
        case LOWER_TRIANGLE:
            for (uint32_t x = 0, y = 0; x < this->width; x += 1, y += 1) {
                for (uint32_t next = x; next < this->height; next += 1) {
                    this->set(next, y, 0);
                }
            }
            break;
        case UPPER_TRIANGLE:
            for (uint32_t x = 0, y = 0; x < this->width; x += 1, y += 1) {
                for (int32_t next = x; next >= 0; next -= 1) {
                    this->set(next, y, 0);
                }
            }
            break;
        case IDENTITY:
            for (uint32_t x = 0; x < this->width; x += 1) {
                for (uint32_t y = 0; y < this->height; y += 1) {
                    if (x == y) {
                        this->set(x, y, 1);
                    } else {
                        this->set(x, y, 0);
                    }
                }
            }
            break;
        case DIAGONAL:
            for (uint32_t x = 0; x < this->width; x += 1) {
                for (uint32_t y = 0; y < this->height; y += 1) {
                    if (x != y) {
                        this->set(x, y, 0);
                    }
                }
            }
            break;
    }
}


void Matrix::mutiply(MatrixReadOnly *matrix) {
    uint32_t m = this->height;
    uint32_t n = matrix->width;
    uint32_t p = matrix->height;

    float64_t buffer;

    Matrix *copy = this->clone();

    if (this->height != m || this->width != n) {   
        this->height = m;
        this->width = n;
        this->content = (float64_t *)realloc(this->content, m * n * sizeof(float64_t));
    }

    for (uint32_t x = 0; x < m; x += 1) {
        for (uint32_t y = 0; y < n; y += 1) {
            buffer = 0;
            for (uint32_t k = 0; k < p; k += 1) {
                buffer += copy->get(x, k) * matrix->get(k, y);
            }
            this->set(x, y, buffer);
        }
    }
    delete copy;
}

void MatrixReadOnly::mutiply_into(MatrixReadOnly *matrix, Matrix *dest) {
    uint32_t m = this->height;
    uint32_t n = matrix->width;
    uint32_t p = matrix->height;

    float64_t buffer;

    if (dest->height != m || dest->width != n) {    
        dest->height = m;
        dest->width = n;
        dest->content = (float64_t *)realloc(dest->content, m * n * sizeof(float64_t));
    }

    for (uint32_t x = 0; x < m; x += 1) {
        for (uint32_t y = 0; y < n; y += 1) {
            buffer = 0;
            for (uint32_t k = 0; k < p; k += 1) {
                buffer += this->get(x, k) * matrix->get(k, y);
            }
            dest->set(x, y, buffer);
        }
    }
}

void Matrix::dot(float64_t val) {
    for (uint32_t x = 0; x < this->width; x += 1) {
        for (uint32_t y = 0; y < this->height; y += 1) {
            this->set(x, y, this->get(x, y) * val);
        }
    }
}

void Matrix::swap_rows(uint32_t y_row1, uint32_t y_row2) {
    float64_t *buffer = new float64_t[this->width];

    for (uint32_t x = 0; x < this->width; x += 1) {
        buffer[x] = this->get(x, y_row1);
    }

    for (uint32_t x = 0; x < this->width; x += 1) {
        this->set(x, y_row1, this->get(x, y_row2));
    }

    for (uint32_t x = 0; x < this->width; x += 1) {
        this->set(x, y_row2, buffer[x]);
    }

    delete[] buffer;
}

Matrix *MatrixReadOnly::solution(MatrixReadOnly *matrix) {
    Matrix *gauss = this->gauss_matrix(matrix);

    Matrix *solution = new Matrix(matrix->width, this->height);

    for (uint32_t x = 0; x < matrix->width; x += 1) {
        for (uint32_t y = 0; y < this->height; y += 1) {
            solution->set(x, y, gauss->get(x + this->width, y));
        }
    }

    delete gauss;
    return solution;
}

Matrix *MatrixReadOnly::solution_lu(MatrixReadOnly *matrix) {


    Matrix *L = new Matrix(this->width, this->height);
    Matrix *U = new Matrix(this->width, this->height);

    float64_t sum;

    for (uint32_t p = 0; p < this->width; p += 1){
        for (uint32_t k = p; k < this->width; k += 1) {
            sum = 0.0;
            for (uint32_t l = 0; l < p; l += 1) {
                sum += (L->get(l, p) * U->get(k, l));
            }
            U->set(k, p, this->get(k, p) - sum);
        }

        for (uint32_t k = p; k < this->width; k += 1) {
            if (p == k) {
                L->set(p, p, 1); 
            } else {
                sum = 0.0;
                for (uint32_t l = 0; l < p; l += 1) {
                    sum += (L->get(l, k) * U->get(p, l));
                }
                L->set(p, k, (this->get(p, k) - sum) / U->get(p, p));
            }
        }
    }
    
    
    Matrix *solution_l = new Matrix(this->width, this->height);
    float64_t divider;
    for (uint32_t y = 0; y < this->width; y += 1){
        for (int32_t x = y - 1; x >= 0 ; x -= 1) {
            for (int32_t k = 0; k < matrix->width ; k += 1) {
                solution_l->set(k, y, solution_l->get(k, y) + solution_l->get(k, x) * L->get(x, y));
            }
        }
        
        divider = L->get(y, y);
        for (int32_t k = 0; k < matrix->width ; k += 1) {
            solution_l->set(k, y, (matrix->get(k, y) - solution_l->get(k, y)) / divider);
        }
    }

    Matrix *solution_u = new Matrix(this->width, this->height);
    for (int32_t y = this->width - 1; y >= 0; y -= 1){
        for (uint32_t x = y + 1; x < this->width ; x += 1) {
            for (int32_t k = 0; k < matrix->width ; k += 1) {
                solution_u->set(k, y, solution_u->get(k, y) + solution_u->get(k, x) * U->get(x, y));
            }
            
        }
        divider = U->get(y, y);
        for (int32_t k = 0; k < matrix->width ; k += 1) {
            solution_u->set(k, y, (solution_l->get(k, y) - solution_u->get(k, y)) / divider);
        }
    }


    delete L;
    delete U;
    return solution_u;
}

DirectSolution *MatrixReadOnly::solution_direct(MatrixReadOnly *matrix) {

    time_point_t start = get_time();
    
    Matrix *solution = this->solution_lu(matrix);

    time_point_t end = get_time();
    
    DirectSolution *output = new DirectSolution(solution);
    output->execution_time = calc_time_diff_in_second(start, end); 

    return output;
}

float64_t MatrixReadOnly::norm() {
    float64_t sum = 0.0;
    for (uint32_t y = 0; y < this->height; y += 1) {
        sum += this->get(0, y) * this->get(0, y);
    }
    return sqrt(sum);
}

IterativeSolution *MatrixReadOnly::solution_jacobi(MatrixReadOnly *input, uint32_t iter, float64_t err_norm) {
    time_point_t start, end;

    start = get_time();

    // D, U, L
    MatrixShapeLookup diag(this, DIAGONAL);
    MatrixShapeLookup tri_upper(this, UPPER_TRIANGLE);
    MatrixShapeLookup tri_lower(this, LOWER_TRIANGLE);

    // L + U
    Matrix *tri_sum = tri_upper.clone(); 
    tri_sum->add(&tri_lower);
    // M = -(D \ (L + U))
    Matrix *m_matrix = diag.solution_lu(tri_sum);
    m_matrix->dot(-1);
    // bm = D \ b
    Matrix *bm_vector = diag.solution_lu(input);

    Matrix *solution = new Matrix(1, this->height);
    Matrix *next_solution = solution->clone();
    Matrix *norm_matrix = this->clone();
    float64_t norm_val;

    end = get_time();
    IterativeSolution *output = new IterativeSolution(solution, iter);
    output->prep_time = calc_time_diff_in_second(start, end);

    start = end;
    uint32_t i = 0;
    for (; i < iter; i += 1) {
        // x = M * x + bm
        m_matrix->mutiply_into(solution, next_solution);
        next_solution->add(bm_vector);
        
        std::swap(solution, next_solution);

        // A * x - b
        this->mutiply_into(solution, norm_matrix);
        norm_matrix->sub(input);

        norm_val = norm_matrix->norm();
        output->norm_history[i] = norm_val;

        if (norm_val <= err_norm) {
            break;
        }
    }
    end = get_time();
    output->execution_time = calc_time_diff_in_second(start, end);
    output->matrix = solution;
    output->iterations = i;

    delete tri_sum;
    delete m_matrix;
    delete bm_vector;
    delete next_solution;
    delete norm_matrix;
    return output;
}


IterativeSolution *MatrixReadOnly::solution_gauss_seidel(MatrixReadOnly *input, uint32_t iter, float64_t err_norm) {
    time_point_t start, end;

    start = get_time();
    
    // D, U, L
    MatrixShapeLookup diag(this, DIAGONAL);
    MatrixShapeLookup tri_upper(this, UPPER_TRIANGLE);
    MatrixShapeLookup tri_lower(this, LOWER_TRIANGLE);
    // D + L
    Matrix *diag_tri_lower_sum = diag.clone();
    diag_tri_lower_sum->add(&tri_lower);
    // M = -((D + L) \ U)
    Matrix *m_matrix = diag_tri_lower_sum->solution_lu(&tri_upper);
    m_matrix->dot(-1);
    // bm = (D + L) \ b
    Matrix *bm_vector = diag_tri_lower_sum->solution_lu(input);

    Matrix *solution = new Matrix(1, this->height);
    Matrix *next_solution = solution->clone();
    Matrix *norm_matrix = this->clone();
    float64_t norm_val;

    end = get_time();
    IterativeSolution *output = new IterativeSolution(solution, iter);
    output->prep_time = calc_time_diff_in_second(start, end);

    start = end;
    uint32_t i = 0;
    for (; i < iter; i += 1) {
        m_matrix->mutiply_into(solution, next_solution);
        next_solution->add(bm_vector);
        
        std::swap(solution, next_solution);
        
        this->mutiply_into(solution, norm_matrix);
        norm_matrix->sub(input);

        norm_val = norm_matrix->norm();
        output->norm_history[i] = norm_val;

        if (norm_val <= err_norm) {
            break;
        }
    }

    end = get_time();
    output->execution_time = calc_time_diff_in_second(start, end);
    output->matrix = solution;
    output->iterations = i;

    delete diag_tri_lower_sum;
    delete m_matrix;
    delete bm_vector;
    delete next_solution;
    delete norm_matrix;
    return output;
}

void Matrix::invert() {
    Matrix identity(this->width, this->height);

    identity.transform(IDENTITY);

    Matrix *gauss = gauss_matrix(&identity);
    
    for (uint32_t x = 0; x < this->width; x += 1) {
        for (uint32_t y = 0; y < this->height; y += 1) {
            this->set(x, y, gauss->get(x + this->width, y));
        }
    }
    delete gauss;
}

void Matrix::transpose() {
    Matrix transposition(this->height, this->width);

    for (uint32_t x = 0; x < this->width; x += 1) {
        for (uint32_t y = 0; y < this->height; y += 1) {
            transposition.set(y, x, this->get(x, y));
        }
    }

    this->copy(&transposition);
}

bool Matrix::is_invertible() {
    Matrix identity(this->width, this->height);

    identity.transform(IDENTITY);

    Matrix *gauss = gauss_matrix(&identity);

    bool is_invertiable = true;
    for (uint32_t diag = 0; diag < gauss->height; diag +=1) {
        if (gauss->get(diag, diag) != 1.0) {
            is_invertiable = false;
            break;
        }
    }

    return is_invertiable;
}

void Matrix::divide(MatrixReadOnly *matrix) {
    Matrix *copy = new Matrix(matrix);
    copy->invert();
    this->mutiply(copy);
    delete copy;
}

Matrix *MatrixReadOnly::gauss_matrix(MatrixReadOnly *matrix) {

    Matrix *gauss = new Matrix(this->width + matrix->width, this->height);

    for (uint32_t x = 0; x < this->width; x += 1) {
        for (uint32_t y = 0; y < this->height; y += 1) {
            gauss->set(x, y, this->get(x, y));
        }
    }
    for (uint32_t x = 0; x < matrix->width; x += 1) {
        for (uint32_t y = 0; y < matrix->height; y += 1) {
            gauss->set(this->width + x, y, matrix->get(x, y));
        }
    }

    uint32_t selected_row;
    for (uint32_t x = 0; x < gauss->width - 1 && x < gauss->height; x += 1) {
        selected_row = UINT32_MAX;
        for (uint32_t y = x; y < gauss->height; y += 1) {
            if (gauss->get(x,y) != 0) {
                selected_row = y;
                break;
            }
        }

        if(selected_row == UINT32_MAX) {
            continue;
        }

        float64_t multiplier;
        for (uint32_t y = 0; y < gauss->height; y += 1) {
            if (y == selected_row) {
                continue;
            }

            multiplier = gauss->get(x, y) / gauss->get(x, selected_row);

            for (uint32_t p = 0; p < gauss->width; p += 1) {
                gauss->set(p, y, gauss->get(p, y) - gauss->get(p, selected_row) * multiplier);
            }
        }
        gauss->swap_rows(selected_row, x);
    }

    for(uint32_t diag = 0; diag < gauss->height; diag += 1) {
        for(uint32_t x = this->width; x <= gauss->width; x += 1) {
            gauss->set(x, diag, gauss->get(x, diag) / gauss->get(diag, diag));
        }
        gauss->set(diag, diag, 1.0);
    }
    return gauss;
}




MatrixShapeLookup::MatrixShapeLookup(MatrixReadOnly *matrix, MatrixType type) {
    this->height = matrix->height;
    this->width = matrix->width;
    this->content = matrix->content;
    this->type = type;
}

float64_t MatrixShapeLookup::get(uint32_t x, uint32_t y) {
   switch (this->type) {
        case LOWER_TRIANGLE:
            if (x > y) {
                return 0.0;
            } else {
                return MatrixReadOnly::get(x, y);
            }
            break;
        case UPPER_TRIANGLE:
            if (x < y) {
                return 0.0;
            } else {
                return MatrixReadOnly::get(x, y);
            }
            break;
        case DIAGONAL:
            if (x != y) {
                return 0.0;
            } else {
                return MatrixReadOnly::get(x, y);
            }
            break;
   }
   return 0.0;
}

void MatrixShapeLookup::print() {
    for (uint32_t y = 0; y < this->height; y += 1) {
        for (uint32_t x = 0; x < this->width; x += 1) {
            printf("%lf,", this->get(x, y));
        }
        printf("\n");
    }
}

void MatrixShapeLookup::print(uint32_t x_limit, uint32_t y_limit) {
    for (uint32_t y = 0; y < y_limit; y += 1) {
        for (uint32_t x = 0; x < x_limit; x += 1) {
            printf("%lf,", this->get(x, y));
        }
        printf("\n");
    }
}