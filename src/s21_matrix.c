#include "s21_matrix.h"

double s21_determinant_calc(matrix_t *A);
void s21_minor_matrix(matrix_t *A, matrix_t *res, int str, int stb);
int s21_sign(int i, int j);

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error = OK;
  if (result == NULL || rows <= 0 || columns <= 0) {
    error = incorrect_matrix;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)malloc(rows * sizeof(double *));
    if (result->matrix != NULL) {
      for (int i = 0; i < rows; i++)
        result->matrix[i] = (double *)malloc(columns * sizeof(double));
    }
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = 0;
      }
    }
  }
  return error;
}

void s21_remove_matrix(matrix_t *A) {
  if (A != NULL) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
    A->matrix = NULL;
    A->columns = 0;
    A->rows = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int eq = SUCCESS;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0 ||
      B == NULL || B->matrix == NULL || B->rows <= 0 || B->columns <= 0) {
    eq = FAILURE;
  } else {
    if (A->rows == B->rows && A->columns == B->columns) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-07) {
            eq = FAILURE;
            i = A->columns;
            break;
          }
        }
      }
    } else
      eq = FAILURE;
  }
  return eq;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0 ||
      B == NULL || B->matrix == NULL || B->rows <= 0 || B->columns <= 0) {
    error = incorrect_matrix;
  } else {
    if (A->rows == B->rows && A->columns == B->columns) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else
      error = calculation_error;
  }
  return error;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0 ||
      B == NULL || B->matrix == NULL || B->rows <= 0 || B->columns <= 0) {
    error = incorrect_matrix;
  } else {
    if (A->rows == B->rows && A->columns == B->columns) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    } else
      error = calculation_error;
  }

  return error;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error = OK;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0) {
    error = incorrect_matrix;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return error;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0) {
    error = incorrect_matrix;
  } else {
    if (A->columns == B->rows) {
      s21_create_matrix(A->rows, B->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          for (int k = 0; k < A->columns; k++) {
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
          }
        }
      }
    } else
      error = calculation_error;
  }
  return error;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int error = OK;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0) {
    error = incorrect_matrix;
  } else {
    s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = A->matrix[j][i];
      }
    }
  }
  return error;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error = OK;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0) {
    error = incorrect_matrix;
  } else {
    if (A->columns == 1 && A->rows == 1) {
      s21_create_matrix(A->rows, A->columns, result);
      result->matrix[0][0] = A->matrix[0][0];
    } else {
      if (A->columns == A->rows) {
        matrix_t tmp_matrix = {0};
        int size = A->columns;
        int newsize = size - 1;
        s21_create_matrix(A->rows, A->columns, result);
        for (int x = 0; x < A->rows; x++) {
          for (int y = 0; y < A->columns; y++) {
            s21_create_matrix(newsize, newsize, &tmp_matrix);
            s21_minor_matrix(A, &tmp_matrix, x, y);
            result->matrix[x][y] =
                s21_sign(x, y) * s21_determinant_calc(&tmp_matrix);
            s21_remove_matrix(&tmp_matrix);
          }
        }

      } else {
        error = calculation_error;
      }
    }
  }
  return error;
}

void s21_minor_matrix(matrix_t *A, matrix_t *res, int str, int stb) {
  int tmp_j = 0;
  int tmp_i = 0;
  for (int x = 0; x < A->rows; x++) {
    for (int y = 0; y < A->columns; y++) {
      if (str != x && stb != y) {
        res->matrix[tmp_i][tmp_j] = A->matrix[x][y];
        tmp_j++;
      }
    }
    tmp_j = 0;
    if (str != x) tmp_i++;
  }
}

int s21_sign(int i, int j) {
  int sign = 1;
  if ((i + j) % 2) {
    sign = -1;
  }
  return sign;
}

int s21_determinant(matrix_t *A, double *result) {
  int error = OK;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0) {
    error = incorrect_matrix;
  } else if (A->columns != A->rows) {
    error = calculation_error;
  } else
    *result = s21_determinant_calc(A);
  return error;
}

double s21_determinant_calc(matrix_t *A) {
  double result = 0;
  matrix_t tmp_matrix = {0};
  if (A->columns == 2) {
    result = (A->matrix[0][0] * A->matrix[1][1]) -
             (A->matrix[0][1] * A->matrix[1][0]);
  } else if (A->columns == 1) {
    result = A->matrix[0][0];
  } else {
    int sign = 1;
    int size = A->columns;
    int newsize = size - 1;
    s21_create_matrix(newsize, newsize, &tmp_matrix);
    for (int x = 0; x < A->columns; x++) {
      int tmp_i = 0;
      for (int i = 1; i < A->columns; i++) {
        int tmp_j = 0;
        for (int j = 0; j < A->columns; j++) {
          if (j == x) continue;
          tmp_matrix.matrix[tmp_i][tmp_j] = A->matrix[i][j];
          tmp_j++;
        }
        tmp_i++;
      }
      result += sign * A->matrix[0][x] * s21_determinant_calc(&tmp_matrix);
      sign = -sign;
    }
    s21_remove_matrix(&tmp_matrix);
  }
  return result;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error = OK;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0) {
    error = incorrect_matrix;
  } else {
    if (A->rows == 1 && A->columns == 1) {
      s21_create_matrix(A->rows, A->columns, result);
      result->matrix[0][0] = 1 / A->matrix[0][0];
    } else {
      if (A->columns == A->rows) {
        double det = 0;
        det = 1 / s21_determinant_calc(A);
        matrix_t tmp_matrix_minor = {0};
        matrix_t tmp_matrix_transport = {0};
        s21_calc_complements(A, &tmp_matrix_minor);
        s21_transpose(&tmp_matrix_minor, &tmp_matrix_transport);
        s21_mult_number(&tmp_matrix_transport, det, result);
        s21_remove_matrix(&tmp_matrix_transport);
        s21_remove_matrix(&tmp_matrix_minor);
      } else
        error = calculation_error;
    }
  }
  return error;
}
