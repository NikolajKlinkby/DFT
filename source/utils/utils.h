#ifndef DFT_UTILS_H
#define DFT_UTILS_H

template <typename scalar_t>
void vector_set(slate::Matrix<scalar_t> vector, int index, scalar_t value);

template <typename scalar_t>
scalar_t vector_get(slate::Matrix<scalar_t> vector, int index);

template <typename scalar_t>
scalar_t vector_get_max(slate::Matrix<scalar_t> vector);

template <typename scalar_t>
scalar_t vector_get_min(slate::Matrix<scalar_t> vector);

template <typename scalar_t>
void interp_deriv(slate::Matrix<scalar_t> x,
                  slate::Matrix<scalar_t> y,
                  slate::Matrix<scalar_t> y_deriv,
                  slate::Matrix<scalar_t> x_mat,
                  slate::Matrix<scalar_t> x_mat_deriv);

#endif //DFT_UTILS_H
