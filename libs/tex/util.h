/*
 * Copyright (C) 2015, Nils Moehrle, Michael Waechter
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_UTIL_HEADER
#define TEX_UTIL_HEADER

#include <string>
#include <vector>
#include <fstream>
#include <cstring>
#include <cerrno>
#include <cassert>

#include <Eigen/Core>

#include "util/exception.h"
#include "util/file_system.h"

#include "math/vector.h"
#include "math/matrix.h"
#include "math/functions.h"

/**
 * Converts an MVE matrix into an Eigen matrix.
 */
template <typename T, int M, int N> Eigen::Matrix<T,M,N>
mve_to_eigen(math::Matrix<T, M, N> const & mat) {
    Eigen::Matrix<T, M, N> ret;
    for (int m = 0; m < M; ++m)
        for (int n = 0; n < N; ++n)
            ret(m, n) = mat(m, n);
    return ret;
}

/**
 * Converts an MVE vector into an Eigen row vector.
 */
template <typename T, int N> Eigen::Matrix<T, 1, N>
mve_to_eigen(math::Vector<T, N> const & vec) {
    Eigen::Matrix<T, 1, N> ret;
    for (int n = 0; n < N; ++n)
        ret(0, n) = vec(n);
    return ret;
}

/**
 * A multi-variate normal distribution which is NOT normalized such that the integral is 1.
 *
 * @param X is the vector for which the function is to be evaluated. Must have size 1xN.
 * @param mu is the mean around which the distribution is centered. Must have size 1xN.
 * @param covariance_inv is the INVERSE of the covariance matrix. Must have size NxN.
 * @return \f$\exp(-\frac{1}{2} (X-\mbox{mu})^T \cdot \mbox{covariance\_inv} \cdot (X-\mbox{mu}))\f$
 */
template <typename T, int N> T const
multi_gauss_unnormalized(Eigen::Matrix<T, 1, N> const & X, Eigen::Matrix<T, 1, N> const & mu,
    Eigen::Matrix<T, N, N> const & covariance_inv) {

    Eigen::Matrix<T, 1, N> mean_removed = X - mu;
    return std::exp(T(-0.5) * mean_removed * covariance_inv * mean_removed.adjoint());
}

/** Return the number suffix for n, e.g. 3 -> "rd". */
inline
std::string number_suffix(int n) {
    if (n % 100 == 11 || n % 100 == 12 || n % 100 == 13)
        return "th";
    if (n % 10 == 1)
        return "st";
    if (n % 10 == 2)
        return "nd";
    if (n % 10 == 3)
        return "rd";
    return "th";
}

/**
  * Write vector to csv file with "Index, $header" as header.
  * @throws util::FileException
  */
template <typename T> void
write_vector_to_csv(std::string const & filename, std::vector<T> const & vector, std::string const & header) {
    std::ofstream out(filename.c_str());
    if (!out.good())
        throw util::FileException(filename, std::strerror(errno));

    out << "Index, "<< header << std::endl;
    for (std::size_t i = 0; i < vector.size(); i++){
        out << i << ", " << vector[i] << std::endl;
    }
    out.close();
}


/**
  * Write vector to binary file.
  * @throws util::FileException
  */
template <typename T> void
vector_to_file(std::string const & filename, std::vector<T> const & vector) {
    std::ofstream out(filename.c_str());
    if (!out.good())
        throw util::FileException(filename, std::strerror(errno));

    out.write(reinterpret_cast<const char*>(&vector[0]), vector.size()*sizeof(T));
    out.close();
}

/**
  * Loads vector from binary file.
  * @throws util::FileException
  */
template <typename T> std::vector<T>
vector_from_file(std::string const & filename) {
    std::ifstream in(filename.c_str());
    if (!in.good())
        throw util::FileException(filename, std::strerror(errno));
    in.seekg (0, in.end);
    const size_t filesize = in.tellg();
    in.seekg (0, in.beg);
    const size_t num_elements = filesize / sizeof(T);
    std::vector<T> vector(num_elements);
    in.read(reinterpret_cast<char*>(&vector[0]), num_elements*sizeof(T));
    in.close();
    return vector;
}

/**
  * Writes the given string to the file given by filename.
  * @throws util::FileException
  */
inline void
write_string_to_file(std::string const & filename, std::string const & string) {
    std::ofstream out(filename.c_str());
    if (!out.good())
        throw util::FileException(filename, std::strerror(errno));

    out << string;
    out.close();
}

/**
  * returns the corresponding jet color encoding for the given value.
  * @warning asserts values within the interval [0, 1]
  */
inline math::Vec4f
get_jet_color(float value) {
    assert(0.0f <= value && value <= 1.0f);
    float mvalue = 4 * value;
    float red   = math::clamp(std::min(mvalue - 1.5f, -mvalue + 4.5f));
    float green = math::clamp(std::min(mvalue - 0.5f, -mvalue + 3.5f));
    float blue  = math::clamp(std::min(mvalue + 0.5f, -mvalue + 2.5f));
    return math::Vec4f(red, green, blue, 1.0f);
}

#endif /* TEX_UTIL_HEADER */
