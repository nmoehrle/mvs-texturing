/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <cstdint>
#include <iostream>

#include <math/vector.h>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include "poisson_blending.h"

typedef Eigen::SparseMatrix<float> SpMat;

math::Vec3f simple_laplacian(int i, mve::FloatImage::ConstPtr img){
    const int width = img->width();
    assert(i > width + 1 && i < img->get_pixel_amount() - width -1);

    return -4.0f * math::Vec3f(&img->at(i, 0))
        + math::Vec3f(&img->at(i - width, 0))
        + math::Vec3f(&img->at(i - 1, 0))
        + math::Vec3f(&img->at(i + 1, 0))
        + math::Vec3f(&img->at(i + width, 0));
}

bool valid_mask(mve::ByteImage::ConstPtr mask){
    const int width = mask->width();
    const int height = mask->height();

    for (int x = 0; x < width; ++x)
        if (mask->at(x, 0, 0) == 255 || mask->at(x, height - 1, 0) == 255)
            return false;

    for (int y = 0; y < height; ++y)
        if (mask->at(0, y, 0) == 255 || mask->at(width - 1, y, 0) == 255)
            return false;

    //TODO check for sane boundary conditions...

    return true;
}

void
poisson_blend(mve::FloatImage::ConstPtr src, mve::ByteImage::ConstPtr mask,
    mve::FloatImage::Ptr dest, float alpha) {

    assert(src->width() == mask->width() && mask->width() == dest->width());
    assert(src->height() == mask->height() && mask->height() == dest->height());
    assert(src->channels() == 3 && dest->channels() == 3);
    assert(mask->channels() == 1);
    assert(valid_mask(mask));

    const int n = dest->get_pixel_amount();
    const int width = dest->width();
    const int height = dest->height();
    const int channels = dest->channels();

    mve::Image<int>::Ptr indices = mve::Image<int>::create(width, height, 1);
    indices->fill(-1);
    int index = 0;
    for (int i = 0; i < n; ++i) {
        if (mask->at(i) != 0) {
            indices->at(i) = index;
            index += 1;
        }
    }
    const int nnz = index;

    std::vector<math::Vec3f> coefficients_b;
    coefficients_b.resize(nnz);

    std::vector<Eigen::Triplet<float, int> > coefficients_A;
    coefficients_A.reserve(nnz); //TODO better estimate...

    for (int i = 0; i < n; ++i) {
        const int row = indices->at(i);
        if (mask->at(i) == 128 || mask->at(i) == 64) {
            Eigen::Triplet<float, int> t(row, row, 1.0f);
            coefficients_A.push_back(t);

            coefficients_b[row] = math::Vec3f(&dest->at(i, 0));
        }

        if (mask->at(i) == 255) {
            const int i01 = indices->at(i - width);
            const int i10 = indices->at(i - 1);
            const int i11 = indices->at(i);
            const int i12 = indices->at(i + 1);
            const int i21 = indices->at(i + width);

            /* All neighbours should be eighter border conditions or part of the optimization. */
            assert(i01 != -1 && i10 != -1 && i11 != -1 && i12 != -1 && i21 != -1);

            Eigen::Triplet<float, int> t01(row, i01, 1.0f);

            Eigen::Triplet<float, int> t10(row, i10, 1.0f);
            Eigen::Triplet<float, int> t11(row, i11, -4.0f);
            Eigen::Triplet<float, int> t12(row, i12, 1.0f);

            Eigen::Triplet<float, int> t21(row, i21, 1.0f);

            Eigen::Triplet<float, int> triplets[] = {t01, t10, t11, t12, t21};

            coefficients_A.insert(coefficients_A.end(), triplets, triplets + 5);

            math::Vec3f l_d = simple_laplacian(i, dest);
            math::Vec3f l_s = simple_laplacian(i, src);

            coefficients_b[row] = (alpha * l_s + (1.0f - alpha) * l_d);
        }
    }

    SpMat A(nnz, nnz);
    A.setFromTriplets(coefficients_A.begin(), coefficients_A.end());

    Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int> > solver;
    solver.compute(A);

    for (int channel = 0; channel < channels; ++channel) {
        Eigen::VectorXf b(nnz);
        for (std::size_t i = 0; i < coefficients_b.size(); ++i)
            b[i] = coefficients_b[i][channel];

        Eigen::VectorXf x(n);
        x = solver.solve(b);

        for (int i = 0; i < n; ++i) {
            int index = indices->at(i);
            if (index != -1) dest->at(i, channel) = x[index];
        }
    }
}
