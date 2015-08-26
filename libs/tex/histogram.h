/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_HISTOGRAM_HEADER
#define TEX_HISTOGRAM_HEADER

#include <vector>
#include <string>

/**
  * Class representing a histogram with a fixed number of bins
  * optimized to calculate approximate permilles.
  */
class Histogram {
    private:
        std::vector<unsigned int> bins;
        float min;
        float max;
        int num_values;

    public:
        /** Constructs a histogram with num_bins bins which clamps values to [_min, _max]. */
        Histogram(float _min, float _max, std::size_t num_bins);

        /** Adds a value to the histogram. The value is clamped to [min, max]. */
        void add_value(float value);

        /**
          * Saves the histogram to a .csv file.
          * @throws util::FileException
          */
        void save_to_file(std::string const & filename) const;

        /**
          * Returns the approximate permille.
          */
        float get_approx_percentile(float percentile) const;
};

#endif /* TEX_HISTOGRAM_HEADER */
