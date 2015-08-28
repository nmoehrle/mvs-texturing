/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <cmath>

#include "rectangular_bin.h"

RectangularBin::RectangularBin(unsigned int width, unsigned int height)
    : width(width), height(height) {
    rects.push_back(Rect<int>(0, 0, width, height));
}

bool RectangularBin::insert(Rect<int> * rect) {
    /* The best score is 0 so we initialize with the worst. */
    unsigned int best_score = width * height;
    std::list<Rect<int> >::iterator best_rect_it = rects.end();
    std::list<Rect<int> >::iterator it = rects.begin();
    for (; it != rects.end(); ++it) {
        Rect<int> free_rect = *it;
        if (rect->width() <= free_rect.width()
            && rect->height() <= free_rect.height() ) {
            unsigned int score = free_rect.size() - rect->size();
            if (score < best_score){
                best_score = score;
                best_rect_it = it;
            }
        }
    }

    /* Fits? */
    if (best_rect_it != rects.end()) {
        Rect<int> best_rect(&(*best_rect_it));
        rects.erase(best_rect_it);

        /* Update the rect. */
        rect->move(best_rect.min_x, best_rect.min_y);

        /* Decide split axis. */
        Rect<int> hsplit_top(best_rect.min_x, rect->max_y, best_rect.max_x, best_rect.max_y);
        Rect<int> hsplit_bottom(rect->max_x, best_rect.min_y, best_rect.max_x, rect->max_y);
        Rect<int> vsplit_left(best_rect.min_x, rect->max_y, rect->max_x, best_rect.max_y);
        Rect<int> vsplit_right(rect->max_x, best_rect.min_y, best_rect.max_x, best_rect.max_y);

        float hsplit_ratio = 1.0f;
        float vsplit_ratio = 1.0f;

        if (hsplit_top.size() != 0 && hsplit_bottom.size() != 0)
            hsplit_ratio = static_cast<float>(hsplit_top.size()) / hsplit_bottom.size();
        if (vsplit_left.size() != 0 && vsplit_right.size() != 0)
            vsplit_ratio = static_cast<float>(vsplit_left.size()) / vsplit_right.size();

        if (std::abs(1.0f - hsplit_ratio) < std::abs(1.0f - vsplit_ratio)){
            if (vsplit_left.size() != 0) rects.push_back(vsplit_left);
            if (vsplit_right.size() != 0) rects.push_back(vsplit_right);
        } else {
            if (hsplit_top.size() != 0) rects.push_back(hsplit_top);
            if (hsplit_bottom.size() != 0) rects.push_back(hsplit_bottom);
        }

        return true;
    } else {
        return false;
    }
}
