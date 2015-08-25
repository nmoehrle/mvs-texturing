/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_RECT_HEADER
#define TEX_RECT_HEADER

#include <cassert>

/**
  * Simple class representing a rectangle.
  */
template <typename T>
class Rect {
public:
    T min_x;
    T min_y;
    T max_x;
    T max_y;

    Rect<T> (void);
    Rect<T> (Rect<T> * rect);
    Rect<T> (T min_x, T min_y, T max_x, T max_y);

    T width(void) const;
    T height(void) const;
    T size(void) const;

    /** Returns true if the rectangle is within or on the border of the given rectangle. */
    bool is_inside(Rect const * orect) const;
    /** Returns true if the rectangle intersects with the given rectangle. */
    bool intersect(Rect const * orect) const;

    void update(T min_x, T min_y, T max_x, T max_y);
    /** Moves the rectangle and updates its position. */
    void move(T x, T y);
};

template <typename T>
Rect<T>::Rect(void) {
    update(0, 0, 1, 1);
}

template <typename T>
Rect<T>::Rect(Rect<T> * rect)
    : min_x(rect->min_x), min_y(rect->min_y), max_x(rect->max_x), max_y(rect->max_y) {
    update(rect->min_x, rect->min_y, rect->max_x, rect->max_y);
}

template <typename T>
Rect<T>::Rect(T min_x, T min_y, T max_x, T max_y) {
    update(min_x, min_y, max_x, max_y);
}

template <typename T>
inline T
Rect<T>::width(void) const{
    return max_x - min_x;
}

template <typename T>
inline T
Rect<T>::height(void) const {
    return max_y - min_y;
}

template <typename T>
inline T
Rect<T>::size(void) const {
    return width() * height();
}

template <typename T>
inline bool
Rect<T>::is_inside(Rect const * rect) const {
    return
        min_x >= rect->min_x &&
        max_x <= rect->max_x &&
        min_y >= rect->min_y &&
        max_y <= rect->max_y;
}

template <typename T>
inline void
Rect<T>::update(T min_x, T min_y, T max_x, T max_y) {
    this->min_x = min_x;
    this->min_y = min_y;
    this->max_x = max_x;
    this->max_y = max_y;
    assert(min_x <= max_x);
    assert(min_y <= max_y);
}

template <typename T>
inline bool
Rect<T>::intersect(Rect<T> const * rect) const {
    return
        min_x < rect->max_x && max_x > rect->min_x &&
        min_y < rect->max_y && max_y > rect->min_y;
}

template <typename T>
inline void
Rect<T>::move(T x, T y) {
    update(x, y, x + width(), y + height());
}

#endif /* TEX_RECT_HEADER */
