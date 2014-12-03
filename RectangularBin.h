#pragma once

#include "Rect.h"
#include <list>
#include <cmath>

/**
  * Implementation of the binpacking algorithm GUILLUTINE from
  * <A HREF="http://clb.demon.fi/files/RectangleBinPack.pdf">A Thousand Ways to Pack the Bin - A Practical Approach to Two-Dimensional Rectangle Bin Packing</A>
  */
class RectangularBin {
    private:
        int width;
        int height;
        std::list<Rect<int> > rects;

    public:
        /**
          * Initializes the rectangular binpacking algorithm to fill a rectangle of the given size.
          */
        RectangularBin(int width, int height);

        /** Returns true and changes the position of the given rect if it fits into the bin. */
        bool insert(Rect<int> * rect);
};
