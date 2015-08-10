#include "mve/image.h"

void
poisson_blend(mve::FloatImage::ConstPtr src, mve::ByteImage::ConstPtr mask,
    mve::FloatImage::Ptr dest, float alpha);
