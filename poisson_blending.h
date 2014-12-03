#include "mve/image.h"

void
poisson_blend(mve::ByteImage::ConstPtr src, mve::ByteImage::ConstPtr mask,
    mve::ByteImage::Ptr dest, float alpha);
