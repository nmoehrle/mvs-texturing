#pragma once

#include <vector>
#include "math/vector.h"
#include "math/functions.h"
#include "texturing.h"

/** Generates and fills the given vector with distinctive colors for debug purpose. */
void
generate_debug_colors(std::vector<math::Vec4f> & colors);

/** Replaces the encapsuled image of the texture_views with images containing the view id on a distinctive color. */
void
generate_debug_embeddings(std::vector<TextureView> * texture_views);
