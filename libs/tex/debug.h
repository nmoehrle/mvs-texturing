#pragma once

#include <vector>
#include "texturing.h"

/** Replaces the encapsuled image of the texture_views with images containing the view id on a distinctive color. */
void
generate_debug_embeddings(std::vector<TextureView> * texture_views);
