/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_DEBUG_HEADER
#define TEX_DEBUG_HEADER

#include <vector>
#include "texturing.h"

/** Replaces the encapsuled image of the texture_views with images containing the view id on a distinctive color. */
void
generate_debug_embeddings(std::vector<TextureView> * texture_views);

#endif /* TEX_DEBUG_HEADER */
