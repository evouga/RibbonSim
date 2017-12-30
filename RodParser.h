#ifndef RODPARSER_H
#define RODPARSER_H

#include "RodConfig.h"

/*
 * Reads a weave definition from file. The file must contain three section: a header, list of rods, and list of constraints.
 *
 * The header consists of:
 * (int)     number of rods
 * (int)     number of constraints
 * (double)  rod thickness (in meters)
 * (double)  Young's modulus (in Newtons per square meter)
 * (double)  rod density (in kilograms per cubic meter)
 *
 * Each rod consists of:
 * (int)                   number of vertices NV
 * (0 or 1)                whether or not the rod is closed. If the rod is closed, the number of segments NS is equal to NV in the spec below. Otherwise, NS = NV-1.
 * (list of 3*NV doubles)  centerline vertex positions
 * (list of 3*NS doubles)  orientation of the centerline segments (segment i connects vertex i and i+1). This vector specifies the "thickness" direction of the rod. It does not need to be normalized but must be perpendicular to its segment.
 * (list of NS double)     width of each rod segment (in meters)
 *
 * Each constraint consists of
 * (int)     first rod coupled by the constraint (zero-indexed)
 * (int)     second rod coupled by the constraint (zero-indexed). This can be the same as the first rod.
 * (int)     segment of first rod (zero-indexed) on which the constraint lies
 * (int)     ditto, for the second rod
 * (double)  location of the constraint on the first rod's segment: if the segment connects verties i and i+1, the constrained point is at (1-alpha) v_i + alpha v_{i+1}
 * (double)  ditto, for the second rod segment
 * (double)  stiffness (strength) of the constraint (N/m^2).
 */
RodConfig *readRod(const char *filename);

#endif
