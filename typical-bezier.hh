#pragma once

#include <geometry.hh>

// `a` should contain exactly 3 points
Geometry::Point2DVector typicalBezier(const Geometry::Point2DVector &a, size_t degree);
