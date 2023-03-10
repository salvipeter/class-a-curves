#include "typical-bezier.hh"

using namespace Geometry;

int main(int argc, char **argv) {
  Point2DVector a = { { 2, 3 }, { 5, 5 }, { 7, 2 } };
  auto cpts = typicalBezier(a, 3);
  for (auto p : cpts)
    std::cout << p << std::endl;
}
