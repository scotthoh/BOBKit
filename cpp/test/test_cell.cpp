#include <clipper/core/cell.h>
#include <iostream>

using namespace clipper;

int main(int argc, char **argv)
{
  Cell_descr cdes(100.0, 100.0, 100.0);
  Cell cell(cdes);
  std::cout << cell.debug() << std::endl;
  std::cout << cell.volume() << std::endl;
  std::cout << cell.format() << std::endl;

  return 0;
}