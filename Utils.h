
#ifndef DELAUNAY_UTILS_H
#define DELAUNAY_UTILS_H

#include <algorithm>


template<class container, class element>
bool contains(const container& c, const element& e)
{
  return (std::find(c.begin(), c.end(), e) != c.end());
}


template<class container>
void erase_dupes(container& c)
{
  std::sort(c.begin(), c.end());
  auto last = std::unique(c.begin(), c.end());
  c.erase(last, c.end());
}


#endif // DELAUNAY_UTILS_H
