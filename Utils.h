
#ifndef INTERP_UTILS_H
#define INTERP_UTILS_H


#include <algorithm>

template<class container, class element>
bool contains(const container& c, const element& e)
{
  return (std::find(c.begin(), c.end(), e) != c.end());
}


#endif // INTERP_UTILS_H
