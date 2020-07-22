#ifndef __COMMON_H__
#define __COMMON_H__

#include <iostream>
#include <string>

class Common {
  public:

    static bool isNucl(char c);
    static bool isAA(char c);
    static std::string reverseComplement(const std::string seq);
    static char revComp(const char* c);

};

#endif
