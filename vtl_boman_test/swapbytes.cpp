#include "swapbytes.h"

bool isCpuLittleEndian()
{
    static const int __one__ = 1;
    static const bool ret = 1 == *(char *)(&__one__); // CPU endianness
    return ret;
}
