#ifndef COVERAGESIM_RAW_READ_ID_H
#define COVERAGESIM_RAW_READ_ID_H

#include <cstring>
#include <string>

class RawReadId {
public:
  RawReadId(const std::string& prefix, int width) : prefix(prefix), width(width) {}

  size_t size() const { return prefix.size() + width; }

  void write(char* out, int rank) const {
    std::memcpy(out, prefix.data(), prefix.size());
    char* digit = out + size();
    for (int i = 0; i < width; ++i) {
      *--digit = '0' + rank % 10;
      rank /= 10;
    }
  }

private:
  const std::string prefix;
  const int width;
};

#endif
