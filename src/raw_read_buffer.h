#ifndef COVERAGESIM_RAW_READ_BUFFER_H
#define COVERAGESIM_RAW_READ_BUFFER_H

#include <cstddef>
#include <memory>

class RawReadBuffer {
public:
  void resize(size_t length) {
    if (length > capacity) {
      bytes.reset(new unsigned char[length]);
      capacity = length;
    }
    used = length;
  }

  unsigned char* data() { return bytes.get(); }
  const unsigned char* data() const { return bytes.get(); }
  size_t size() const { return used; }

private:
  std::unique_ptr<unsigned char[]> bytes;
  size_t capacity = 0;
  size_t used = 0;
};

#endif
