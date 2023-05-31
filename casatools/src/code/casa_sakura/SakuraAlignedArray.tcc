#include <cassert>
#include <iostream>

#include <libsakura/sakura.h>

#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogOrigin.h>

namespace casa {
template<typename T>
inline SakuraAlignedArray<T>::SakuraAlignedArray(size_t num_data) :
  num_data_(num_data) {
#if !defined(NDEBUG)
  casacore::LogIO logger(casacore::LogOrigin("SakuraAlignedArray", "SakuraAlignedArray", WHERE));
  logger << casacore::LogIO::DEBUGGING << "Constructing SakuraAlignedArray..." << casacore::LogIO::POST;
#endif

  initialize();

#if !defined(NDEBUG)
  logger << casacore::LogIO::DEBUGGING << "  Initial Address = " << storage_ << casacore::LogIO::POST;
  logger << casacore::LogIO::DEBUGGING << "  Aligned Address = " << data_ << casacore::LogIO::POST;
#endif
}

template<typename T>
inline SakuraAlignedArray<T>::SakuraAlignedArray(casacore::Vector<T> const &in_vector) :
  num_data_(in_vector.nelements()) {
#if !defined(NDEBUG)
  casacore::LogIO logger(casacore::LogOrigin("SakuraAlignedArray", "SakuraAlignedArray", WHERE));
  logger << casacore::LogIO::DEBUGGING << "Constructing SakuraAlignedArray..." << casacore::LogIO::POST;
#endif

  initialize();

  T *ptr = data_;
  for (size_t i = 0; i < num_data_; ++i) {
    ptr[i] = in_vector(i);
  }

#if !defined(NDEBUG)
  logger << casacore::LogIO::DEBUGGING << "  Initial Address = " << storage_ << casacore::LogIO::POST;
  logger << casacore::LogIO::DEBUGGING << "  Aligned Address = " << data_ << casacore::LogIO::POST;
#endif
}

template<typename T>
inline void SakuraAlignedArray<T>::initialize() {
  storage_ = nullptr;
  data_ = nullptr;

  size_t size_required = sizeof(T) * num_data_;
  size_t size_of_arena = size_required + LIBSAKURA_SYMBOL(GetAlignment)() - 1;
  storage_ = malloc(size_of_arena);
  if (storage_ == nullptr) {
    data_ = nullptr;
    throw std::bad_alloc();
  }
  data_ = reinterpret_cast<T *>(LIBSAKURA_SYMBOL(AlignAny)(
				size_of_arena, storage_, size_required));
  assert(data_ != nullptr);
  assert(LIBSAKURA_SYMBOL(IsAligned)(data_));
}

template<typename T>
inline SakuraAlignedArray<T>::~SakuraAlignedArray() {
#if !defined(NDEBUG)
  casacore::LogIO logger(casacore::LogOrigin("SakuraAlignedArray", "~SakuraAlignedArray", WHERE));
  logger << casacore::LogIO::DEBUGGING << "Destructing SakuraAlignedArray..." << casacore::LogIO::POST;
#endif

  free(storage_);
}

}  // End of casa namespace.
