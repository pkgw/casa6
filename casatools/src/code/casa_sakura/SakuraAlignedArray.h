#ifndef _CASA_SAKURA_ALIGNED_ARRAY_H_
#define _CASA_SAKURA_ALIGNED_ARRAY_H_

#include <iostream>
#include <memory>
#include <string>

//#include <libsakura/sakura.h>

#include <casacore/casa/aipstype.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>

namespace casa { //# NAMESPACE CASA - BEGIN

template<typename T>
class SakuraAlignedArray {
public:
  SakuraAlignedArray(size_t num_data);
  SakuraAlignedArray(casacore::Vector<T> const &in_vector);
  ~SakuraAlignedArray();

  T *data() const {return data_;}
  casacore::Vector<T> casaVector() const {
      return casacore::Vector<T>(casacore::IPosition(1, num_data_), data_,
                                 casacore::SHARE);
  }
private:
  void initialize();
  size_t num_data_;      // number of data to be stored
  void *storage_;        // starting address of allocated memory (unaligned)
  T *data_;               // pointer to aligned data
};

} //# NAMESPACE CASA - END

#include <casa_sakura/SakuraAlignedArray.tcc>

#endif /* _CASA_SAKURA_ALIGNED_ARRAY_H_ */
