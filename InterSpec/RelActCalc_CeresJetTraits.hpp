#ifndef RelActCalc_CeresJetTraits_hpp
#define RelActCalc_CeresJetTraits_hpp

#include "InterSpec_config.h"

#include <limits>

#include <Eigen/Core>
#include <ceres/jet.h>

// Current Ceres provides this complete Eigen trait itself.  Retain the local
// specialization only as a compatibility adapter for older dependencies that
// lack infinity()/quiet_NaN(); otherwise a more-specialized duplicate can
// silently drift from Ceres' generic Jet<T,N> implementation.
#if !defined(INTERSPEC_CERES_JET_NUMTRAITS_HAS_NONFINITE) \
    || !INTERSPEC_CERES_JET_NUMTRAITS_HAS_NONFINITE

namespace Eigen {
  
  /** To use `Jets<>` in Eigen matrices, `ceres` defines
   `template <T,int N> struct NumTraits<ceres::Jet<T, N>>{...}`, however, it doesnt include
   `infinity()` or `quiet_NaN()`, which we need, so need to create an even more specialized version of the
   struct, with these functions.
   */
  template <int N>
  struct NumTraits<ceres::Jet<double, N>> {
    using Real = ceres::Jet<double, N>;
    using NonInteger = ceres::Jet<double, N>;
    using Nested = ceres::Jet<double, N>;
    using Literal = ceres::Jet<double, N>;
    
    static typename ceres::Jet<double, N> dummy_precision() { return ceres::Jet<double, N>(1e-12); }
    static inline Real epsilon() { return Real(std::numeric_limits<double>::epsilon()); }
    static inline Real infinity(){ return Real( std::numeric_limits<double>::infinity() ); }
    static inline Real quiet_NaN(){ return Real( std::numeric_limits<double>::quiet_NaN() ); }
    static inline int digits10() { return std::numeric_limits<double>::digits10; }
    static inline int max_digits10() { return std::numeric_limits<double>::max_digits10; }
    static inline Real highest() { return Real((std::numeric_limits<double>::max)()); }
    static inline Real lowest() { return Real(-(std::numeric_limits<double>::max)()); }
    
    enum {
      IsComplex = 0,
      IsInteger = 0,
      IsSigned,
      ReadCost = 1,
      AddCost = 1,
      // For Jet types, multiplication is more expensive than addition.
      MulCost = 3,
      HasFloatingPoint = 1,
      RequireInitialization = 1
    };
    
    template <bool Vectorized>
    struct Div {
      enum {
#if defined(EIGEN_VECTORIZE_AVX)
        AVX = true,
#else
        AVX = false,
#endif
        
        // Assuming that for Jets, division is as expensive as
        // multiplication.
        Cost = 3
      };
    };
  };//struct NumTraits<ceres::Jet<double, N>>
}//namespace Eigen

#endif // old Ceres Eigen::NumTraits<Jet> compatibility

#endif // RelActCalc_CeresJetTraits_hpp
