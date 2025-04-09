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