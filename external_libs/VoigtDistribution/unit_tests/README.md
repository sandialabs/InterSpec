# Unit Tests for VoigtDistribution

This directory contains unit tests for the VoigtDistribution library.

## Test Files

- `test_faddeeva.cpp` - Tests for Faddeeva function implementations
- `test_voigt_functions.cpp` - Tests for Voigt profile functions
- `test_distribution_combinations.cpp` - Tests for combined Voigt + exponential tail distributions

## Running Tests

```bash
cd build
make test_voigt_distribution
./test_voigt_distribution
```

## Test Coverage

- Faddeeva function accuracy against reference implementations
- Voigt profile normalization and edge cases
- Exponential tail parameter validation
- Combined distribution integration
- Automatic differentiation compatibility