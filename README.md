Description
==========

Parallel-DivSufSort is a parallel lightweight suffix array construction algorithm for byte alphabets written in C++.
It is a implementation based on:
+ [divsufsort](https://github.com/y-256/libdivsufsort) implementation of induced sorting by Yuta Mori.
+ [parallel-range-light](https://github.com/jlabeit/parallel-range-light) implementation of prefix doubling for integer alphabets. 

A detailed description and benchmarks of the algorithm can be found in the following work.
> Julian Labeit, Julian Shun, and Guy E. Blelloch. Parallel Lightweight Wavelet Tree, Suffix Array and FM-Index Construction. DCC 2015.

Installation
==========
The following steps have been tested on Ubuntu 14.04 with gcc 5.3.0 and cmake 2.8.12.
```shell
git clone https://github.com/jlabeit/parallel-divsufsort.git
cd parallel-divsufsort
mkdir build
cd build
cmake ..
make
make install
```
Note that in the default version the cilkplus implementation by gcc is used for parallelization.
To change this setting edit parallelization settings in the CMakeLists.txt file.

Getting Started
=========
An example application can be found in examples/main.cpp.
The library provides two basic functions to build the suffix array over a text.

```c++
// 32 bit version.
uint8_t divsufsort(const uint8_t *T, int32_t *SA, int32_t n);
// 64 bit version.
uint8_t divsufsort(const uint8_t *T, int64_t *SA, int64_t n);
```

To use the library include the header `divsufsort.h`, link against the library `divsufsort` and `libprange`.

