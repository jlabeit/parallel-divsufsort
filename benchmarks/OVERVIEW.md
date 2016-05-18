Specifications
==========

  * Processors:       4 X 12 X 1.9Ghz AMD Opteron(tm) Processor 6168 (48 Processors)
  * L1 Cache:         64 KB
  * L2 Cache:         512 KB 
  * RAM:              256 Gb main memory
  * Operating system: Ubuntu 12.04.5 LTS
  * Compiler:         GCC version 5.2.0
  * Opt-Flags:       '-O2 -ffast-math -funroll-loops'

The times are the average of three runs, in seconds.
Memory usage is the high water usage measured by getrusage syscall in mega bytes.
For details see 'examples/timings.cpp'.
For parallelization the gcc cilkplus support was used.
Runs with more than 8 active threads were executed with interleaved memory policy (numactl -i all).

Test Files
==========

The text collection from the [Pizza & Chilli](http://pizzachili.dcc.uchile.cl) corpus was used for the benchmarks.

*Remarks:* For all but the input file 'english' the full version of the text file was used. 
For the input 'english' the 1024Mb version was used.

Results
==========

The following tables compare running times and memory consumption of the original divsufsort with parallel-divsufsort with different thread count.

|  Files         | File Size |  divsufsort  |
|:---------------|----------:|-------------:|
|  dblp.xml      |   295     |     65.1     |
|  dna           |   385     |    145.0     | 
|  english       |  1024     |    417.5     | 
|  pitches       |    53     |      9.1     | 
|  proteins      |  1129     |    457.6     | 
|  sources       |   201     |     40.4     | 

Running times here stated in seconds.

|  Files         | File Size |  divsufsort  |  1   |  2   |  4   |  8   | 16   | 32   | 48   |
|:---------------|----------:|-------------:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
|  dblp.xml      |   295     |  1413        | 1853 | 1859 | 1870 | 1891 | 1915 | 1926 | 1939 |
|  dna           |   385     |  1927        | 2504 | 2508 | 2511 | 2519 | 2526 | 2545 | 2555 |
|  english       |  1024     |  5121        | 6878 | 6883 | 6892 | 6904 | 6921 | 6930 | 6956 |
|  pitches       |    53     |   267        |  350 |  351 |  351 |  352 |  354 |  359 |  363 |
|  proteins      |  1129     |  5647        | 7634 | 7635 | 7636 | 7638 | 7641 | 7649 | 7656 |
|  sources       |   201     |  1006        | 1314 | 1315 | 1317 | 1319 | 1323 | 1331 | 1339 |

Memory consumption stated in mega bytes.
