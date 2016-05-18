Specifications
==========

  * Processors:       4 X 12 X 1.9Ghz AMD Opteron(tm) Processor 6168 (48 Processors)
  * L1 Cache:         64 KB
  * L2 Cache:         512 KB 
  * RAM:              2 Gb main memory
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
Remarks: For all but the input file 'english' the full version of the text file was used. 
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

Runnin times here measured in seconds.

|  Files         | File Size |  divsufsort  |
|:---------------|----------:|-------------:|
|  dblp.xml      |   295     |  1413        |
|  dna           |   385     |  1927        |
|  english       |  1024     |  5121        |
|  pitches       |    53     |   267        |
|  proteins      |  1129     |  5647        |
|  sources       |   201     |  1006        |

Memory consumption is stated in MB.
