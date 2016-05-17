#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <unistd.h>
#include <sys/resource.h>
#include <stdio.h>

#include <divsufsort.h>


using namespace std;

typedef int32_t num_type; // Currently only int32_t and int64_t are supported.
int times = 3; // How often the time measurement is repeated.

size_t getPeakRSS() {
	struct rusage rusage;
	getrusage( RUSAGE_SELF, &rusage);
	return (size_t)(rusage.ru_maxrss * 1024L);
}

int main(int argc, char* args[]) {
	if (argc != 2) {
		cout << "Expected one argument (input file)."
			<< endl;	
		return -1;
	}
	std::cout.precision(4);
	string text;
	{ // Read input file.
		ifstream input_file(args[1]);
		input_file.seekg(0, ios::end);   
		text.reserve(input_file.tellg());
		input_file.seekg(0, ios::beg);
		text.assign((istreambuf_iterator<char>(input_file)),
				istreambuf_iterator<char>());
	}
	num_type size = text.size();
	num_type *SA = new num_type[size];
	for (int i = 0; i < times; ++i) {
		auto start = chrono::steady_clock::now();
		divsufsort((sauchar_t*)text.data(), SA, size);
		auto end = chrono::steady_clock::now();
		auto diff = end - start;
		cout <<	chrono::duration <double, milli> (diff).count() / 1000.0 << ", ";
		cout.flush();
	}
	cout << getPeakRSS() / (1024*1024)<< endl;
	if (sufcheck((sauchar_t*)text.data(), SA, size, false)) {
		cout << "Sufcheck failed!" << endl;
		return -1;
	}
	return 0;
}
