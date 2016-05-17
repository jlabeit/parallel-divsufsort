#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <divsufsort.h>


using namespace std;

typedef int32_t num_type; // Currently only int32_t and int64_t are supported.
int times = 3; // How often the time measurement is repeated.

int main(int argc, char* args[]) {
	if (argc != 2) {
		cout << "Expected one argument (input file)."
			<< endl;	
		return -1;
	}

	string text;
	{ // Read input file.
		ifstream input_file(args[1]);
		input_file.seekg(0, ios::end);   
		text.reserve(input_file.tellg());
		input_file.seekg(0, ios::beg);
		text.assign((istreambuf_iterator<char>(input_file)),
				istreambuf_iterator<char>());
	}
	//text += '\0';
	num_type size = text.size();
	num_type *SA = new num_type[size];
	for (int i = 0; i < times; ++i) {
		auto start = chrono::steady_clock::now();
		divsufsort((sauchar_t*)text.data(), SA, size);
		auto end = chrono::steady_clock::now();
		auto diff = end - start;
		cout <<	chrono::duration <double, milli> (diff).count();
		cout.flush();
		if (i < times -1) cout << ", ";
	}
	cout << endl;
	if (sufcheck((sauchar_t*)text.data(), SA, size, false)) {
		cout << "Sufcheck failed!" << endl;
		return -1;
	}
	return 0;
}
