#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <divsufsort.h>


using namespace std;

void printPos(int pos, string& text) {
	for (int i = 0; i < 10 && i + pos < (int)text.size(); ++i) {
		cout << text[i+pos];
	} cout << endl;
}

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
	int size = text.size();
	int *SA = new int[size];
	auto start = chrono::steady_clock::now();
	divsufsort((sauchar_t*)text.data(), SA, size);
	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	cout << "Parallel Divsufsort time: " <<
		chrono::duration <double, milli> (diff).count()<< " ms" << endl;
	return 0;

}
