#include <string>
#include <fstream>
#include <iostream>
#include <chrono>

#include "database.h"

#include "fastWY.h"

using namespace std;
using uint = unsigned int;

inline void exit_with_help(){
	cout << "-- FastWY method for sequential pattern mining --\n"
			"Usage: train [-options] input_file\n"
			"options:\n"
			"    -m : minimum supportSum (default 1)\n"
			"         If you want to decide on the percentage of the whole, please enter a value less than 1.\n"
			"    -M : minimum length of pattern (default 1)\n"
			"    -L : maximum length of pattern (default 10)\n"
			"    -F : name of reslut file (default output/result.csv)\n"
			"    -p : maximum interval of event (default 0 | -1:none)\n"
			"    -C : whether to do CloSpan (default 0:do not | 1:do)\n"
			"    -R : Repeat resampling in FastWY (default 10000)\n"
			"    -a : significance level (default 0.05)\n"
			"    -s : the Mode of counting supportSum(default 0)\n"
			"     	: 0 is 0 or 1 per record,1 is the number of pattern\n"
			"    -S : double or single-side test (default 0:double | 1:upper-side | 2:lower-side)\n" << endl;

	exit(1);
}

int main(int argc, char **argv){

	//default
	double tMinsup = 1;
	uint tMinpat = 1;
	uint tMaxpat = 10;
	uint tMinItem = 1;
	uint tMaxItem = 10;
	uint tR = 10000;
	double tAlpha = 0.05;
	uint tWild = 0;
	uint tSide = 0;

	int tInterval = 0;
	uint tSupMode = 0;
	uint tCloSpan = 0;
	string tFilename = "result.csv";

	int tI;
	for(tI = 1; tI < argc; tI++){
		if(argv[tI][0] != '-'){
			break;
		}
		if(++tI >= argc){
			exit_with_help();
		}
		switch(argv[tI - 1][1]){
			case 'm':
				tMinsup = atof(argv[tI]);
				break;
			case 'M':
				tMinpat = atoi(argv[tI]);
				break;
			case 'L':
				tMaxpat = atoi(argv[tI]);
				break;
			case 'i':
				tMinItem = atoi(argv[tI]);
				break;
			case 'I':
				tMaxItem = atoi(argv[tI]);
				break;
			case 'F':
				tFilename = argv[tI];
				break;
			case 'p':
				tInterval = atoi(argv[tI]);
				break;
			case 's':
				tSupMode = atoi(argv[tI]);
				break;
			case 'C':
				tCloSpan = atoi(argv[tI]);
				break;
			case 'R':
				tR = atoi(argv[tI]);
				break;
			case 'a':
				tAlpha = atof(argv[tI]);
				break;
			case 'w':
				tWild = atof(argv[tI]);
				break;
			case 'S':
				tSide = atof(argv[tI]);
				break;
			default:
				cerr << "Error unknown option: -" << argv[tI - 1][1] << endl;
				exit_with_help();
				break;
		}
	}

	if(tI >= argc){
		cerr << "Error please input filename" << endl;
		exit_with_help();
	}

	if(tWild == 1 && tInterval != 0){
		cerr << "Error: If you want to consider WildCard, please set Interval value to 0" << endl;
		exit_with_help();
	}

	//read data
	Database tDatabase;
	tDatabase.read(argv[tI]);
	vector<vector<Event> > tTransaction = tDatabase.get_transaction();
	vector<double> tY = tDatabase.get_y();

	chrono::system_clock::time_point start, end;
	start = chrono::system_clock::now();

	//make fastWY
	FastWY tFastWY(tMinsup, tMinpat, tMaxpat, tMinItem, tMaxItem, tR, tAlpha, tCloSpan, tSupMode, tInterval, tWild, tSide);

	tFastWY.init(tTransaction, tY);

	cout << "init finish" << endl;
	tFastWY.main();
	cout << "FastWY finish" << endl;

	end = chrono::system_clock::now();
	double elapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
	elapsed *= 0.001;

	cout << "****************************" << endl;
	cout << "time(sec):" << elapsed << endl;

	//write outcome
	tFastWY.printSigPattern(tFilename);

	tFilename.replace(tFilename.find(".csv"), 4, "_all.csv");
	tFastWY.printTree(tFilename);
	cout << "print finish" << endl;

	return 0;
}
