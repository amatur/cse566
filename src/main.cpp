#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdlib.h>
 #include <stdint.h>
#include <unordered_set>

using namespace std;

typedef unsigned char uchar;

typedef struct {
	string nodes;
	string sequence;
} unitig_struct_t;


double readTimer(){
    return clock() / (double)CLOCKS_PER_SEC;
}

// Get current date/time, format is YYYY-MM-DD HH:mm:ss
inline string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d %X\n", &tstruct);

    return buf;
}

int get_data(const string& unitigFileName,
	uchar*& data,
	vector<unitig_struct_t>& unitigs,
	uint64_t& char_count
	)
{
	try
	{
		ifstream unitigFile;
		unitigFile.open(unitigFileName);
	   	char_count = 0;
	   	string line;

	   	string genome;
	   	// counting total number of reads and their length
	   	while (getline(unitigFile , line))
	   	{
	   		line = line.substr(2,line.length() - 2); // removing the leading "> "
	   		char_count += line.length() + 1; // +1 for the '\0' separator
	   		unitig_struct_t unitig_struct;
	   		unitig_struct.nodes = line;
	   		getline(unitigFile , unitig_struct.sequence); // this is the actual string
	   		unitigs.push_back(unitig_struct);
	   	}

	   	unitigFile.close();

	   	cout << "Input file " << unitigFileName << " contains " << unitigs.size() << " unitigs." << endl;
	   	cout << "The temporary data array will have size " << char_count/(double)1000000 << "MB." << endl;

	   	uint64_t i = 0;
		data = new uchar[char_count];
		for (auto unitig_struct : unitigs)
		{
	    	for (uint64_t j = 0; j < unitig_struct.nodes.length(); j++)
	    	{
	    		data[i] = (uchar)unitig_struct.nodes[j];
	    		i++;

	    	}
	    	data[i] = '\0';
	    	i++;
		}
	   	cout << "Created the data array" << endl;

	   	return EXIT_SUCCESS;
	} catch (exception& error)
	{ // check if there was any error
		std::cerr << "Error: " << error.what() << std::endl;
		return EXIT_FAILURE;
	}

}


void print_maximal_unitigs(const string& unitigFileName)
{
	unordered_set<string> already_printed_unitigs;

	try
	{
		ofstream unitigFile;
		unitigFile.open(unitigFileName);
        unitigFile << "It prints " << endl;
		unitigFile.close();
	}
	catch (exception& error)
	{ // check if there was any error
		std::cerr << "Error: " << error.what() << std::endl;
	}
}

int main(int argc, char** argv)
{
    uint64_t char_count;
    string unitigFileName;

    uchar *data = NULL;


	vector<unitig_struct_t> unitigs;

 	double startTime = readTimer();


	// if (argc <= 1)
	// {
	// 	cerr << "Usage: ./maximality <unitig file> " << endl;
	//  	return 0;
	// }
    //
	// unitigFileName = argv[1];
    //
    //unitigFileName = "test/test3.unitigs.fa";
    unitigFileName = "test/test3.fa";

	cout << unitigFileName << endl;

	if (EXIT_FAILURE == get_data(unitigFileName, data, unitigs, char_count))
	{
		return EXIT_FAILURE;
	}

	// report some information.
	cout << "Time for loading the data: " << readTimer() - startTime << "sec" << endl;

    print_maximal_unitigs("aa");
 	// we sample internal nodes and check their abundances
 	// startTime = readTimer();
 	// print_maximal_unitigs(unitigFileName, unitigs, rlcsa);
 	// cout << "Time for printing maximal unitigs: " << readTimer() - startTime << "sec" << endl;
 	// cout << "Time: " << currentDateTime();

	return EXIT_SUCCESS;
}
