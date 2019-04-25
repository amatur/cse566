#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdlib.h>
 #include <stdint.h>
#include <unordered_set>
#include<map>
#include<set>
#include <sstream>
#include<algorithm>
#include <cstdlib>



using namespace std;

typedef unsigned char uchar;

typedef struct {
	int serial;
	string sequence;
    int ln;
    int kc;
    float km;
} unitig_struct_t;


typedef struct {
    //1 means +, 0 means -
	bool left;
    bool right;
    int toNode;
} edge_t;

vector<int> orphanedLink;
vector<int> branchitigs;
vector<vector<edge_t> > adjList;


inline string plus_strings(const string& a, const string& b, size_t kmersize)
{
	if (a == "") return b;
	if (b == "") return a;
	string ret = a + b.substr(kmersize - 1, b.length() - (kmersize - 1));
	return ret;
}

string delSpaces(string &str)
{
   str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
   return str;
}

bool charToBool(char c)
{
    if(c=='+'){
        return true;
    }else{
        if(c!='-') cout<<"ERRRRRROOR!"<<endl;
        return false;
    }
}

void printGraph(map<string, set<string> > adjList){

        for (map<string, set<string> > ::iterator it = adjList.begin(); it != adjList.end(); it++ )
        {
            //cout<<it->first<<" :::: ";
            set<string> v = it->second;

            if(v.empty()){
                cout << '(' << it->first <<", " << " --- " << ')'<<endl;
            }
            for (set<string> :: iterator i = v.begin(); i != v.end(); ++i){
                cout << '(' << it->first <<", " << *i << ')'<<endl;
            }
            //cout<<endl;
        }
}


string reverseComplement(string base){
    int len =  base.length();
    char* out = new char[len+1];
    out[len] = '\0';
    for (int i = 0; i < len; i++) {
        if (base[i]=='A') out[len-i-1]='T';
        else if (base[i]=='C') out[len-i-1]='G';
        else if (base[i]=='G') out[len-i-1]='C';
        else if (base[i]=='T') out[len-i-1]='A';
    }
    string outString(out);
    free(out);
    return outString;
}


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
map<string, set<string> > makeDBGfromKmers(set<string> kmer_set, int k){
    map<string, set<string> > adjList; //simple cycle

    for(set<string>:: iterator it=kmer_set.begin(); it!=kmer_set.end(); ++it){
        string kmer = *it;
        //cout<<"main: "<<k1mer<<" ";

        string kmerFrom, kmerTo;

        kmerFrom = kmer.substr(0, k-1);
        kmerTo = kmer.substr(1, k-1);

        if(adjList.count(kmerFrom)==0){
            adjList[kmerFrom] = set<string>();
        }
        adjList[kmerFrom].insert(kmerTo);


        //added
        if(adjList.count(kmerTo)==0){
            adjList[kmerTo] = set<string>();
        }
    }

    return adjList;
}


int countInArcs(int node){
    int count;
    string line;
    string countFile = "incount.txt";
    string inputFile = "unitigs.fa";
    ostringstream stringStream;
    stringStream << "grep -o -i :"<<node<<": "<<inputFile<<" | wc -l > "<<countFile;
    string copyOfStr = stringStream.str();
    system(copyOfStr.c_str());
    ifstream cf;
    cf.open(countFile);
    getline(cf , line);
    line = delSpaces(line);
    sscanf (line.c_str(),"%d",  &count);
    cf.close();
    return count;
    //system("grep -o -i :2: unitigs.fa | wc -l > incount.txt");
}

int countOutArcs(int node){
    return (adjList.at(node)).size();
}
void computeYtoV(){



}



int get_data(const string& unitigFileName,
	uchar*& data,
	vector<unitig_struct_t>& unitigs,
	uint64_t& char_count
	)
{
        ifstream unitigFile;
        try
    	{
    		unitigFile.open(unitigFileName);

            if(unitigFile==NULL){
                throw "ERROR: File does not exist!!!";
            }else{
                cout<<"File opened successfully."<<endl;
            }

        } catch (const char* msg) {
            cout << msg << endl;
        }

	   	string line;

        int nodeNum;
        char lnline[20];
        char kcline[20];
        char kmline[20];
        char edgesline[100000];
        char c;
        int kc;
        int km;

	   	// counting total number of reads and their length
	   	while (getline(unitigFile , line))
	   	{
	   	    sscanf (line.c_str(),"%*c %d %s  %s  %s %[^\n]s",  &nodeNum, lnline, kcline, kmline, edgesline);
            //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-

            unitig_struct_t unitig_struct;

            sscanf (lnline,"%*5c %d", &unitig_struct.ln);
            sscanf (kcline,"%*5c %d", &unitig_struct.kc);
            sscanf (kmline,"%*5c %f", &unitig_struct.km);

            char c1, c2;
            stringstream ss(edgesline);

            vector<edge_t> edges;
            while(getline(ss, line, ' ')) {
                if(delSpaces(line).length()!=0){
                    sscanf (line.c_str(),"%*2c %c %*c %d  %*c  %c", &c1,  &nodeNum, &c2);
                    edge_t newEdge;
                    //L:-:0:-
                    //cout<<nodeNum<<" "<<line<<endl;
                    newEdge.left = charToBool(c1);
                    newEdge.right = charToBool(c2);
                    newEdge.toNode = nodeNum;
                    edges.push_back(newEdge);
                }

            }
            adjList.push_back(edges);
            cout<<((adjList.at(0)).size())<<endl;
            // for (edge_t e: edges){
            //     cout<<int(e.left)<<","<<int(e.right)<< " "<<e.toNode<<endl;
            // }

	   		unitig_struct.serial = nodeNum;
	   		getline(unitigFile , unitig_struct.sequence); // this is the actual string
	   		unitigs.push_back(unitig_struct);
	   	}

	   	unitigFile.close();

	   	cout << "Created the data array" << endl;

	   	return EXIT_SUCCESS;


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

    //grep -o -i iphone Tweet_Data | wc -l


    //system("grep -rIhEo \"\\b[a-zA-Z0-9.-]+\\@[a-zA-Z0-9.-]+\\/[a-zA-Z0-9.-]+\\b\" /home/*.txt > email.txt");

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
    unitigFileName = "unitigs.fa";


	if (EXIT_FAILURE == get_data(unitigFileName, data, unitigs, char_count))
	{
		return EXIT_FAILURE;
	}


    //cout<<"HHHH "<<countInArcs(0)<<endl;
    cout<<"Out arc of 0 "<<countOutArcs(2)<<endl;
	cout << unitigFileName << endl;


	// report some information.
	//cout << "Time for loading the data: " << readTimer() - startTime << "sec" << endl;

    //print_maximal_unitigs("aa");


 	// we sample internal nodes and check their abundances
 	// startTime = readTimer();
 	// print_maximal_unitigs(unitigFileName, unitigs, rlcsa);
 	// cout << "Time for printing maximal unitigs: " << readTimer() - startTime << "sec" << endl;
 	// cout << "Time: " << currentDateTime();

	return EXIT_SUCCESS;
}
