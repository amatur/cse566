#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <unordered_set>
#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <list>
#include <stack>
#include <unordered_map>
using namespace std;

typedef unsigned char uchar;
typedef struct {
    int serial;
    int startPos;
    int endPos;
} new_node_info_t;
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
typedef struct {
    //1 means +, 0 means -
    bool left;
    bool right;
    int toNode;
    int kmerStartIndex;
    int kmerEndIndex;
} newEdge_t;


// ------- PARAMETERS -------- //
int k = 21;
string  unitigFileName  = "data/list_reads.unitigs.fa";


vector<vector<edge_t> > adjList;
vector<vector<edge_t> > newAdjList;
vector<edge_t> resolveLaterEdges;
vector<unitig_struct_t> unitigs;


inline string plus_strings(const string& a, const string& b, size_t kmersize) {
    if (a == "") return b;
    if (b == "") return a;
    string ret = a + b.substr(kmersize - 1, b.length() - (kmersize - 1));
    return ret;
}

string delSpaces(string &str) {
    str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
    return str;
}

bool charToBool(char c) {
    if (c == '+') {
        return true;
    } else {
        if (c != '-') cout << "ERRRRRROOR!" << endl;
        return false;
    }
}

string reverseComplement(string base) {
    int len = base.length();
    char* out = new char[len + 1];
    out[len] = '\0';
    for (int i = 0; i < len; i++) {
        if (base[i] == 'A') out[len - i - 1] = 'T';
        else if (base[i] == 'C') out[len - i - 1] = 'G';
        else if (base[i] == 'G') out[len - i - 1] = 'C';
        else if (base[i] == 'T') out[len - i - 1] = 'A';
    }
    string outString(out);
    free(out);
    return outString;
}

double readTimer() {
    return clock() / (double) CLOCKS_PER_SEC;
}

inline string currentDateTime() {
    // Get current date/time, format is YYYY-MM-DD HH:mm:ss
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof (buf), "%Y-%m-%d %X\n", &tstruct);
    return buf;
}

int countInArcs(int node) {
    //system("grep -o -i :2: unitigs.fa | wc -l > incount.txt");
    int count;
    string line;
    string countFile = "incount.txt";
    //string inputFile = "unitigs.fa";
    string inputFile = unitigFileName;

    ostringstream stringStream;
    stringStream << "grep -o -i :" << node << ": " << inputFile << " | wc -l > " << countFile;
    string copyOfStr = stringStream.str();
    system(copyOfStr.c_str());
    ifstream cf;
    cf.open(countFile);
    getline(cf, line);
    line = delSpaces(line);
    sscanf(line.c_str(), "%d", &count);
    cf.close();
    return count;
}

int countOutArcs(int node) {
    return (adjList.at(node)).size();
}

inline char boolToCharSign(bool sign){
   return (sign == true) ? '+' : '-';
}

void printGraph(vector<vector<edge_t> > adjList){
    for (int i = 0; i < adjList.size(); i++) {
        cout<<i<<"# ";
        for(edge_t edge : adjList.at(i)){
            cout<<boolToCharSign(edge.left)<<":"<<edge.toNode<<":"<<boolToCharSign(edge.right)<<", ";
        }
        cout<<endl;
    }
}


class Graph{
public:
    int V = adjList.size();
    char* color;
    int* p;
    bool* nodeSign;
    new_node_info_t* oldToNew;
    int time;

    Graph(){
        color = new char[V];
        p = new int[V];
        nodeSign = new bool[V];
        oldToNew = new new_node_info_t[V];
        for(int i = 0; i<V ; i++){
            oldToNew[i].serial = -1;
        }
    }

    //grey hole create representative
    // ta na hole get representative



    // The function to check if edge u-v can be considered as next edge in
    // Euler Tout
//    bool isValidNextEdge(int u, int v)
//    {
//      // The edge u-v is valid in one of the following two cases:
//
//      // 1) If v is the only adjacent vertex of u
//        return edgeFrom[u] = adjX.at(1).to;
//   }


    void DFS_visit(int u){
        stack<edge_t> s;

        edge_t uEdge;
        uEdge.toNode = u;
        s.push(uEdge);

        while(!s.empty()){
            edge_t xEdge = s.top();
            int x = xEdge.toNode;
            s.pop();
            if(color[x] == 'w'){
                time = time + 1;
                color[x] = 'g';
                s.push(xEdge);

                vector<edge_t> adjx = adjList.at(x);
                int t = x;
                if(adjx.size() == 0){
                    while(p[t]!=-1){
                         cout<<p[t]<<" ";
                         t = p[t];
                    }

                }
                cout<<endl;


                for(edge_t edge: adjx ){
                    int y = edge.toNode;
                    edge_t yEdge;
                    yEdge.toNode = y;

                    if(color[y] == 'w'){
                        //adding modification
                        //check if it is a valid walk
                        //valid if
                        if(p[x] == -1){
                            nodeSign[x] = edge.left;
                            nodeSign[y] = edge.right;
                            p[y] = x;
                            s.push(yEdge);
                        }else if(nodeSign[x] == edge.left){
                            nodeSign[y] = edge.right;
                            p[y] = x;
                            s.push(yEdge);
                        }


                    }
                }

            }else if(color[x] == 'g'){
                time = time + 1;
                color[x] = 'b';
            }

        }
    }

    void DFS(){

         for(int i=0; i<V; i++){
             color[i] = 'w';
             p[i] = -1;
         }

         for(int i=0; i<V; i++){
             if(color[i] == 'w'){
                 DFS_visit(i);
             }
         }

    }

    ~Graph(){
        delete [] color;
        delete [] p;
        delete [] nodeSign;
        delete [] oldToNew;
    }
};




int getNodeNumFromFile(string filename){
    //grep '>' list_reads.unitigs.fa | tail -n 1
    int nodeCount;
    string countFile = "incount.txt";
    ostringstream stringStream;
    stringStream << "grep '>' "<< filename << " | tail -n 1" << countFile;
    string copyOfStr = stringStream.str();
    system(copyOfStr.c_str());
    ifstream cf;
    cf.open(countFile);
    string line;
    getline(cf, line);
    sscanf(line.c_str(), "%*c %d", &nodeCount);
    cf.close();
    return nodeCount;
}

int get_data(const string& unitigFileName,
        uchar*& data,
        vector<unitig_struct_t>& unitigs,
        uint64_t& char_count
        ) {
    ifstream unitigFile;
    try {
        unitigFile.open(unitigFileName);
        //            if(unitigFile==NULL){
        //                throw "ERROR: File does not exist!!!";
        //            }else{
        //                cout<<"File opened successfully."<<endl;
        //            }

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
    bool doCont = false;

    getline(unitigFile, line);


    do {
        cout << line << endl;
        sscanf(line.c_str(), "%*c %d %s  %s  %s %[^\n]s", &nodeNum, lnline, kcline, kmline, edgesline);
        //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-

        unitig_struct_t unitig_struct;

        sscanf(lnline, "%*5c %d", &unitig_struct.ln);
        sscanf(kcline, "%*5c %d", &unitig_struct.kc);
        sscanf(kmline, "%*5c %f", &unitig_struct.km);

        char c1, c2;
        stringstream ss(edgesline);

        vector<edge_t> edges;
        while (getline(ss, line, ' ')) {
            if (delSpaces(line).length() != 0) {
                sscanf(line.c_str(), "%*2c %c %*c %d  %*c  %c", &c1, &nodeNum, &c2); //L:-:0:-
                edge_t newEdge;

                newEdge.left = charToBool(c1);
                newEdge.right = charToBool(c2);
                newEdge.toNode = nodeNum;
                edges.push_back(newEdge);
            }

        }
        adjList.push_back(edges);

        // for (edge_t e: edges){
        //     cout<<int(e.left)<<","<<int(e.right)<< " "<<e.toNode<<endl;
        // }

        unitig_struct.serial = nodeNum;

        doCont = false;
        while (getline(unitigFile, line)) {
            if (line.substr(0, 1).compare(">")) {
                unitig_struct.sequence = unitig_struct.sequence + line;
                unitigs.push_back(unitig_struct);
            } else {
                doCont = true;
                break;
            }
        }
    } while (doCont);
    unitigFile.close();

    cout << "Complete reading input." << endl;
    return EXIT_SUCCESS;
}


int main(int argc, char** argv) {
    uint64_t char_count;
    uchar *data = NULL;

    double startTime = readTimer();

    cout << "Starting reading " << unitigFileName << endl;
    if (EXIT_FAILURE == get_data(unitigFileName, data, unitigs, char_count)) {
        return EXIT_FAILURE;
    }
    
    Graph G;
    cout<<G.V;
    
    printGraph(adjList);
    


    return EXIT_SUCCESS;
}
