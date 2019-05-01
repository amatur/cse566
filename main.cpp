// --- VERSION 1.4 ----
// Sequence stitching done, not sure if it is okay, works for toy example
// Statistics with correct edge sign label, sequence stitching still not done, a lot of bugs fixed

#include<cmath>
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
    edge_t edge;
    int fromNode;
} edge_both_t;

typedef struct {
    edge_t edge;
    int kmerStartIndex;
    int kmerEndIndex;
} newEdge_t;

int DEBUG = 2;

// ------- PARAMETERS -------- //
int K = 5;
string UNITIG_FILE = "data/list_reads.unitigs.fa";


vector<vector<edge_t> > adjList;
vector<vector<newEdge_t> > newAdjList;
vector<edge_both_t> resolveLaterEdges;
vector<unitig_struct_t> unitigs;
vector<string> newSequences;

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
    string inputFile = UNITIG_FILE;
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

int maximumUnitigLength() {
    //grep '>' list_reads.unitigs.fa | cut -d : -f3 | cut -d ' ' -f1 | sort -n | tail -n 1 > incount.txt
    int count;
    string line;
    string countFile = "incount.txt";
    ostringstream stringStream;
    stringStream << "grep '>' " << UNITIG_FILE << " | cut -d : -f3 | cut -d ' ' -f1 | sort -n | tail -n 1 > " << countFile;
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

inline char boolToCharSign(bool sign) {
    return (sign == true) ? '+' : '-';
}


// @@ --- ALL PRINTING CODE --- //

void printBCALMGraph(vector<vector<edge_t> > adjList) {
    for (int i = 0; i < adjList.size(); i++) {
        cout << i << "# ";
        for (edge_t edge : adjList.at(i)) {
            cout << boolToCharSign(edge.left) << ":" << edge.toNode << ":" << boolToCharSign(edge.right) << ", ";
        }
        cout << endl;
    }
}

void printAllBCALMSequences(vector<unitig_struct_t> unitigs) {
    for (unitig_struct_t unitig : unitigs) {
        cout << unitig.serial << ": " << unitig.ln << " " << unitig.sequence.length() << endl; //sequence only^
        // full print
        //cout<<unitig.serial<<": "<<unitig.ln<<" "<<unitig.sequence.length()<<":"<<unitig.sequence<<endl;
    }
}


class Graph {
public:
    int V = adjList.size();
    int countNewNode;
    int time;

    char* color;
    int* p;
    bool* nodeSign;
    new_node_info_t* oldToNew;
    bool* saturated;

    Graph() {
        color = new char[V];
        p = new int[V];
        nodeSign = new bool[V];
        oldToNew = new new_node_info_t[V];
        saturated = new bool[V];
        for (int i = 0; i < V; i++) {
            oldToNew[i].serial = -1;
            saturated[i] = false;
        }
    }

    void DFS_visit(int u) {
        stack<edge_t> s;

        edge_t uEdge;
        uEdge.toNode = u;
        s.push(uEdge);

        while (!s.empty()) {
            edge_t xEdge = s.top();
            int x = xEdge.toNode;
            s.pop();

            if (color[x] == 'w') {
                //Original DFS code
                time = time + 1;
                color[x] = 'g';
                s.push(xEdge);
                vector<edge_t> adjx = adjList.at(x);

                // Now our branching code ::

                // For a white x
                // Consider 2 case:
                // Case 1. p[x] = -1, it can happen in two way, x is the first one ever in this connected component, or no one wanted to take x
                // either way, if p[x] = -1, i can be representative of a new node in new graph
                // Case 2. p[x] != -1, so x won't be the representative/head of a newHome. x just gets added to its parent's newHome.
                int u = unitigs.at(x).ln; //unitig length
                if (p[x] == -1) {
                    oldToNew[x].serial = countNewNode++; // countNewNode starts at 0, then keeps increasing   

                    //make the sequence
                    //NOT CORRECT? I am not sure
                    newSequences.push_back(unitigs.at(x).sequence);

                    oldToNew[x].startPos = 1;
                    if (u <= K) {
                        oldToNew[x].endPos = 1; // do we actually see this? yes
                        //cout<< "u<=k???"<<endl;
                    } else {
                        oldToNew[x].endPos = u - K + 1;
                    }

                } else {
                    oldToNew[x].serial = oldToNew[p[x]].serial;
                    oldToNew[x].startPos = oldToNew[p[x]].endPos + 1;
                    if (u <= K) {
                        oldToNew[x].endPos = oldToNew[x].startPos + 1; // do we actually see this? yes
                        //cout<< "u<=k???"<<endl;
                    } else {
                        oldToNew[x].endPos = u - K + (oldToNew[x].startPos); //check correctness
                    }

                    // x says: Now that I know where my newHome is: I can extend my parent's sequence
                    // Is it more complicated than this?
                    string parentSeq = newSequences.at(oldToNew[x].serial);
                    string childSeq = unitigs.at(x).sequence;

                    // Is it CORRECT? just for testing now
                    newSequences.at(oldToNew[x].serial) = plus_strings(parentSeq, childSeq, K);
                }

                // x->y is the edge, x is the parent we are extending
                for (edge_t yEdge : adjx) { //edge_t yEdge = adjx.at(i);
                    int y = yEdge.toNode;

                    if (DEBUG > 0) {
                        cout << "Edge " << x << "->" << y << endl;
                    }

                    //Normal DFS
                    if (color[y] == 'w') {
                        s.push(yEdge);
                    }

                    if(DEBUG > 4){
                        // DEBUGGING a particular edge
                        if (y == 2 && x == 0) {
                            cout << "Saturated? " << saturated[x] << endl;
                        }
                    }
                    
                    //handle self-loop, self-loop will always be an extra edge
                    if (y == x) {
                        edge_both_t e;
                        e.edge = yEdge;
                        e.fromNode = x;
                        resolveLaterEdges.push_back(e);
                    } else if (saturated[x]) {
                        // Since x is saturated, we only add resolveLater edges
                        // no need to check for consistency
                        if (y != p[x]) {
                            edge_both_t e;
                            e.edge = yEdge;
                            e.fromNode = x;
                            resolveLaterEdges.push_back(e);
                        }
                    } else {
                        // If x has space to take a child, meaning x is not saturated
                        // hunting for potential child

                        if (color[y] == 'w' && p[y] == -1) {
                            // y has white color & got no parent => means it's homeless, so let's see if we can take it as a child of x
                            //But just see if it is eligible to be a child, i.e. is it consistent (sign check)?
                            
                            //2 case, Does x's child have grandparent?
                                // No.
                            if (p[x] == -1) {
                                // case 1: child has no grandparent
                                // so extend path without checking any sign
                                nodeSign[x] = yEdge.left;
                                nodeSign[y] = yEdge.right;
                                p[y] = x;
                                saturated[x] = true; //found a child
                                if (y == 2 && x == 0) {
                                    cout << "HUNT " << saturated[x] << endl;
                                }

                                //TESTED NOT YET
                                if (nodeSign[y] == false) {
                                    unitigs.at(y).sequence = reverseComplement(unitigs.at(y).sequence);
                                }
                                
                                //Yes.
                            } else if (nodeSign[x] == yEdge.left) {
                                // case 2: child (=y) has grandparent, i.e. x's parent exists
                                nodeSign[y] = yEdge.right;
                                p[y] = x;
                                saturated[x] = true; //found a child

                                //TESTED NOT YET
                                if (nodeSign[y] == false) {
                                    unitigs.at(y).sequence = reverseComplement(unitigs.at(y).sequence);
                                }

                            } else {
                                // do we reach this case? 
                                edge_both_t e;
                                e.edge = yEdge;
                                e.fromNode = x;
                                resolveLaterEdges.push_back(e);
                            }
                            
                        } else {
                            if (y != p[x]) {
                                edge_both_t e;
                                e.edge = yEdge;
                                e.fromNode = x;
                                resolveLaterEdges.push_back(e);
                                if (DEBUG > 4) {
                                    // DEBUGGING a particular edge
                                    if (y == 2 && x == 0) {
                                        cout << "Saturated? " << saturated[x] << endl;
                                    }
                                }
                            } else {
                                if (DEBUG > 4) {
                                    // DEBUGGING a particular edge
                                    if (y == 2 && x == 0) {
                                        cout << "Saturated? " << saturated[x] << endl;
                                    }
                                }
                            }


                        }
                    }
                }

            } else if (color[x] == 'g') {
                time = time + 1;
                color[x] = 'b';
            }

        }
    }

    void DFS() {

        for (int i = 0; i < V; i++) {
            color[i] = 'w';
            p[i] = -1;
        }

        for (int i = 0; i < V; i++) {
            if (color[i] == 'w') {
                DFS_visit(i);
            }
        }

        for (int i = 0; i < countNewNode; i++) {
            newAdjList.push_back(vector<newEdge_t>());
        }

        for (edge_both_t e : resolveLaterEdges) {
            // e.start -> e.end : they have different newHome
            // look at the sign of both end points. if from sign and edge from are not equal, then revert the edge label
            // We consider that ACT and CTT has an edge.
            // ACT(+)  (+)->(-)  AAG(+) : this is fine, just add the edge
            // ACT(-)  (+)->(-)  AAG(+) : you need to convert it, because => AGT  (+)->(-)  AAG => this is not correct:
            // to fix this AGT (-)->(-) AAG
            // don't touch the node sign: just fix the edge signs 
            if (nodeSign[e.fromNode] != e.edge.left) {
                e.edge.left = !e.edge.left;
            }
            if (nodeSign[e.edge.toNode] != e.edge.right) {
                e.edge.right = !e.edge.right;
            }

            int x = e.fromNode;
            int u = unitigs.at(x).ln;
            newEdge_t newEdge;
            newEdge.edge = e.edge;
            newEdge.kmerEndIndex = oldToNew[e.edge.toNode].startPos;
            newEdge.kmerStartIndex = oldToNew[x].startPos;

            newAdjList[oldToNew[x].serial].push_back(newEdge);
            cout << "old: " << x << "->" << e.edge.toNode << ", new:" << " (" << oldToNew[x].serial << "->" << newEdge.edge.toNode << ")" << endl;
        }
    }


    ~Graph() {
        delete [] color;
        delete [] p;
        delete [] nodeSign;
        delete [] oldToNew;
        delete [] saturated;
    }
};


void printNewGraph(Graph &G){
    list<int> *newToOld = new list<int>[G.countNewNode];

    // PRINT NEW GRAPH
    for (int i = 0; i < G.V; i++) {
        newToOld[G.oldToNew[i].serial].push_back(i);
        cout << "old " << i << "-> new" << G.oldToNew[i].serial << endl;
    }
    delete [] newToOld;

}

void formattedOutput(Graph &G){
     
    string stitchedUnitigs = "stitchedUnitigs.txt";
     ofstream myfile;
    myfile.open (stitchedUnitigs);
    //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:- 
    for (int newNodeNum = 0; newNodeNum<G.countNewNode; newNodeNum++){
         myfile << '>' << newNodeNum <<" LN:i:"<<newSequences.at(newNodeNum).length()<<" ";
         vector<newEdge_t> edges = newAdjList.at(newNodeNum);
         for(newEdge_t edge: edges){
             myfile<<"L:" << edge.kmerStartIndex << ":" << boolToCharSign(edge.edge.left) << ":" << edge.edge.toNode << ":" << boolToCharSign(edge.edge.left) <<":"<< edge.kmerEndIndex <<" ";
         }
        myfile<<endl;
        myfile<<newSequences.at(newNodeNum);
        myfile<<endl;
    }
    //myfile << '>' << newNodeNum <<">0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:- " ;
    myfile.close();
    
}

int getNodeNumFromFile(string filename) {
    //grep '>' list_reads.unitigs.fa | tail -n 1
    int nodeCount;
    string countFile = "incount.txt";
    ostringstream stringStream;
    stringStream << "grep '>' " << filename << " | tail -n 1" << countFile;
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
    bool doCont = false;

    getline(unitigFile, line);


    do {
        unitig_struct_t unitig_struct;

        sscanf(line.c_str(), "%*c %d %s  %s  %s %[^\n]s", &unitig_struct.serial, lnline, kcline, kmline, edgesline);
        //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-

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

    cout << "Starting reading file: " << UNITIG_FILE << endl;
    if (EXIT_FAILURE == get_data(UNITIG_FILE, data, unitigs, char_count)) {
        return EXIT_FAILURE;
    }

    double TIME_READ_SEC = readTimer() - startTime;

    Graph G;
    cout << G.V << endl;

    G.DFS();

   
    
    
    //count total number of edges
    int E = 0;
    for (int i = 0; i < G.V; i++) {
        E += adjList[i].size();
    }

    int E_new = resolveLaterEdges.size();
    int V = G.V;
    int V_new = G.countNewNode;

    int C = 0;
    int C_new = 0;
    for (unitig_struct_t unitig : unitigs) {
        C += unitig.ln;
    }


    /** --DEBUG: PRINT THE LENGTHS **/
    if (DEBUG > 1) {
        int pp = 0;
        for (string s : newSequences) {
            C_new += s.length();
            cout << pp++ << " -> Length = " << s.length() << endl;
        }
    }




    for (int i = 0; i < newSequences.size(); i++) {

        cout << i << " " << newSequences.at(i) << endl;
        //cout<<i<<" "<<newSequences.at(i).length()<<endl;
    }


    double TIME_TOTAL_SEC = readTimer() - startTime;


    // For collecting stats
    int U_MAX = maximumUnitigLength();
    int EDGE_INT_DTYPE_SIZE;
    if (U_MAX - K + 1 > 0) {
        EDGE_INT_DTYPE_SIZE = log2(U_MAX - K + 1);
    } else {
        EDGE_INT_DTYPE_SIZE = log2(K);
    }
    EDGE_INT_DTYPE_SIZE = ceil(EDGE_INT_DTYPE_SIZE / 8.0);


    int ACGT_DTYPE_SIZE = 1; // 1 byte to store each char
    int NODENUM_DTYPE_SIZE = 8; // 1 byte to store each char
    int SIGN_DTYPE_SIZE = 1; // 1 byte for sign information (+,- can be stored in 1 byte)
    int spaceBefore = C * ACGT_DTYPE_SIZE + E * (NODENUM_DTYPE_SIZE + SIGN_DTYPE_SIZE);
    int save = (C - C_new) * ACGT_DTYPE_SIZE + (E - E_new)*(NODENUM_DTYPE_SIZE + SIGN_DTYPE_SIZE);
    int overhead = (E_new)*(2 * EDGE_INT_DTYPE_SIZE);
    float persaved = ((save - overhead)*1.0 / spaceBefore) * 100.0;

    //    cout << "Time for loading the data: " << TIME_READ_SEC << " sec" << endl;
    //    cout << "Total construction time: " << TIME_TOTAL_SEC<< " sec" << endl;
    //    cout << "V: " << V << endl;
    //    cout << "V_new: " << V_new << endl;
    //    cout << "E: " << E << endl;
    //    cout << "E_new: " << E_new << endl;
    //    cout << "C: " << C << endl;
    //    cout << "C_new: " << C_new << endl;
    //    cout << "Number of Bytes required for storing one edge start or end: " << EDGE_INT_DTYPE_SIZE << endl;
    //    cout << "Space saved: " << save - overhead << " bytes." << endl;
    //    cout << "Space before: " << spaceBefore << " bytes." << endl;
    //    cout << "Percent saved: " << ((save - overhead)*1.0 / spaceBefore) * 100.0 << "%" << endl;
    formattedOutput(G);
    printf("%d \t %d \t %d \t %d \t %d \t %d \t %f \t %f \t %.2f%% \t %d \t %d \t %f \t %f\n", V, V_new, E, E_new, C, C_new, spaceBefore / 1024.0, (save - overhead) / 1024.0, persaved, U_MAX, K, TIME_READ_SEC, TIME_TOTAL_SEC);
    //printGraph(adjList);
    //printAllSequences(unitigs);
    return EXIT_SUCCESS;
}