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
    edge_t edge;
    int fromNode;
} edge_both_t;

typedef struct {
    //1 means +, 0 means -
    bool left;
    bool right;
    int toNode;
    int kmerStartIndex;
    int kmerEndIndex;
} newEdge_t;


// ------- PARAMETERS -------- //
int k = 11;
string  unitigFileName  = "data/list_reads.unitigs.fa";


vector<vector<edge_t> > adjList;
vector<vector<edge_t> > newAdjList;
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


// @@ --- ALL PRINTING CODE --- //
void printGraph(vector<vector<edge_t> > adjList){
    for (int i = 0; i < adjList.size(); i++) {
        cout<<i<<"# ";
        for(edge_t edge : adjList.at(i)){
            cout<<boolToCharSign(edge.left)<<":"<<edge.toNode<<":"<<boolToCharSign(edge.right)<<", ";
        }
        cout<<endl;
    }
}
void printAllSequences(vector<unitig_struct_t> unitigs){
    for (unitig_struct_t unitig : unitigs) {
        cout<<unitig.serial<<": "<<unitig.ln<<" "<<unitig.sequence.length()<<endl; //sequence only^
        // full print
        //cout<<unitig.serial<<": "<<unitig.ln<<" "<<unitig.sequence.length()<<":"<<unitig.sequence<<endl;
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
    bool* saturated;
    int countNewNode;
    
    Graph(){
        color = new char[V];
        p = new int[V];
        nodeSign = new bool[V];
        oldToNew = new new_node_info_t[V];
        saturated = new bool[V];
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
                //Original DFS code
                time = time + 1;
                color[x] = 'g';
                s.push(xEdge);
                vector<edge_t> adjx = adjList.at(x);
                
                
                // For a white x
                // p[x] = -1, it can happen in two way, I am the first one ever in this con.component, or no one wanted to take me
                // either way, if p[x] = -1, i can be representative of a new node in new graph
                // p[x] != -1, so I won't be the representative/head of a newHome. I just get added to my parent's newHome.
                
                
                
                int u = unitigs.at(x).ln; //unitig length
                if(p[x] == -1){
                    
                    oldToNew[x].serial = countNewNode++; //start at 0   
                    //make the sequence
                    
                    //NOT CORRECT
                    newSequences.push_back(unitigs.at(x).sequence);
                    
                    oldToNew[x].startPos = 1;
                    if(u <= k){
                        oldToNew[x].endPos = 1; // do we actually see this?
                        cout<< "u<=k???"<<endl;
                    }else{
                        oldToNew[x].endPos = u - k + 1;
                    }
                    
                }else{
                    
                    oldToNew[x].serial = oldToNew[p[x]].serial;
                    oldToNew[x].startPos = oldToNew[p[x]].endPos + 1;
                     if(u <= k){
                        oldToNew[x].endPos = oldToNew[x].startPos + 1; // do we actually see this?
                        cout<< "u<=k???"<<endl;
                    }else{
                        oldToNew[x].endPos =  u - k + (oldToNew[x].startPos); //check correctness
                    }
                    
                    // I know what's my new home now
                    // more complicated than this
                    string parentSeq = newSequences.at(oldToNew[x].serial);
                    string childSeq = unitigs.at(x).sequence;
//                    if(xEdge.left = false){
//                        parentSeq = reverseComplement(parentSeq);
//                    }
                    // NOT CORRECT, just for testing now
                    newSequences.at(oldToNew[x].serial) = plus_strings(parentSeq, childSeq, k);
                }
                
                
                // x->y is the edge, x is the parent we are extending
                for(edge_t yEdge: adjx ){
                    int y = yEdge.toNode;
                    
                    //Normal DFS
                    if(color[y] == 'w'){
                        s.push(yEdge);
                    }
                    
                    //Code for making new graph
                    if(saturated[x]){
                        // ami saturated, ekhon just resolve later edge ber kora
                        // no need to check for consistency
                        
                        edge_both_t e;
                        e.edge = yEdge;
                        e.fromNode = x;
                        resolveLaterEdges.push_back(e);
                        
                    
                    }else{
                        // amar nijer jayga thakle, mane ami jodi saturated na hoi 
                        // hunting for potential child

                        if(color[y] == 'w'){
                        // white mane pobitro homeless obosshoi, eke nite CHABO amar child hishebe
                        // nite chacchi, but dekhte hobe she ki neyar moto kina
                        // is it consistent?
                        //2 case, my child will has grandparent or not
                            if(p[x] == -1){
                                // case 1: child has no grandparent
                                nodeSign[x] = yEdge.left;
                                nodeSign[y] = yEdge.right;
                                p[y] = x;
                                saturated[x] = true;    //found a child
                                
                                
                            }else if(nodeSign[x] == yEdge.left){
                                 // case 2: child has grandparent, my parent exists
                                nodeSign[y] = yEdge.right;
                                p[y] = x;
                                saturated[x] = true;  //found a child
                            }
                            // ei 2 case er konotatei porenai, tahole oke home dite parchi na
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

        for (edge_both_t e : resolveLaterEdges) {
             
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

    cout << "Starting reading file: " << unitigFileName << endl;
    if (EXIT_FAILURE == get_data(unitigFileName, data, unitigs, char_count)) {
        return EXIT_FAILURE;
    }
    
    
    
    Graph G;
    cout<<G.V<<endl;
    
    G.DFS();
    
    list<int> *newToOld = new list<int>[G.countNewNode]; 
    
    //int newToOld = new int[G.V];
    for (int i = 0; i < G.V; i++) {
        newToOld[G.oldToNew[i].serial].push_back(i); 
    }
    
    cout<<"Number of OLD nodes: "<<G.V<<endl;
    cout<<"Number of new nodes: "<<G.countNewNode<<endl;
    
    delete [] newToOld;

//    for (int i = 0; i < newSequences.size(); i++) {
//        //cout<< i<< " " <<newSequences.at(i)<<endl;
//        cout<<i<<" "<<newSequences.at(i).length()<<endl;
//    }
    cout<<"Num edges: " <<resolveLaterEdges.size();

    
    //printGraph(adjList);
    //printAllSequences(unitigs);
    return EXIT_SUCCESS;
}