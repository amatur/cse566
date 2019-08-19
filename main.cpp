// --- VERSION 4.1 ----
//upperbound with divide by 2

//Caution:
//removed all self-loops


#include <cmath>
#include<cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
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
#include <utility>
#include <queue>
#include <deque>
#include <unistd.h>

using namespace std;

int K = 0;
string UNITIG_FILE;

enum DEBUGFLAG_T { NONE = 0,  UKDEBUG = 0, VERIFYINPUT = 1, INDEGREEPRINT = 2, DFSDEBUGG = 3, PARTICULAR = 4, NODENUMBER_DBG = 5, OLDNEWMAP = 9, PRINTER = 10, SINKSOURCE = 12};

enum ALGOMODE_T { BASIC = 0, INDEGREE_DFS = 1, INDEGREE_DFS_1 = 2, OUTDEGREE_DFS = 3, OUTDEGREE_DFS_1 = 4, INDEGREE_DFS_INVERTED = 5, PLUS_INDEGREE_DFS = 6, RANDOM_DFS = 7, NODEASSIGN = 8, SOURCEFIRST = 9, NEWMETHOD = 10, PROFILE_ONLY = 14, EPPRIOR=12, GRAPHPRINT = 13, TIGHTUB = 14};

DEBUGFLAG_T DBGFLAG = NODENUMBER_DBG;
ALGOMODE_T ALGOMODE = PROFILE_ONLY;

string mapmode[] = {"basic", "indegree_dfs", "indegree_dfs_initial_sort_only", "outdegree_dfs", "outdegree_dfs_initial_sort_only", "inverted_indegree_dfs", "plus_indegree_dfs", "random_dfs", "node_assign", "source_first", "new_method", "profile_only", "endpoint_priority", "graph_print", "tight_ub"
};

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

int isolated_node_count = 0;
int sink_count = 0;
int source_count = 0;
int sharedparent_count = 0;
int sharparentCntRefined = 0;
int onecount = 0;


struct node_sorter {
    int node;
    int sortkey;
    //bool operator() (struct node_sorter  i, struct node_sorter  j) { return (i.sortkey<j.sortkey);}
};
bool sort_by_key (struct node_sorter i, struct node_sorter j) { return (i.sortkey<j.sortkey); }
bool sort_by_key_inverted (struct node_sorter i, struct node_sorter j) { return (i.sortkey>j.sortkey); }

int* global_indegree;
int* global_outdegree;
int* global_plusindegree;
int* global_plusoutdegree;
//bool* global_selfloop;

int* global_issinksource;
int* global_priority;

map<pair <int, int>, int> inOutCombo;

vector<vector<edge_t> > adjList;
vector<vector<edge_t> > reverseAdjList;
vector<vector<newEdge_t> > newAdjList;
vector<edge_both_t> resolveLaterEdges;
vector<unitig_struct_t> unitigs;
map<int, string> newSequences;
vector<list<int> > newToOld;


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



class GroupMerger {
public:
    map<int, bool> fwdVisited;
    map<int, bool> bwdVisited;
    map<int, int> bwdWalkId;
    map<int, int> fwdWalkId;
    GroupMerger() {
    }
    void connectGroups(int from, int to){
        fwdVisited[from] = false;
        bwdVisited[to] = false;
        fwdWalkId[from] = to;
        bwdWalkId[to] = from;
    }
    ~GroupMerger() {
    }
};


class DisjointSet {
    unordered_map<int, int> parent;

public:
    DisjointSet() {
    }
    
    
    
    void make_set(int id) {
        this->parent[id] = -1;
    }
    
    void Union(int xId, int yId) {
       int xset = find_set(xId);
       int yset = find_set(yId);
        
        if(xset != yset)
        {
            parent[xset] = yset;
        }
        
        
//        // Attach smaller rank tree under root of high rank tree
//        // (Union by Rank)
//        if (subsets[xroot].rank < subsets[yroot].rank)
//            subsets[xroot].parent = yroot;
//        else if (subsets[xroot].rank > subsets[yroot].rank)
//            subsets[yroot].parent = xroot;
//
//        // If ranks are same, then make one as root and increment
//        // its rank by one
//        else
//        {
//            subsets[yroot].parent = xroot;
//            subsets[xroot].rank++;
//        }
    }
    
    int find_set(int id) {
//        if(parent.count(id)==0){
//            exit(EXIT_FAILURE);
//            return -1;
//        }
        if (parent[id] == -1)
            return id;
        return find_set(parent[id]);
        
        
//        // find root and make root as parent of i (path compression)
//        if (subsets[i].parent != i)
//            subsets[i].parent = find(subsets, subsets[i].parent);
//
//        return subsets[i].parent;
    }
    
    
    
    ~DisjointSet(){
    }
    
};




class Graph {
public:
    int V = adjList.size();
    int countNewNode = 0;
    int time = 0;
    
    char* color;
    int* p;
    bool* nodeSign;
    new_node_info_t* oldToNew;
    bool* saturated;
    struct node_sorter * indegree;
    struct node_sorter * outdegree;
    bool* counted;
    DisjointSet disSet;
    GroupMerger gmerge;
    
    Graph() {
        color = new char[V];
        p = new int[V];
        nodeSign = new bool[V];
        oldToNew = new new_node_info_t[V];
        saturated = new bool[V];
        indegree = new struct node_sorter[V];
        global_indegree = new int[V];
        global_outdegree = new int[V];
        global_plusindegree = new int[V];
        global_plusoutdegree = new int[V];
        //global_selfloop = new bool[V];
        global_issinksource = new int[V];
        global_priority = new int[V];
        counted = new bool[V];
        
        
        
        for (int i = 0; i < V; i++) {
            
            if(ALGOMODE == NEWMETHOD){
                disSet.make_set(i);
            }
            
            
            oldToNew[i].serial = -1;
            saturated[i] = false;
            indegree[i].sortkey = 0;
            indegree[i].node = i;
            global_indegree[i] = 0;
            global_outdegree[i] = 0;
            global_plusindegree[i] = 0;
            global_plusoutdegree[i] = 0;
            //global_selfloop[i] = false;
            global_issinksource[i] = 0;
            global_priority[i] = 0;
            counted[i] = false;
        }
    }
    
    
    inline bool sRight(edge_t plusminusedge){
        return !(plusminusedge.right == true);
    }
    
    inline bool sLeft(edge_t plusminusedge){
        return (plusminusedge.left == true);
    }
    
    void indegreePopulate(){
        //vector<int> selflooper;
        //vector<int> selflooper_pos;
        
//        int xc = 0;
//        for(vector<edge_t> elist: adjList){
//            int pos = 0;
//            for(edge_t e: elist){
//
//                if(e.toNode == xc){
//                    global_selfloop[xc] = true;
//                    cout<<"self-loop: "<<boolToCharSign(e.left)<<" "<<xc<<" "<<boolToCharSign(e.right)<<endl;
//                    selflooper.push_back(xc);
//                    selflooper_pos.push_back(pos);
//                }
//                pos++;
//
//            }
//            xc++;
//        }
        //cout<<"idnegree: "<<global_indegree[12]<<" "<<"outdegree: "<<global_outdegree[12]<<endl;
        
        
        // REMOVE ALL SELF-LOOPS
        
//        for(int i = 0; i<selflooper.size(); i++){
//            cout<<"before size"<<" "<<selflooper[i]<<":"<<adjList[selflooper[i]].size()<<endl;
//            adjList[selflooper[i]].erase(adjList[selflooper[i]].begin()+selflooper_pos[i]);
//            cout<<"after size"<<" "<<selflooper[i]<<":"<<adjList[selflooper[i]].size()<<endl;
//        }
//        cout<<"All self-loops removed."<<endl;
        
        
        int xc = 0;
        for(vector<edge_t> elist: adjList){
            for(edge_t e: elist){
                global_indegree[e.toNode] += 1;
                indegree[e.toNode].sortkey = indegree[e.toNode].sortkey + 1;
                if(e.right == true){
                    global_plusindegree[e.toNode] += 1;
                }
                if(e.left == true){
                    global_plusoutdegree[xc] += 1;
                }
                
            }
            global_outdegree[xc] = elist.size();
            xc++;
        }
        
        for(int i = 0; i<5; i++){
            for(int j = 0; j<5; j++){
                inOutCombo[make_pair(i,j)] = 0;
            }
        }
        for(int i = 0; i<V; i++){
            int minusindegree = (global_indegree[i] - global_plusindegree[i] );
            int minusoutdegree = (global_outdegree[i] - global_plusoutdegree[i] );
            
            pair<int, int> a;
            a = make_pair(global_plusindegree[i], global_plusoutdegree[i]);
            inOutCombo[a] = (inOutCombo.count(a)  ? inOutCombo[a] + 1 : 1  );
            
            
            
            if(DBGFLAG == SINKSOURCE){
                cout<<i<<"is ";
            }
            
//            if(minusoutdegree == 0 && minusindegree != 0){
//                sink_count++;
//                if(DBGFLAG == SINKSOURCE){
//                    cout<<"sink, ";
//                }
//
//            }
//            if(minusindegree == 0 && minusoutdegree != 0){
//                source_count++;
//                if(DBGFLAG == SINKSOURCE){
//                    cout<<"source, ";
//                }
//            }
            
            
            
            
            if(global_plusoutdegree[i] == 0 && global_plusindegree[i] != 0){
                sink_count++;
                global_issinksource[i] = 1;
                global_priority[i] = 5;
                counted[i] = true;
                
                if(DBGFLAG == SINKSOURCE){
                    cout<<"sink, ";
                }

            }
            if(global_plusindegree[i] == 0 && global_plusoutdegree[i] != 0){
                source_count++;
                global_issinksource[i] = 1;
                global_priority[i] = 5;
                counted[i] = true;
                
                if(DBGFLAG == SINKSOURCE){
                    cout<<"source, ";
                }
            }

            
//                if(global_plusoutdegree[i] == 0){
//                    sink_count++;
//                    if(DBGFLAG == SINKSOURCE){
//                        cout<<"sink, ";
//                    }
//
//                }
//                if(global_plusindegree[i] == 0){
//                    source_count++;
//                    if(DBGFLAG == SINKSOURCE){
//                        cout<<"source, ";
//                    }
//                }
            
            
            //global_outdegree[i] += global_indegree[i];
            if(global_indegree[i] == 0){
                global_issinksource[i] = 1;
                isolated_node_count++;
                
                if(DBGFLAG == SINKSOURCE){
                    cout<<"isolated, ";
                }
            }
            if(global_indegree[i] == 1){
                onecount++;
            }
            
            if(DBGFLAG == SINKSOURCE){
                cout<<endl;
            }
            
        }
        
        xc = 0;
        
        
        for(vector<edge_t> elist: adjList){
            stack<int> countedNodes;
            
            
            for(edge_t e_xy: elist){
                int y = e_xy.toNode;
            }
            
            
            /*

            int neighborCount = 0;
            for(edge_t e_xy: elist){
                int y = e_xy.toNode;
                
//                if(!counted[y] && global_plusindegree[y] == 1 && global_plusoutdegree[y] == 1 && e_xy.right == true ){
//                    neighborCnt++;
//                    counted[y] = true;
//                }
                
                if(!counted[y]){
                    vector<edge_t> adjY = adjList[y];
                    bool eligible = true;
                    for(edge_t e_from_y : adjY){
                        if(e_from_y.toNode!=xc){
                            if(sRight(e_xy) == sLeft(e_from_y) ){
                                eligible = false;
                                break;
                            }
                        }
                    }
                    if(eligible){
                        counted[y] = true;
                        global_priority[y] = 4;
                        neighborCount++;
                        countedNodes.push(y);
                    }
                    
                }
            }
            

            if(global_issinksource[xc] == 1){
                if(neighborCount>1){
                    sharedparent_count += neighborCount - 1 ;
                }else{
                    while(!countedNodes.empty()){
                        counted[countedNodes.top()] = false;
                        countedNodes.pop();
                    }
                }
            }else{
                if(neighborCount>2){
                    sharedparent_count += neighborCount - 2 ;
                }else{
                    while(!countedNodes.empty()){
                        counted[countedNodes.top()] = false;
                        countedNodes.pop();
                    }
                }
            }
            */
            
            
            
            xc++;
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
                if(ALGOMODE == RANDOM_DFS){
                    random_shuffle ( adjx.begin(), adjx.end() );
                }
                
                if(ALGOMODE == EPPRIOR){
                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                         {
                             return global_priority[lhs.toNode]   <  global_priority[rhs.toNode]  ;
                         });
                }
                
                if(ALGOMODE == INDEGREE_DFS){
                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                         {
                             return global_indegree[lhs.toNode] < global_indegree[rhs.toNode];
                         });
                }
                
                if(ALGOMODE == PLUS_INDEGREE_DFS){
                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                         {
                             return global_indegree[lhs.toNode]  - global_plusindegree[lhs.toNode] > global_indegree[lhs.toNode]  - global_plusindegree[rhs.toNode];
                         });
                }
                
                if(ALGOMODE == INDEGREE_DFS_INVERTED){
                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                         {
                             return global_indegree[lhs.toNode] > global_indegree[rhs.toNode];
                         });
//                    for(edge_t eeee: adjx){
//                        cout<< global_indegree[eeee.toNode]<< " ";
//                    }
//                    cout<<endl;
                }
                if (ALGOMODE == OUTDEGREE_DFS){
                    if(p[x] == -1){
                        sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                             {
                                 return global_outdegree[lhs.toNode] < global_outdegree[rhs.toNode];
                             });
                    }else if (nodeSign[x] == false){
                        sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                             {
                                 return global_outdegree[lhs.toNode] - global_plusoutdegree[lhs.toNode] < global_outdegree[rhs.toNode] - global_plusoutdegree[rhs.toNode];
                             });
                    }else if (nodeSign[x] == true){
                        sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                             {
                                 return global_plusoutdegree[lhs.toNode] < global_plusoutdegree[rhs.toNode];
                             });
                    }
                }
                
                
                // Now our branching code ::
                
                // For a white x
                // Consider 2 case:
                // Case 1. p[x] = -1, it can happen in two way, x is the first one ever in this connected component, or no one wanted to take x
                // either way, if p[x] = -1, i can be representative of a new node in new graph
                // Case 2. p[x] != -1, so x won't be the representative/head of a newHome. x just gets added to its parent's newHome.
                int u = unitigs.at(x).ln; //unitig length
                
                if (p[x] == -1) {
                    
                    list<int> xxx;
                    xxx.push_back(x);
                    newToOld.push_back(xxx);
                    oldToNew[x].serial = countNewNode++; // countNewNode starts at 0, then keeps increasing
                    
                    
                    
                    
                    
                    //make the sequence
                    //NOT CORRECT? I am not sure
                    if(nodeSign[x]==false){
                        newSequences[oldToNew[x].serial] = reverseComplement(unitigs.at(x).sequence);
                    }else{
                        newSequences[oldToNew[x].serial] = (unitigs.at(x).sequence);
                    }
                    
                    
                    oldToNew[x].startPos = 1;
                    if (u < K) {
                        oldToNew[x].endPos = 1; // do we actually see this? yes
                        if(DBGFLAG == UKDEBUG){
                            cout<< "node: "<< x<<"u< k ***** u = "<<u<<endl;
                        }
                    } else {
                        oldToNew[x].endPos = u - K + 1;
                    }
                    
                } else {
                    
                    newToOld[oldToNew[p[x]].serial].push_back(x);
                    oldToNew[x].serial = oldToNew[p[x]].serial;
                    
                    
                    if(ALGOMODE==NEWMETHOD){
                        disSet.Union(x, p[x]);
                    }
                    
                    
                    oldToNew[x].startPos = oldToNew[p[x]].endPos + 1;
                    if (u < K) {
                        oldToNew[x].endPos = oldToNew[x].startPos + 1; // do we actually see this? yes
                        if(DBGFLAG == UKDEBUG){
                            cout<< "node: "<< x<<"u< k ***** u = "<<u<<endl;
                        }
                    } else {
                        oldToNew[x].endPos = u - K + (oldToNew[x].startPos); //check correctness
                    }
                    
                    // x says: Now that I know where my newHome is: I can extend my parent's sequence
                    // Is it more complicated than this?
                    string parentSeq = newSequences[oldToNew[x].serial];
                    string childSeq = unitigs.at(x).sequence;
                    
                    // Is it CORRECT? just for testing now
                    if(nodeSign[x]==false){
                        childSeq = reverseComplement(childSeq);
                    }
                    newSequences[oldToNew[x].serial] = plus_strings(parentSeq, childSeq, K);
                }
                
                // x->y is the edge, x is the parent we are extending
                for (edge_t yEdge : adjx) { //edge_t yEdge = adjx.at(i);
                    int y = yEdge.toNode;
                    // cout << "Edge " << x << "->" << y << " "<<global_indegree[y] << endl;
                    if (DBGFLAG == DFSDEBUGG) {
                        cout << "Edge " << x << "->" << y << endl;
                    }
                    
                    //Normal DFS
                    if (color[y] == 'w') {
                        s.push(yEdge);
                    }
                    
                    if(DBGFLAG == PARTICULAR){
                        // DEBUGGGING a particular edge
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
                            if (p[x] == -1 && ALGOMODE != NODEASSIGN) {
                                // case 1: child has no grandparent
                                // so extend path without checking any sign
                                
                                nodeSign[x] = yEdge.left;
                                nodeSign[y] = yEdge.right;
                                p[y] = x;
                                saturated[x] = true; //found a child
                                
                                
                                //TESTED NOT YET
                                //                                if (nodeSign[y] == false) {
                                //                                    unitigs.at(y).sequence = reverseComplement(unitigs.at(y).sequence);
                                //                                }
                                
                                //Yes.
                            } else if (nodeSign[x] == yEdge.left) {
                                // case 2: child (=y) has grandparent, i.e. x's parent exists
                                nodeSign[y] = yEdge.right;
                                p[y] = x;
                                saturated[x] = true; //found a child
                                
                                //TESTED NOT YET
                                //                                if (nodeSign[y] == false) {
                                //                                    unitigs.at(y).sequence = reverseComplement(unitigs.at(y).sequence);
                                //                                }
                                
                            } else {
                                // do we reach this case?
                                edge_both_t e;
                                e.edge = yEdge;
                                e.fromNode = x;
                                resolveLaterEdges.push_back(e);
                            }
                            
                        } else {
                            
                            //merger
                            if(ALGOMODE == NEWMETHOD){
                                // y is not white
                                bool consistentEdge = (nodeSign[y] == yEdge.right && (p[x]==-1 || (p[x]!=-1&& nodeSign[x] == yEdge.left)) );
                                if(p[y]==-1 && consistentEdge && oldToNew[x].serial != oldToNew[y].serial){
                                    
                                    //cout<<"x: "<<x<<":" <<disSet.find_set(x)<<" ";
                                    //cout<<"y: "<<y<<":" <<disSet.find_set(y) <<endl;
                                    
                                    //not in same group already, prevent cycle
                                    if(disSet.find_set(x)!=disSet.find_set(y)){
                                        
                                       nodeSign[x] = yEdge.left;
                                       nodeSign[y] = yEdge.right;
                                       p[y] = x;
                                       saturated[x] = true; //found a child
                                       // oldToNew[y].serial
                                        
                                        disSet.Union(x, y);
                                        //cout<<"x: "<<disSet.find_set(x);
                                        //cout<<"y: "<<disSet.find_set(y);
                                        cout<<endl;
                                        
                                        //cout<<oldToNew[x].serial<<"+"<<oldToNew[y].serial<<endl;
                                        gmerge.connectGroups(oldToNew[x].serial,oldToNew[y].serial );
                                        
                                    }
                                    
                                }

                            }
                            
                            if (y != p[x]) {
                                edge_both_t e;
                                e.edge = yEdge;
                                e.fromNode = x;
                                resolveLaterEdges.push_back(e);
                                if (DBGFLAG == PARTICULAR) {
                                    // DEBUGGGING a particular edge
                                    if (y == 2 && x == 0) {
                                        cout << "Saturated? " << saturated[x] << endl;
                                    }
                                }
                            } else {
                                if (DBGFLAG == PARTICULAR) {
                                    // DEBUGGGING a particular edge
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
   
    
    bool canReachSinkSource(int v, bool visited[], bool sign)
    {
        // Mark the current node as visited and
        // print it
        
        
        
        
        visited[v] = true;
        //cout << v << " ";
        bool reachable = false;
        
        if(global_plusoutdegree[v] == 0 && global_plusindegree[v] != 0){
            //cout<<v<<"is sink.";
            return true;//sink
            
        }
        if(global_plusindegree[v] == 0 && global_plusoutdegree[v] != 0){
             //cout<<v<<"is source.";
            return true;//source
        }
        if(global_indegree[v] == 0){
            //cout<<v<<"is isolated.";
            return true;//isolated
        }
        
        
        // Recur for all the vertices adjacent
        // to this vertex
        vector<edge_t>::iterator i;
        for (i = adjList[v].begin(); i != adjList[v].end(); ++i){
            
            if (!visited[(*i).toNode] && sign==(*i).left){
                reachable = canReachSinkSource((*i).toNode, visited, (*i).right);
                if(reachable==true){
                    return true;
                }
            }
            
        }
        

        return reachable;
        
    }
    
    void DFS() {
       
        
        if(ALGOMODE == NODEASSIGN){
            for (int i=0; i<V; i++) {
                nodeSign[i] = true;
//                if(global_plusindegree[i] +  global_plusoutdegree[i] < global_outdegree[i] - global_plusoutdegree[i] + global_indegree[i] - global_plusindegree[i] ){
                if(global_plusindegree[i]< global_indegree[i] - global_plusindegree[i]){
                    nodeSign[i] = true;
                }
            }
        }
//        for(int i = 0; i<V; i++){
//            cout<<global_outdegree[i] << " ";
//        }
//
        
        if(ALGOMODE == SOURCEFIRST){
            for (int i = 0; i < V; i++) {
                indegree[i].node = i;
                indegree[i].sortkey = global_issinksource[i];
            }
            vector<struct node_sorter> myvector (indegree, indegree+V);
            sort (myvector.begin(), myvector.end(), sort_by_key_inverted);
            copy(myvector.begin(), myvector.end(), indegree);
        }
        
        if(ALGOMODE == EPPRIOR){
            for (int i = 0; i < V; i++) {
                indegree[i].node = i;
                indegree[i].sortkey = global_priority[i];
            }
            vector<struct node_sorter> myvector (indegree, indegree+V);
            sort (myvector.begin(), myvector.end(), sort_by_key_inverted);
            copy(myvector.begin(), myvector.end(), indegree);
        }
        
        if(ALGOMODE == INDEGREE_DFS_INVERTED){
            for (int i = 0; i < V; i++) {
                indegree[i].node = i;
                indegree[i].sortkey = global_indegree[i];
            }
            vector<struct node_sorter> myvector (indegree, indegree+V);
            sort (myvector.begin(), myvector.end(), sort_by_key_inverted);
            //random_shuffle ( myvector.begin(), myvector.end() );
            copy(myvector.begin(), myvector.end(), indegree);
            
        }
        
        if (ALGOMODE == INDEGREE_DFS || ALGOMODE == INDEGREE_DFS_1 ){
            for (int i = 0; i < V; i++) {
                indegree[i].node = i;
                indegree[i].sortkey = global_indegree[i];
            }
            vector<struct node_sorter> myvector (indegree, indegree+V);
            sort (myvector.begin(), myvector.end(), sort_by_key);
            copy(myvector.begin(), myvector.end(), indegree);
            
            if(DBGFLAG == INDEGREEPRINT){
                cout<<"print in degrees"<<endl;
                for(int i = 0; i<V; i++){
                    cout<<indegree[i].node<<"->"<<indegree[i].sortkey<<endl;
                }
            }
        }
        
        
        if(ALGOMODE == RANDOM_DFS){
            for (int i = 0; i < V; i++) {
                indegree[i].node = i;
                indegree[i].sortkey = global_indegree[i];
            }
            vector<struct node_sorter> myvector (indegree, indegree+V);
            sort (myvector.begin(), myvector.end(), sort_by_key);
            random_shuffle ( myvector.begin(), myvector.end() );
            copy(myvector.begin(), myvector.end(), indegree);
            
        }
        
        
        if (ALGOMODE == OUTDEGREE_DFS || ALGOMODE == OUTDEGREE_DFS_1){
            for (int i = 0; i < V; i++) {
                indegree[i].node = i;
                indegree[i].sortkey = global_outdegree[i];
            }
            vector<struct node_sorter> myvector (indegree, indegree+V);
            sort (myvector.begin(), myvector.end(), sort_by_key);
            copy(myvector.begin(), myvector.end(), indegree);
        }
        
        double time_a = readTimer();
        for (int i = 0; i < V; i++) {
            color[i] = 'w';
            p[i] = -1;
        }
        cout<<"Basic V loop time: "<<readTimer() - time_a<<" sec"<<endl;
        
        
        
        
        time_a = readTimer();
     
        for (int j = 0; j < V; j++) {

            
            
            int i;
            if(ALGOMODE == OUTDEGREE_DFS || ALGOMODE == OUTDEGREE_DFS_1 || ALGOMODE == INDEGREE_DFS || ALGOMODE == INDEGREE_DFS_1 || ALGOMODE == SOURCEFIRST){
                i = indegree[j].node;
            }else{
                i = j;
            }
            
            if (color[i] == 'w') {
                if(DBGFLAG == DFSDEBUGG ){
                    cout<<"visit start of node: "<<i<<endl;
                }
                DFS_visit(i);
            }
        }
        cout<<"DFS time: "<<readTimer() - time_a<<" sec"<<endl;
        
        
        //reuse the colors
        //~ for (int i = 0; i < V; i++) {
            //~ color[i] = 'w';
        //~ }
        
        
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
            
            int x = e.fromNode;
            int u = unitigs.at(x).ln;
            newEdge_t newEdge;
            newEdge.kmerEndIndex = oldToNew[e.edge.toNode].startPos;
            newEdge.kmerStartIndex = oldToNew[x].endPos;
            
            if(e.fromNode!=e.edge.toNode) {
                if (nodeSign[e.fromNode] != e.edge.left) {
                    e.edge.left = !e.edge.left;
                    newEdge.kmerStartIndex = oldToNew[x].startPos;
                }
                if (nodeSign[e.edge.toNode] != e.edge.right) {
                    e.edge.right = !e.edge.right;
                    newEdge.kmerEndIndex = oldToNew[e.edge.toNode].endPos;
                }
            }
            
            
            newEdge.edge = e.edge;
            newEdge.edge.toNode = oldToNew[newEdge.edge.toNode].serial;
            
            newAdjList[oldToNew[x].serial].push_back(newEdge);
            
            if(DBGFLAG == OLDNEWMAP){
                cout << "old: " << x << "->" << e.edge.toNode << ", new:" << " (" << oldToNew[x].serial << "->" << newEdge.edge.toNode << ")" << endl;
                
            }
        }
    }
    
    
    ~Graph() {
        delete [] color;
        delete [] p;
        delete [] nodeSign;
        delete [] oldToNew;
        delete [] saturated;
        delete [] indegree;
        delete [] counted;
//        delete [] global_indegree;
//        delete [] global_outdegree;
//        delete [] global_plusindegree;
//        delete [] global_plusoutdegree;
        //delete [] global_selfloop;
    }
};


void printNewGraph(Graph &G){
    list<int> *newToOld = new list<int>[G.countNewNode];
    
    // PRINT NEW GRAPH
    for (int i = 0; i < G.V; i++) {
        newToOld[G.oldToNew[i].serial].push_back(i);
        //cout << "old " << i << "-> new" << G.oldToNew[i].serial << endl;
    }
    
    for (int i = 0; i < G.countNewNode; i++) {
        list<int> adj = newToOld[i];
        
        cout<<"new " << i<<": old (index) ";
        for(int val : adj){
            cout<<val<<" ";
        }
        cout<<endl;
        
    }
    
    delete [] newToOld;
    
}


void tableDegreeDist(Graph &G){
    string plainOutput = "plainOutput.fa";
    ofstream plainfile;
    plainfile.open(plainOutput);
    
    string stitchedUnitigs = "stitchedUnitigs.txt";
    ofstream myfile;
    myfile.open (stitchedUnitigs);
    //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-
    for (int newNodeNum = 0; newNodeNum<G.countNewNode; newNodeNum++){
        myfile << '>' << newNodeNum <<" LN:i:"<<newSequences[newNodeNum].length()<<" ";
        //plainfile << '>' << newNodeNum;
        
        vector<newEdge_t> edges = newAdjList.at(newNodeNum);
        for(newEdge_t edge: edges){
            myfile<<"L:" << edge.kmerStartIndex << ":" << boolToCharSign(edge.edge.left) << ":" << edge.edge.toNode << ":" << boolToCharSign(edge.edge.left) <<":"<< edge.kmerEndIndex <<" ";
        }
        //plainfile<<endl;
        myfile<<endl;
        
        plainfile<<newSequences[newNodeNum];
        myfile<<newSequences[newNodeNum];
        
        plainfile<<endl;
        myfile<<endl;
    }
    //myfile << '>' << newNodeNum <<">0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:- " ;
    myfile.close();
    plainfile.close();
    
}



void formattedOutput(Graph &G){
    string plainOutput = "plainOutput.fa";
    ofstream plainfile;
    plainfile.open(plainOutput);
    
    string stitchedUnitigs = "stitchedUnitigs.txt";
    ofstream myfile;
    myfile.open (stitchedUnitigs);
    //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-
    for (int newNodeNum = 0; newNodeNum<G.countNewNode; newNodeNum++){
        myfile << '>' << newNodeNum <<" LN:i:"<<newSequences[newNodeNum].length()<<" ";
        //plainfile << '>' << newNodeNum;
        
        vector<newEdge_t> edges = newAdjList.at(newNodeNum);
        for(newEdge_t edge: edges){
            myfile<<"L:" << edge.kmerStartIndex << ":" << boolToCharSign(edge.edge.left) << ":" << edge.edge.toNode << ":" << boolToCharSign(edge.edge.left) <<":"<< edge.kmerEndIndex <<" ";
        }
        //plainfile<<endl;
        myfile<<endl;
        
        plainfile<<newSequences[newNodeNum];
        myfile<<newSequences[newNodeNum];
        
        plainfile<<endl;
        myfile<<endl;
    }
    //myfile << '>' << newNodeNum <<">0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:- " ;
    myfile.close();
    plainfile.close();
    
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
    unitigFile.open(unitigFileName);
    
    
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
        edgesline[0] = '\0';
        sscanf(line.c_str(), "%*c %d %s  %s  %s %[^\n]s", &unitig_struct.serial, lnline, kcline, kmline, edgesline);
//        if(unitig_struct.serial == 1241914){
//            cout<<line<<endl;
//            cout<<edgesline<<endl;
//        }
        //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-
        
        sscanf(lnline, "%*5c %d", &unitig_struct.ln);
        sscanf(kcline, "%*5c %d", &unitig_struct.kc);
        sscanf(kmline, "%*5c %f", &unitig_struct.km);
        
        char c1, c2;
        stringstream ss(edgesline);
        
        vector<edge_t> edges;
        while (getline(ss, line, ' ')) {
            if (delSpaces(line).length() != 0) {
                if(DBGFLAG==VERIFYINPUT){
                    cout<<line<<endl;
                }
//                if(unitig_struct.serial == 1241914){
//                    cout<<line<<endl;
//                }
                sscanf(line.c_str(), "%*2c %c %*c %d  %*c  %c", &c1, &nodeNum, &c2); //L:-:0:-
                edge_t newEdge;
                
                bool DELSELFLOOP=true;
                if(DELSELFLOOP){
                    if((unitig_struct.serial)!= nodeNum){
                        newEdge.left = charToBool(c1);
                        newEdge.right = charToBool(c2);
                        newEdge.toNode = nodeNum;
                        edges.push_back(newEdge);
                    }
                }else{
                    newEdge.left = charToBool(c1);
                    newEdge.right = charToBool(c2);
                    newEdge.toNode = nodeNum;
                    edges.push_back(newEdge);
                }
                
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

int getFileSizeBits(string fn = "plainOutput.fa.gz"){
    system("gzip plainOutput.fa");
    int count;
    string line;
    string countFile = "incount.txt";
    ostringstream stringStream;
    //du -h plainOutput.fa | cut -f1
    stringStream << "du -k " << fn << " | cut -f1 > " << countFile;
    string copyOfStr = stringStream.str();
    system(copyOfStr.c_str());
    ifstream cf;
    cf.open(countFile);
    getline(cf, line);
    line = delSpaces(line);
    sscanf(line.c_str(), "%d", &count);
    cf.close();
    return count*1024*8;
}

set<int> extractIntegerWords(string str)
{
    set<int> retset;
    stringstream ss;
    
    /* Storing the whole string into string stream */
    ss << str;
    
    /* Running loop till the end of the stream */
    string temp;
    int found;
    while (!ss.eof()) {
        
        /* extracting word by word from stream */
        ss >> temp;
        
        /* Checking the given word is integer or not */
        if (stringstream(temp) >> found)
            retset.insert(found);
        
        /* To save from space at the end of string */
        temp = "";
    }
    return retset;
}

void makeGraphDot(string ipstr){
    FILE * fp;
    
    fp = fopen ("/Users/Sherlock/Downloads/graphviz-2.40.1/graph.gv", "w+");
   
    fprintf(fp, "digraph d {\n");
    //string ipstr = "20 19 18";
    set<int> verticesMarked;
    set<int> vertices = extractIntegerWords(ipstr) ;
    set<int> vMarked(vertices.begin(), vertices.end());
    //set<pair<int, int> > edges;
    
    for(int x: vertices){
        if(x>=adjList.size()){
            cout<<"wrong, do again"<<endl;
            return;
        }
        vector<edge_t> adjX = adjList[x];
        for(edge_t ex: adjX){
            vertices.insert(ex.toNode);
            fprintf(fp, "%d -> %d[taillabel=\"%d\", headlabel=\"%d\", arrowhead=\"none\"]\n", x, ex.toNode, ex.left, !ex.right);
            
            //            pair<int, int> p;
            //            if(x < ex.toNode){
            //                p.first = x;
            //                p.second = ex.toNode;
            //            }else{
            //                p.second = x;
            //                p.first = ex.toNode;
            //            }
            //            edges.insert(p);
        }
    }
    for(int x: vertices){
        if(vMarked.count(x)>0){
            fprintf(fp, "%d [label=\"%d\", color=\"red\"]\n", x, x);
        }else{
            fprintf(fp, "%d [label=\"%d\"]\n", x, x);
        }
        
    }
    
    //for all int in list
    //make a list of neighbors add them
    
    
    fprintf(fp, "}\n");
    
    fclose(fp);
}

int main(int argc, char** argv) {
    const char* nvalue = "" ;
    
    int c ;
    while( ( c = getopt (argc, argv, "i:k:m:d:") ) != -1 )
    {
        switch(c)
        {
            case 'i':
                if(optarg) nvalue = optarg;
                break;
            case 'm':
                if(optarg) {
                    ALGOMODE = static_cast<ALGOMODE_T>(std::atoi(optarg));
                }
                break;
            case 'd':
                if(optarg) {
                    DBGFLAG = static_cast<DEBUGFLAG_T>(std::atoi(optarg));
                }
                break;
            case 'k':
                if(optarg) {
                    K = std::atoi(optarg) ;
                    if(K<=0){
                        fprintf(stderr, "Error: Specify a positive k value.\n");
                        exit(EXIT_FAILURE);
                    }
                }else{
                    fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name>\n",
                            argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
             default: /* '?' */
                fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name>\n",
                        argv[0]);
             exit(EXIT_FAILURE);
             
        }
    }

    if(K==0 || strcmp(nvalue, "")==0){
        fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name>\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }
    
    
    
    UNITIG_FILE = string(nvalue);
    
    
    ifstream infile(UNITIG_FILE);
    if(!infile.good()){
        fprintf(stderr, "Error: File named \"%s\" cannot be opened.\n", UNITIG_FILE.c_str());
        exit(EXIT_FAILURE);
    }
    
    
    uint64_t char_count;
    uchar *data = NULL;
    
    double startTime = readTimer();
    cout << "## START reading file: " << UNITIG_FILE << ": K = "<<K<<endl;
    if (EXIT_FAILURE == get_data(UNITIG_FILE, data, unitigs, char_count)) {
        return EXIT_FAILURE;
    }
    double TIME_READ_SEC = readTimer() - startTime;
    cout<<"TIME to read file "<<TIME_READ_SEC<<" sec."<<endl;
    

    Graph G;
    
    //count total number of edges
    int E = 0;
    for (int i = 0; i < G.V; i++) {
        E += adjList[i].size();
    }
    int V = G.V;
    int numKmers = 0;
    int C = 0;
    
    for (unitig_struct_t unitig : unitigs) {
        C += unitig.ln;
        numKmers +=  unitig.ln - K + 1;
    }
    
    
    if(DBGFLAG == NODENUMBER_DBG){
        cout<<"Total Nodes: "<<V<<" Edges: "<<E<<" K-mers: "<<numKmers<<endl;
    }
    
    
    cout<<"## START gathering info about upper bound. "<<endl;
    double time_a = readTimer();
    G.indegreePopulate();
    //        bool visitedForReachable[V];
    //        for (int i = 0; i< V; i++) {
    //            visitedForReachable[i]  = false;
    //        }
    //        canReachSinkSource(int v, bool *visited)
    //bool visitedForReachable[V];
    cout<<"TIME for information gather: "<<readTimer() - time_a<<" sec."<<endl;
    
    delete [] global_indegree;
    delete [] global_outdegree;
    delete [] global_plusindegree;
    delete [] global_plusoutdegree;
    
    if(ALGOMODE == GRAPHPRINT){
        char sss[1000];
        cout<<"say something: "<<endl;
        while(true){
            //string ipstr = "20 19 18";
            gets(sss);
            string ipstr(sss);
            if(ipstr=="stop"){
                break;
            }
            makeGraphDot(ipstr);
            cout<<"done print, say again:"<<endl;
        }
    }
    
    
    if(ALGOMODE == PROFILE_ONLY){
        return 0;
    }
    
    
    
    int walkstarting_node_count = ceil((sharedparent_count + sink_count + source_count)/2.0) + isolated_node_count;
    float upperbound = (1-((C-(K-1)*(G.V - walkstarting_node_count*1.0))/C))*100.0;
    float upperboundLoose = (1-((C-(K-1)*(G.V - 1.0))/C))*100.0;
   
    
    printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.6f%%\t%.6f%%\n", K, numKmers, C, V, E, isolated_node_count, onecount, sink_count, source_count, sharedparent_count, upperbound, upperboundLoose);
    
    // Iterating the map and printing ordered values
    for (auto i = inOutCombo.begin(); i != inOutCombo.end(); i++) {
        cout << "(" << i->first.first<< ", "<< i->first.second << ")" << " := " << i->second << '\n';
    }
    
    if(ALGOMODE == PROFILE_ONLY){
        return 0;
    }
    
    
    //STARTING DFS
    cout<<"## START DFS: "<<endl;
    G.DFS();
    
    
    if(DBGFLAG == PRINTER){
        printBCALMGraph(adjList);
        printNewGraph(G);
        for(int i = 0; i< G.countNewNode; i++){
            cout<<"new ->" <<i<<" ";
            for(int x: newToOld[i]){
                cout<<x<<" ";
            }
            cout<<endl;
        }
    }
    
    
    cout<<"## START stitching strings: "<<endl;
    time_a = readTimer();
    //fix sequences
    for(int i = 0; i< G.countNewNode; i++){
        string s = "";
        for(int x: newToOld[i]){
            if(G.nodeSign[x] == false){
                s = plus_strings(s, reverseComplement(unitigs.at(x).sequence), K);
            }else{
                s = plus_strings(s, (unitigs.at(x).sequence), K);
            }
        }
        newSequences[i] = s;
        //cout<<endl;
    }
    cout<<"TIME to stitch: "<<readTimer() - time_a<<" sec."<<endl;

    
    //Stats after all done
    time_a = readTimer();
    int E_new = resolveLaterEdges.size();
    int V_new = G.countNewNode;
    int C_new = 0;
 
    map<int, string>::iterator it;
    int maxlen = 0;
    for (it = newSequences.begin(); it != newSequences.end(); it++)
    {
        
        C_new += (it->second).length();
        if((it->second).length() >maxlen){
            maxlen =(it->second).length();
        }
        
    }
    cout<<"TIME to edge resolve: "<<readTimer() - time_a<<" sec."<<endl;

    cout<<"GROUP PRINT"<<endl;
    bool* merged = new bool[G.countNewNode];
    for (int i = 0; i<G.countNewNode; i++) {
        merged[i] = false;
    }
    
    
    /***MERGE START***/
    int C_better = 0;
    int Vcounttttt = 0;
    ofstream betterfile;
    betterfile.open("betterUnitigs.fa");

    
    for ( const auto& p: G.gmerge.fwdWalkId)
    {
        if(G.gmerge.fwdVisited[p.first] == false){
            
            int fromnode =p.first;
            int tonode = p.second;
            deque<int> lst;
            
            lst.push_back(fromnode);
            lst.push_back(tonode);
            
            G.gmerge.fwdVisited[fromnode] = true;
            G.gmerge.bwdVisited[tonode] = true;
            
            
            if(G.gmerge.fwdVisited.count(tonode)>0){
                while(G.gmerge.fwdVisited[tonode] == false){
                    G.gmerge.fwdVisited[tonode] = true;
                    tonode = G.gmerge.fwdWalkId[tonode];
                    G.gmerge.bwdVisited[tonode] = true;
                    
                    lst.push_back(tonode);
                    if(G.gmerge.fwdVisited.count(tonode)==0)
                        break;
                }
            }
            if(G.gmerge.bwdVisited.count(fromnode)>0){
                while(G.gmerge.bwdVisited[fromnode] == false){
                    G.gmerge.bwdVisited[fromnode] = true;
                    fromnode = G.gmerge.bwdWalkId[fromnode];
                    G.gmerge.fwdVisited[fromnode] = true;
                    
                    lst.push_front(fromnode);
                    if(G.gmerge.bwdVisited.count(fromnode)==0)
                        break;
                }
            }
            
            string mergeString = "";
            for(auto i: lst){
                merged[i] = true;
                mergeString = plus_strings(mergeString, newSequences[i], K);
                //cout<<i<<" ";
            }
            
            
            //cout<<endl;
            Vcounttttt ++;
            C_better+=mergeString.length();
            betterfile << '>' << fromnode <<" LN:i:"<<mergeString.length()<<" ";
            betterfile<<endl;
            betterfile<<mergeString;
            betterfile<<endl;
            
        }
    }
    
    
    for (int newNodeNum = 0; newNodeNum<G.countNewNode; newNodeNum++){
        if(merged[newNodeNum] == false){
            betterfile << '>' << newNodeNum <<" LN:i:"<<newSequences[newNodeNum].length()<<" ";
            betterfile<<endl;
            betterfile<<newSequences[newNodeNum];
            betterfile<<endl;
            C_better+=newSequences[newNodeNum].length();
            Vcounttttt++;
        }
        
       
        
    }
    betterfile.close();
    delete[] merged;
    
//    for(map<int, int> a: mergelist){
//        cout<<a.first<<" -> "<<a.second<<endl;
//        string a = plus_strings(u, childSeq, K);
//    }
    
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
    float persaved = ((save - overhead)*1.0 / spaceBefore) * 100.0; //including edge info
    float saved_c = (1-(C_better*1.0/C))*100.0;
    
    printf("%d \t %d \t %d \t %d \t %d \t %d \t %f \t %f \t %.2f%% \t %d \t %d \t %d \t %f \t %f \t %d \t %d \t %d \t %d \t %.6f%% \t %.6f%% \t %s \t %d \t %d \n", V, V_new, E, E_new, C, C_new, spaceBefore / 1024.0, (save - overhead) / 1024.0, persaved, U_MAX, maxlen, K, TIME_READ_SEC, TIME_TOTAL_SEC, isolated_node_count, onecount, sink_count, source_count, upperbound, saved_c, mapmode[ALGOMODE].c_str(), numKmers, sharedparent_count);
    
    time_a = readTimer();
    formattedOutput(G);
    cout<<"TIME to output: "<<readTimer() - time_a<<" sec."<<endl;
    
    //printf(")
    printf("%d \t %d \t %d \t %d \t %d \t %d \t %f \t %f \t %.2f%% \t %d \t %d \t %d \t %f \t %f \t %d \t %d \t %d \t %d \t %.6f%% \t %.6f%% \t %s \t %d \t %d \t %.6f\n", V, Vcounttttt, E, E_new, C, C_better, spaceBefore / 1024.0, (save - overhead) / 1024.0, persaved, U_MAX, maxlen, K, TIME_READ_SEC, TIME_TOTAL_SEC, isolated_node_count, onecount, sink_count, source_count, upperbound, saved_c, mapmode[ALGOMODE].c_str(), numKmers, sharedparent_count, getFileSizeBits()*1.0/numKmers);

    
    return EXIT_SUCCESS;
}
