

//
//
//int glob = 0;
//vector<string> DFS(int k, map<string, vector<string> > & adjList, string v, vector<string> circuit){
//    glob++;
//
//    vector<string> circuitLocal(circuit.begin(),circuit.end());
//
//    if(adjList[v].size() != 0){
//        int n = adjList[v].size() ;
//        int i = n-1;
//        while(true){
//            if(i==-1) break;
//            string neighborString = adjList[v].at(i);
//
//            //remove edge
//            circuitLocal.push_back(adjList[v].at(i));
//
//            adjList[v].erase(adjList[v].begin() + i);
//            DFS(k, adjList, neighborString, circuitLocal);
//
//            adjList[v].push_back(neighborString);
//            circuitLocal.pop_back();
//
//            i--;
//        }
//    }else{
//        if(circuitLocal.size() == kmerList.size()+1){
//            for (int i = 0; i<circuitLocal.size(); i++) {
//                //cout<<circuitLocal.at(i)<< " ";
//            }
//            string out = "";
//            for (int i = 0; i<kmerList.size(); i++) {
//                out+=circuitLocal.at(i)[0];
//            }
//            results.insert(out);
//
//        }
//    }
//    return circuit;
//}





/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 


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
 * */