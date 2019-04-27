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