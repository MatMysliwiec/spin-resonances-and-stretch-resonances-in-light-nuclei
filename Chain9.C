#include <iostream>
#include <fstream> 

TChain * Chain(int run1, int run2){

static const char * filemask = "/home/finger/Pulpit/rewrite_tree/run%04d.root";

    TChain * chain = new TChain("tree","tree");
	
std::ofstream outfile ("runs_from_chain.txt");
outfile << run1 << " " << run2 << std::endl;
outfile.close();

    for (int i=run1; i<=run2; i++)
        chain->Add(Form(filemask, i));
    return chain;
}
