#include <iostream>
#include <fstream> 

TChain * Chain_full(){

static const char * filemask = "/data2/mbsdaq/e009/root/run%04d.root";

int first = 4383;
int last  = 7838;

std::ofstream outfile ("runs_from_chain.txt");
outfile << first << " " << last << std::endl;
outfile.close();

    TChain * chain = new TChain("tree","tree");

    /*for (int i=4383; i<=4594; i++)
        chain->Add(Form(filemask, i));

    for (int i=4596; i<=4691; i++)
        chain->Add(Form(filemask, i));

    for (int i=4718; i<=4724; i++)
        chain->Add(Form(filemask, i));*/

    for (int i=6984; i<=7245; i++)
        chain->Add(Form(filemask, i));

    for (int i=7368; i<=7395; i++)
        chain->Add(Form(filemask, i));

    for (int i=7568; i<=7647; i++)
        chain->Add(Form(filemask, i));

    for (int i=7649; i<=7706; i++)
        chain->Add(Form(filemask, i));

    for (int i=7709; i<=7838; i++)
        chain->Add(Form(filemask, i));

    return chain;
}
