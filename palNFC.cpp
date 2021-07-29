#include "palabos3D.h"
#include "palabos3D.hh"

#include <cstdlib>
#include <cmath>


using namespace plb;

#define DESCRIPTOR descriptors::D3Q19Descriptor

typedef double T;

std::string outdir(" ");

struct SimulationParameters
{
    std::vector<T> xDomain;
    std::vector<T> yDomain;
    std::vector<T> zDomain;
    
    std::vector<std::string> obstFileName;
    
    int flowDirection;
    T Re;
    T characteristicLength;
    plint resolution;
    T dt;
    
    plint maxIter;
    plint statIter;
    plint outIter;
    plint cpIter;
    
    T rho;
    T nu;
    T cSmago; 


};

void readInputParameters(std::string xmlInputFileName, SimulationParameters& param)
{

}

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir(outdir);
    
    if (argc != 2 && argc != 3) 
    {
        pcout << "Usage: " << argv[0] << " xml-input-file-name [xml-continue-file-name]" << std::endl;
        exit(1);
    }
    
    std::string xmlInputFileName;
    xmlInputFileName = std::string(argv[1]);
    
    std::string xmlRestartFileName;
    bool continueSimulation = false;
    if (argc == 3) 
    {
        xmlRestartFileName = std::string(argv[2]);
        continueSimulation = true;
    }
    
    SimulationParameters param;
    readInputParameters(xmlInputFileName, param);


    return 0;
}
