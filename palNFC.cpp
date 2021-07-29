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
    XMLreader document(xmlInputFileName);
    document["geometry"]["simulationDomain"]["x"].read(param.xDomain);
    PLB_ASSERT(param.xDomain.size() == 2 && param.xDomain[1] > param.xDomain[0]); 
    document["geometry"]["simulationDomain"]["y"].read(param.yDomain);
    PLB_ASSERT(param.yDomain.size() == 2 && param.yDomain[1] > param.yDomain[0]);
    document["geometry"]["simulationDomain"]["z"].read(param.zDomain);
    PLB_ASSERT(param.zDomain.size() == 2 && param.zDomain[1] > param.zDomain[0]);
    
    document["geometry"]["obstacleFileName"].read(param.obstFileName);
    document["geometry"]["flowDirection"].read(param.flowDirection);
    
    document["numerics"]["Re"].read(param.Re);
    PLB_ASSERT(param.Re > 0);
    
    document["numerics"]["characteristicLength"].read(param.characteristicLength);
    PLB_ASSERT(param.characteristicLength > 0);
    
    document["numerics"]["resolution"].read(param.resolution);
    PLB_ASSERT(param.resolution > 0);
    
    document["numerics"]["dt"].read(param.dt);
    PLB_ASSERT(param.dt > 0);
    
    document["numerics"]["maxIter"].read(param.maxIter);
    PLB_ASSERT(param.maxIter > 0);
    
    document["numerics"]["cSmago"].read(param.cSmago);
    PLB_ASSERT(param.cSmago >= 0);
    
    document["fluid"]["nu"].read(param.nu);
    PLB_ASSERT(param.nu > 0);
    
    document["fluid"]["rho"].read(param.rho);
    PLB_ASSERT(param.rho > 0);
    
    document["output"]["statIter"].read(param.statIter);
    document["output"]["outIter"].read(param.outIter);
    document["output"]["cpIter"].read(param.cpIter);    
    
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
