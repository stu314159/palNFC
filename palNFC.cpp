#include "palabos3D.h"
#include "palabos3D.hh"
//#include "palabos2D.h"
//#include "palabos2D.hh"

#include <cstdlib>
#include <cmath>


using namespace plb;

#define DESCRIPTOR descriptors::D3Q19Descriptor

typedef double T;

std::string outdir(" ");

struct SimulationParameters
{
    // 2D vectors storing min/max for each dimension of Cartesian coordinate system.
    std::vector<T> xDomain;
    std::vector<T> yDomain;
    std::vector<T> zDomain;
    
    // [file name].STL with the STL-formatted description of an obstruction
    std::vector<std::string> obstFileName;
    
    // 0 = x+, 1 = y+, 2 = z+  This is a feature that NFC does not have
    int flowDirection;
    
    // Reynolds number of the flow and characteristic length.  
    // combined with fluid properties a representative average flow velocity is determined.
    T Re;
    T characteristicLength;
    
    // for NFC this is N_divs.  How many lattice points along the characteristic dimension
    plint resolution;
    
    // time step size.
    T dt;
    
    //number of initial iterations for smoothly increasing the inlet velocity
    plint startIter;
    
    // total iterations
    plint maxIter;
    
    // iterations between consol output of simulation statistics
    plint statIter;
    
    // iterations between data output (*.vti files)
    plint outIter;
    
    // iterations between creation of "check-point" files
    plint cpIter;
    
    
    // fluid density and kinematic viscosity
    T rho;
    T nu;
    
    // fluid ambient pressure (customarily zero)
    T ambientPressure;
    
    // Smagorinski constant for turbulent flow model
    T cSmago; 
    
    // parameters not set by the user
    
    // geometric length in x-, y-, and z-directions
    T lx, ly, lz;
    
    // distance between lattice points (uniform in all directions)
    T dx; 
    
    // relaxation constant
    T omega;
    
    // number of lattice points in x-, y-, and z-direction
    plint nx, ny, nz;
    
    // non-dimensional density (customarily equal to 1)
    T rho_LB;
    
    // non-dimensionalized ambient pressure (customarily 0)
    T ambientPressure_LB;
    
    // inlet velocity vector in physical units.
    // note NFC would only allow velocity in the z-direction
    Array<T,3> inletVelocity;
    
    // inlet velocity in LBM units
    Array<T,3> inletVelocity_LB;
    
    // Palabos data structures to represent sets of lattice points 
    // on the boundaries.
    Box3D inlet, outlet;
    Box3D lateral1;
    Box3D lateral2;
    Box3D lateral3;
    Box3D lateral4;
    
    plint smallEnvelopeWidth; // standard width
    
};

Array<T,3> getVelocity(Array<T,3> targetValue, plint iIter, plint startIter)
{
    return (targetValue * util::sinIncreasingFunction<T>(iIter,startIter));
}

void readInputParameters(std::string xmlInputFileName, SimulationParameters& param)
{
    // use the XMLreader object to collect SimulationParameters from
    // the designated input file.
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
    
    document["numerics"]["startIter"].read(param.startIter);
    
    document["numerics"]["maxIter"].read(param.maxIter);
    PLB_ASSERT(param.maxIter > 0);
    
    document["numerics"]["cSmago"].read(param.cSmago);
    PLB_ASSERT(param.cSmago >= 0);
    
    document["fluid"]["nu"].read(param.nu);
    PLB_ASSERT(param.nu > 0);
    
    document["fluid"]["rho"].read(param.rho);
    PLB_ASSERT(param.rho > 0);
    
    document["fluid"]["ambientPressure"].read(param.ambientPressure);
    
    document["output"]["statIter"].read(param.statIter);
    document["output"]["outIter"].read(param.outIter);
    document["output"]["cpIter"].read(param.cpIter);    
    
}

void defineOuterDomain(SimulationParameters& param)
{
    if (param.flowDirection == 0)  // flow in positive x-direction
    {
        param.inlet = Box3D(0,0,0,param.ny-1,0,param.nz-1);
        param.outlet = Box3D(param.nx - 1, param.nx-1, 0, param.ny-1,0,param.nz-1);
        param.lateral1 = Box3D(1,param.nx-2,0,0,0,param.nz-1); // "bottom"
        param.lateral2 = Box3D(1,param.nx-2,param.ny-1,param.ny-1,0,param.nz-1); // "top"
        param.lateral3 = Box3D(1,param.nx-2,1,param.ny-2,0,0); //"left"
        param.lateral4 = Box3D(1,param.nx-2,1,param.ny-2,param.nz-1,param.nz-1); // "right"    
    } else if (param.flowDirection == 1) // flow in positive y-direction
    {
        param.inlet = Box3D(0,param.nx-1,0,0,0,param.nz-1);
        param.outlet = Box3D(0,param.nx-1,param.ny-1,param.ny-1,0,param.nz-1);
        param.lateral1 = Box3D(0,param.nx-1,1,param.ny-2,0,0);
        param.lateral2 = Box3D(0,param.nx-1,1,param.ny-2,param.nz-1,param.nz-1);
        param.lateral3 = Box3D(0,0,1,param.ny-2,1,param.nz-2);
        param.lateral4 = Box3D(param.nx-1,param.nx-1,1,param.ny-2,1,param.nz-2);    
    } else // flow in positive z-direction
    {
        param.inlet = Box3D(0,param.nx-1,0,param.ny-1,0,0);
        param.outlet = Box3D(0,param.nx-1,0,param.ny-1,param.nz-1,param.nz-1);
        param.lateral1 = Box3D(0,0,0,param.ny-1,1,param.nz-2);
        param.lateral2 = Box3D(param.nx-1,param.nx-1,0,param.ny-1,1,param.nz-2);
        param.lateral3 = Box3D(1,param.nx-2,0,0,1,param.nz-2);
        param.lateral4 = Box3D(1,param.nx-2,param.ny-1,param.ny-1,1,param.nz-2);    
    }

}

void calculateDerivedSimulationParameters(SimulationParameters& param)
{
    // simple measure of the extent of the domain in each direction
    param.lx = param.xDomain[1] - param.xDomain[0];
    param.ly = param.yDomain[1] - param.yDomain[0];
    param.lz = param.zDomain[1] - param.zDomain[0];
    
    // recall: param.resolution is the number of lattice points along the
    // characteristicLength --> there are (param.resolution - 1) "gaps"
    param.dx = param.characteristicLength / (param.resolution - 1.0);
    
    // set up number of lattice points in each direction
    param.nx = util::roundToInt(param.lx/param.dx) + 1;
    param.ny = util::roundToInt(param.ly/param.dx) + 1;
    param.nz = util::roundToInt(param.lz/param.dx) + 1;
    
    // customary to set rho_LB to 1.0
    param.rho_LB = 1.0;
    
    // ambient pressure scaled to LB units.
    param.ambientPressure_LB = (1.0/param.rho)*(param.dt*param.dt/(param.dx*param.dx))*param.ambientPressure;
    
    T velMag = param.Re*param.nu/param.characteristicLength;
    if (param.flowDirection == 0)
    {
        param.inletVelocity[0]=velMag;
        param.inletVelocity[1]=0.; param.inletVelocity[2]=0.;
    } else if (param.flowDirection == 1)
    {
        param.inletVelocity[0] = 0; param.inletVelocity[2]=0;
        param.inletVelocity[1] = velMag;
    } else
    {
        param.inletVelocity[2] = velMag;
        param.inletVelocity[0] = 0; param.inletVelocity[1] = 0;
    }
    
    
    param.inletVelocity_LB = param.inletVelocity*(param.dt/param.dx);
    
    T nu_LB = param.nu * param.dt / (param.dx * param.dx);
    param.omega = 1.0 / (DESCRIPTOR<T>::invCs2 * nu_LB + 0.5);
    
    defineOuterDomain(param);
    
    param.smallEnvelopeWidth = 1; 
}

void createFluidBlocks(SimulationParameters& param, MultiBlockLattice3D<T,DESCRIPTOR>*& lattice)
{
    Dynamics<T,DESCRIPTOR> *dynamics = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(param.omega, param.cSmago);
    
    pcout << "Dynamics: Smagorinsky BGK. " << std::endl;
    
    Box3D fullDomain(0,param.nx-1,0,param.ny-1,0,param.nz-1);
    lattice = generateMultiBlockLattice<T,DESCRIPTOR>(fullDomain,dynamics->clone(),param.smallEnvelopeWidth).release();
    
    defineDynamics(*lattice, lattice->getBoundingBox(),dynamics->clone());
    delete dynamics;
    lattice->toggleInternalStatistics(false);



}

void outerDomainBoundaryConditions(SimulationParameters const& param, MultiBlockLattice3D<T,DESCRIPTOR> *lattice, OnLatticeBoundaryCondition3D<T,DESCRIPTOR> *bc)
{

    Array<T,3> velocity = getVelocity(param.inletVelocity_LB,0,param.startIter);
    
    pcout << "No-slip lateral boundaries. " << std::endl; //uh...I want *no-slip*
    lattice->periodicity().toggleAll(false);
    
    
    // specified velocity on the inlet
    bc->setVelocityConditionOnBlockBoundaries(*lattice,param.inlet,boundary::dirichlet);
    setBoundaryVelocity(*lattice,param.inlet,velocity);
    
    // no-slip walls
    bc->setVelocityConditionOnBlockBoundaries(*lattice,param.lateral1,boundary::dirichlet);
    bc->setVelocityConditionOnBlockBoundaries(*lattice,param.lateral2,boundary::dirichlet);
    bc->setVelocityConditionOnBlockBoundaries(*lattice,param.lateral3,boundary::dirichlet);
    bc->setVelocityConditionOnBlockBoundaries(*lattice,param.lateral4,boundary::dirichlet);
   
    
    Array<T,3> zeroVel(0.0,0.0,0.0);
    setBoundaryVelocity(*lattice,param.lateral1,zeroVel);
    setBoundaryVelocity(*lattice,param.lateral2,zeroVel);
    setBoundaryVelocity(*lattice,param.lateral3,zeroVel);
    setBoundaryVelocity(*lattice,param.lateral4,zeroVel);
   
    
    // outlet boundary
    if (param.flowDirection == 0)
    {
        bc->addPressureBoundary0P(param.outlet,*lattice);
    } else if (param.flowDirection == 1)
    {
        bc->addPressureBoundary1P(param.outlet,*lattice);
    } else 
    {
        bc->addPressureBoundary2P(param.outlet,*lattice);
    }
    setBoundaryDensity(*lattice,param.outlet,param.rho_LB);
    setBoundaryVelocity(*lattice,param.outlet,velocity);
}

template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice, SimulationParameters const& param, plint iter)
{
    T dx = param.dx;
    T dt = param.dt;
    ParallelVtkImageOutput3D<T> vtkOut(createFileName("vtk",iter,1),3,dx);
    vtkOut.writeData<3,float>(*computeVelocity(lattice),"velocity",dx/dt);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice),"velocityNorm",dx/dt);
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)),"vorticity",1./dt);
}

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir(outdir);
    
    plint numCores = global::mpi().getSize();
    pcout << "Number of MPI processes: " << numCores << std::endl;
    
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
    
    // read input parameters and compute derived parameters
    SimulationParameters param;
    readInputParameters(xmlInputFileName, param);
    calculateDerivedSimulationParameters(param);
    
    pcout << "omega = " << param.omega << std::endl;
    plint nnodes = param.nx*param.ny*param.nz;
    pcout << "Number of lattice points: " << nnodes << std::endl;
    
    // create the lattice
    MultiBlockLattice3D<T,DESCRIPTOR> *lattice = 0;    
    createFluidBlocks(param, lattice);
    
    // Boundary Conditions
    pcout << "Applying boundary conditions..." << std::endl;
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR> *bc = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    outerDomainBoundaryConditions(param,lattice,bc);
    delete bc;
    
    // incorporate the voxelized obstruction
    
    
    // initialize the lattice
    
    // cary out time stepping with periodic data output
    uint vtkNum = 0;
    
    for( plint iT=0; iT<param.maxIter; ++iT)
    {
        if(iT%param.outIter == 0 && iT>0)
        {
            pcout << "step " << iT << std::endl;
            writeVTK(*lattice,param,vtkNum); vtkNum++;
        }
        // execute a time iteration
        //lattice->collideAndStream();
        
    }
    
    delete lattice;
    
    return 0;
}
