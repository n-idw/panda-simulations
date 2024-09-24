// c++ includes
#include <iostream>

// ROOT includes
#include <TROOT.h>
#include <TRandom.h>

// FairRoot includes
#include <FairBoxGenerator.h>
#include <FairLogger.h>

// PandaRoot includes
#include <PndEmcGeoPar.h>
#include <PndEvtGenDirect.h>

int sim_complete(Int_t nEvents = 10, TString prefix = "", TString inputGen = "", Double_t pBeam = 1.642,
                 Int_t seed = 42)
{
    // Print the input parameters
    std::cout << std::endl;
    LOG(info) << "========== Input Parameters ==========";
    LOG(info) << "nEvents: " << nEvents;
    LOG(info) << "prefix: " << prefix.Data();
    LOG(info) << "inputGen: " << inputGen.Data();
    LOG(info) << "pBeam: " << pBeam;
    LOG(info) << "seed: " << seed;
    std::cout << std::endl;

    // Set the random seed
    gRandom->SetSeed(seed);

    //-------------------------------------------------------------------------//
    //                               User Settings                             //
    //-------------------------------------------------------------------------//

    TString parAsciiFile = "all.par"; // File that contains all detector parameters

    TString decayMode = "UserDecayConfig.C";

    //--------------------------------------------------------------//
    //             Create the Simulation Run Manager                //
    //--------------------------------------------------------------//

    PndMasterRunSim *fRun = new PndMasterRunSim();

    if (inputGen.Contains("dec")) // EvtGen Generator
    {
        LOG(info) << "Using the EvtGen generator...";
        std::ifstream decayFile(inputGen.Data());
        if (!decayFile)
        {
            LOG(ERROR) << "File " << inputGen.Data() << " not found!";
            return 1;
        }
        LOG(info) << "Input generator: " << inputGen.Data();
        fRun->SetUserDecay(decayMode);
        fRun->SetInput(inputGen);
    }
    else if (inputGen.Contains("DBoxGEN")) // Double Box Generator
    {
        LOG(info) << "Using Double BoxGenerator...";

        // 1st BoxGenerator
        FairBoxGenerator *boxGen1 = new FairBoxGenerator(13, 5); // 13=muon; 3122=Lambda; multiplicity
        boxGen1->SetName("boxGen for 5 mu-");                    // Name of the generator for the log file
        boxGen1->SetPRange(0.1, 1.5);                            // Momentum Range: 100 MeV to 1.5 GeV
        boxGen1->SetPhiRange(0., 360.);                          // Azimuth angle range [degree]
        boxGen1->SetThetaRange(22., 140.);                       // Polar angle in lab system range [degree], STT
        boxGen1->SetXYZ(0., 0., 0.); // Set vertex position in [cm]
        fRun->AddGenerator(boxGen1);

        // 2nd BoxGenerator
        FairBoxGenerator *boxGen2 = new FairBoxGenerator(-13, 5); // -13=anti-muon; -3122=anti-Lambda; multiplicity
        boxGen2->SetName("boxGen for 5 mu+");                     // Name of the generator for the log file
        boxGen2->SetPRange(0.1, 1.5);                             // Momentum range: 100 MeV to 1.5 GeV
        boxGen2->SetPhiRange(0., 360.);                           // Azimuth angle range [degree]
        boxGen2->SetThetaRange(22., 140.);                        // Polar angle in lab system range [degree], STT
        boxGen2->SetXYZ(0., 0., 0.); // Set vertex position in [cm]
        fRun->AddGenerator(boxGen2);
    }

    // other settings for the generator and propagator
    fRun->SetName("TGeant4");              // Engine for the detector simulation
    fRun->SetParamAsciiFile(parAsciiFile); // Name of the Ascii file that contains all the detector parameters (It
                                           // automatically searches in the $VMCWORKDIR)
    fRun->SetNumberOfEvents(nEvents);      // Number of events to be generated
    fRun->SetBeamMom(pBeam);               // Beam momentum of the anti-protons in GeV/c
    fRun->SetStoreTraj(kTRUE);             // Store the trajectories of the particles

    // -----  Initialization   ------------------------------------------------
    LOG(info) << "Setting up the run...";
    fRun->Setup(prefix);

    // -----   Geometry   -----------------------------------------------------
    LOG(info) << "Creating the geometry...";
    fRun->CreateGeometry();

    // -----   Event generator   ----------------------------------------------
    fRun->SetGenerator();

    // -----   Add tasks   ----------------------------------------------------
    LOG(info) << "Adding simulation tasks...";
    fRun->AddSimTasks();

    // -----   Initialise and run   --------------------------------------------
    LOG(info) << "Initialising and running the simulation...";

    fRun->Init(); // There still seem to be some issues with initializing some EMC functions when the container is
                  // loaded
    fRun->Run(nEvents);
    fRun->Finish();

    return 0;
}
