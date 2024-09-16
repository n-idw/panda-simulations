// Macro for running Panda simulation  with Geant3  or Geant4 (M. Al-Turany)
// This macro is supposed to run the full simulation of the panda detector
// to run the macro:
// root  sim_complete.C  or in root session root>.x  sim_complete.C
// to run with different options:(e.g more events, different momentum, Geant4)
// root  sim_complete.C"(100, "TGeant4", 2)"

int sim_complete(Int_t nEvents=50000, TString prefix="xixi", TString inputGen="XiXibar.DEC", Double_t pBeam=4.6, Int_t seed=42, Int_t TargetMode=0) {
    
    std::cout << "FLAGS: " << nEvents << "," << prefix << "," << inputGen << "," << pBeam << std::endl;
    std::cout << "SEED : " << seed << std::endl;
    std::cout << "Mode : " << TargetMode << std::endl;
    
    // Set Seed for Random Generator
    // gRandom->SetSeed();
    
    // OR from a specific seed.
    gRandom->SetSeed(seed);
    
    
    /* ************************************************************************
    * TString inputGen =
    * EvtGen -> "xxxxxxxx.dec" (parses dec-file for initial particle)
    * 		    OR "xxxxxxx.dec:initial_particle"
    *
    * DPM    -> "dpm_xxxxx"
    * FTF    -> "ftf_xxxxx"
    * PYTHIA -> "pythia_xxxxx" (future replacement of DPM/FTF for Bkg Simulation)
    *
    * BOX    -> "box:type(pdgcode,mult):p(min,max):tht(min,max):phi(min,max)"
    * PIPI   -> "pipi:cosTheta(min,max)"
    * LEP    -> "leplep:pid(value):gegm(value):cosTheta(min,max)"
    ************************************************************************ */

    TString parAsciiFile = "all.par";
    //TString prefix     = "";                // "" (auto from inputGen), "llbar_fwp", "evtcomplete";
    //TString inputGen   = "llbar_fwp.DEC";   // EvtGen: llbar_fwp.DEC, bkg_xixibar.DEC, etc.
    //TString inputGen   = "dpm";             // BkgGen: dpm, ftf, pythia8 (will be default)
    TString decayMode = "/export/home/jenny/XiXibar4_6_DecayByGeant/UserDecayConfig.C";

    //----- Create the Simulation Run Manager
    PndMasterRunSim *fRun = new PndMasterRunSim();
    fRun->SetInput(inputGen);
    fRun->SetName("TGeant4");
    fRun->SetParamAsciiFile(parAsciiFile);
    fRun->SetNumberOfEvents(nEvents);
    fRun->SetBeamMom(pBeam);
    fRun->SetUserDecay(decayMode);

    //----- Set Target Mode
    // TargetMode=0: point-like target (default)
    // TargetMode=1: Cluster-jet Target
    // TargetMode=2: Pellet Target
    // TargetMode=3: Pellet Stations
    // TargetMode=4: Cluster-jet target ??? with measured gas density profile in z
    fRun->SetTargetMode(TargetMode);
    
    //----- Initialization
    fRun->Setup(prefix);

    //----- Geometry
    fRun->CreateGeometry();

    //----- Event Generator
    fRun->SetGenerator();

    //----- Add Sim Tasks
    fRun->AddSimTasks();

    //----- Intialise & Run
    fRun->Init();
    fRun->Run(nEvents); 
    fRun->Finish();
    return 0;
}
