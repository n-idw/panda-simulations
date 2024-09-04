#include "makeRootFilesForML.h"

/// @brief Returns the content of a CSV file as a 2d string vector
/// @param geoCsvFileName Name of the CSV file to read
/// @return 2d string vector containing the content of the CSV file
std::vector<std::vector<std::string>> readGeomFromCsv(const char* geoCsvFileName)
{
    // Open the CSV file
    std::ifstream geoCsvFile(geoCsvFileName);

    // Check if the file was opened successfully
    if (!geoCsvFile.is_open()) 
    {
        throw std::runtime_error(Form("Could not open file: %s", geoCsvFileName));
    }
    
    // Declare the 2d vector to store the CSV data
    std::vector<std::vector<std::string>> geoData;
    
    // Declare a string to store the current line of the CSV file
    std::string line;

    // Skip the first line of the file as it only contains the column names
    std::getline(geoCsvFile, line);

    // Iterate over the lines of the file
    while (std::getline(geoCsvFile, line)) {
        // Declare a vector to store the field of the current line
        std::vector<std::string> row;
        
        // Initialize a string stream with the current line
        std::istringstream s(line);
        
        // Declare a string to store the current field
        std::string field;
        
        // Iterate over the fields of the current line which are separated by commas and store them in the row vector
        while (getline(s, field, ',')) {
            row.push_back(field);
        }

        // Store the row vector in the geoData vector
        geoData.push_back(row);
    }

    // Close the CSV file
    geoCsvFile.close();

    return geoData;
}

int makeRootFilesForML(std::string prefix, std::string outputDir)
{
    // Start timer to measure the elapsed time
    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "Get the STT geometry data..." << std::endl;

    // 2d string vector containing the STT geometry data
    auto sttGeoData = readGeomFromCsv("/home/nikin105/mlProject/data/detectorGeometries/tubePos.csv");

    // Define an vector that contains the tube polarity (0 = straight, 1 = +3° skewed, -1 = -3° skewed)
    std::vector<int> tubePolarity;

    // Fill the vector according to the values from the geometry data
    for(int tube = 0; tube < sttGeoData.size(); tube++)
    {
        // Check if the tube is straight, skewed by +3° or skewed by -3°
        if(std::stoi(sttGeoData.at(tube).at(9)) == 0)
            tubePolarity.push_back(0);
        else
            tubePolarity.push_back(std::stoi(sttGeoData.at(tube).at(9))/abs(std::stoi(sttGeoData.at(tube).at(9)))); // +-1 = a / abs(a)
    }

    std::cout << std::endl;
    std::cout << "Read the simulation and digitization files..." << std::endl;
    std::cout << "Simulation File: " << prefix + "_sim.root" << std::endl;
    std::cout << "Digitalization File: " << prefix + "_digi.root" << std::endl;

    // ROOT::EnableImplicitMT();

    // Create a ROOT data frame for the simulation and digitization files
    RDataFrame dfSim ("pndsim", prefix + "_sim.root" );
    RDataFrame dfDigi("pndsim", prefix + "_digi.root");

    std::cout << std::endl;
    std::cout << "Write the data to ROOT file..." << std::endl;
    std::cout << "Output Directory: " << outputDir << "/mlData.root" << std::endl;

    TFile* outputFile = new TFile(Form("%s/mlData.root",outputDir.c_str()), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) 
    {
        std::cerr << "Failed to open destination file!" << std::endl;
        return 0;
    }

    // Output
    RVecD v_x, v_y, v_z, v_dep_charge, v_isochrone, v_energy_loss, 
          v_tx, v_ty, v_tz, v_tT, v_tpx, v_tpy, v_tpz,
          v_vx, v_vy, v_vz, v_px, v_py, v_pz, v_start_time;
    
    RVecI v_hit_id, v_volume_id, v_module_id, v_layer_id, v_sector_id,
           v_skewed, v_particle_id, v_nhits, v_pdgcode, v_primary, v_particle_num;

    TTree* cells     = new TTree("cells"    , "cells");
    TTree* hits      = new TTree("hits"     , "hits");
    TTree* truth     = new TTree("truth"    , "truth");
    TTree* particles = new TTree("particles", "particles");
    
    TBranch* b_hits_hit_id    = hits->Branch("hit_id"   , &v_hit_id);
    TBranch* b_hits_x         = hits->Branch("x"        , &v_x);
    TBranch* b_hits_y         = hits->Branch("y"        , &v_y);
    TBranch* b_hits_z         = hits->Branch("z"        , &v_z);
    TBranch* b_hits_volume_id = hits->Branch("volume_id", &v_volume_id);
    TBranch* b_hits_layer_id  = hits->Branch("layer_id" , &v_layer_id);
    TBranch* b_hits_module_id = hits->Branch("module_id", &v_module_id);

    TBranch* b_cells_hit_id      = cells->Branch("hit_id"     , &v_hit_id);
    TBranch* b_cells_dep_charge  = cells->Branch("dep_charge" , &v_dep_charge);
    TBranch* b_cells_energy_loss = cells->Branch("energy_loss", &v_energy_loss);
    TBranch* b_cells_volume_id   = cells->Branch("volume_id"  , &v_volume_id);
    TBranch* b_cells_layer_id    = cells->Branch("layer_id"   , &v_layer_id);
    TBranch* b_cells_module_id   = cells->Branch("module_id"  , &v_module_id);
    TBranch* b_cells_sector_id   = cells->Branch("sector_id"  , &v_sector_id);
    TBranch* b_cells_isochrone   = cells->Branch("isochrone"  , &v_isochrone);
    TBranch* b_cells_skewed      = cells->Branch("skewed"     , &v_skewed);

    TBranch* b_truth_hit         = truth->Branch("hit_id"     , &v_hit_id);
    TBranch* b_truth_x           = truth->Branch("tx"         , &v_tx);
    TBranch* b_truth_y           = truth->Branch("ty"         , &v_ty);
    TBranch* b_truth_z           = truth->Branch("tz"         , &v_tz);
    TBranch* b_truth_t           = truth->Branch("tT"         , &v_tT);
    TBranch* b_truth_px          = truth->Branch("tpx"        , &v_tpx);
    TBranch* b_truth_py          = truth->Branch("tpy"        , &v_tpy);
    TBranch* b_truth_pz          = truth->Branch("tpz"        , &v_tpz);
    TBranch* b_truth_particle_id = truth->Branch("particle_id", &v_particle_id);

    TBranch* b_particles_particle_id = particles->Branch("particle_id", &v_particle_num);
    TBranch* b_vx                    = particles->Branch("vx"         , &v_vx);
    TBranch* b_vy                    = particles->Branch("vy"         , &v_vy);
    TBranch* b_vz                    = particles->Branch("vz"         , &v_vz);
    TBranch* b_px                    = particles->Branch("px"         , &v_px);
    TBranch* b_py                    = particles->Branch("py"         , &v_py);
    TBranch* b_pz                    = particles->Branch("pz"         , &v_pz);
    TBranch* b_nhits                 = particles->Branch("nhits"      , &v_nhits);
    TBranch* b_pdgcode               = particles->Branch("pdgcode"    , &v_pdgcode);
    TBranch* b_start_time            = particles->Branch("start_time" , &v_start_time);
    TBranch* b_primary               = particles->Branch("primary"    , &v_primary);

    int numProcSimEvents  = 0;
    int numProcDigiEvents = 0;

    auto fillSimTrees = 
    [
        truth, particles, &numProcSimEvents,
        &v_hit_id, &v_tx, &v_ty, &v_tz, &v_tT, &v_tpx, &v_tpy, &v_tpz, &v_particle_id, &v_particle_num,
        &v_vx, &v_vy, &v_vz, &v_px, &v_py, &v_pz, &v_nhits, &v_pdgcode, &v_start_time, &v_primary
    ]
    (
        RVecD iv_tx, RVecD iv_ty, RVecD iv_tz, RVecD iv_tT, RVecD iv_tpx, RVecD iv_tpy, RVecD iv_tpz, RVecI iv_particle_id,
        RVecD iv_vx, RVecD iv_vy, RVecD iv_vz, RVecD iv_px, RVecD iv_py, RVecD iv_pz, RVecI iv_nhits, RVecI iv_pdgcode, RVecD iv_start_time, RVecI iv_generator_flag
    )
    {
        // Truth
        v_hit_id.clear();
        for (int i = 0; i < iv_tx.size(); i++)
            v_hit_id.push_back(i);

        v_tx          = iv_tx;
        v_ty          = iv_ty;
        v_tz          = iv_tz;
        v_tT          = iv_tT;
        v_tpx         = iv_tpx;
        v_tpy         = iv_tpy;
        v_tpz         = iv_tpz;
        v_particle_id = iv_particle_id;

        // Particles
        v_particle_num.clear();
        for (int i = 0; i < iv_tx.size(); i++)
            v_particle_num.push_back(i);

        v_vx         = iv_vx;
        v_vy         = iv_vy;
        v_vz         = iv_vz;
        v_px         = iv_px;
        v_py         = iv_py;
        v_pz         = iv_pz;
        v_nhits      = iv_nhits;
        v_pdgcode    = iv_pdgcode;
        v_start_time = iv_start_time;

        v_primary.clear();
        for(RVecI::iterator generator_flag = iv_generator_flag.begin(); generator_flag != iv_generator_flag.end(); ++generator_flag)
        {
            if(*generator_flag == 0)
                v_primary.push_back(0);
            else
                v_primary.push_back(1);
        }

        truth    ->Fill();
        particles->Fill();

        numProcSimEvents++;
        if (numProcSimEvents % 1000 == 0)
            std::cout << "Processed " << numProcSimEvents << " simulation events" << std::endl;
    };

    auto fillDigiTrees = 
    [
        hits, cells, sttGeoData, &numProcDigiEvents,
        &v_hit_id, &v_x, &v_y, &v_z, &v_volume_id, &v_layer_id, &v_module_id,
        &v_dep_charge, &v_energy_loss, &v_sector_id, &v_isochrone, &v_skewed
    ]
    (
        RVecI iv_hit_id, RVecD iv_x, RVecD iv_y, RVecD iv_z, RVecI iv_volume_id, RVecI iv_module_id,
        RVecD iv_dep_charge, RVecD iv_isochrone
    )
    {
        // Hits
        v_hit_id    = iv_hit_id;
        v_x         = iv_x;
        v_y         = iv_y;
        v_z         = iv_z;
        v_volume_id = iv_volume_id;
        v_module_id = iv_module_id;

        v_layer_id.clear();
        for(RVecI::iterator moduleID = v_module_id.begin(); moduleID != v_module_id.end(); ++moduleID)
            v_layer_id.push_back(stoi(sttGeoData.at(*moduleID-1).at(1)));

        // Cells
        v_dep_charge  = iv_dep_charge;
        v_energy_loss = iv_dep_charge / 1e6;

        v_sector_id.clear();
        for(RVecI::iterator moduleID = v_module_id.begin(); moduleID != v_module_id.end(); ++moduleID)
            v_sector_id.push_back(stoi(sttGeoData.at(*moduleID-1).at(2)));

        v_isochrone = iv_isochrone;

        v_skewed.clear();
        for(RVecI::iterator moduleID = v_module_id.begin(); moduleID != v_module_id.end(); ++moduleID)
            v_skewed.push_back(stoi(sttGeoData.at(*moduleID-1).at(10)));

        hits     ->Fill();
        cells    ->Fill();

        numProcDigiEvents++;
        if (numProcDigiEvents % 1000 == 0)
            std::cout << "Processed " << numProcDigiEvents << " digitization events" << std::endl;
    };
    
    std::vector<std::string> simParamNames = 
    {
        "STTPoint.fX", "STTPoint.fY", "STTPoint.fZ", "STTPoint.fTime", "STTPoint.fPx", "STTPoint.fPy", "STTPoint.fPz", "STTPoint.fTrackID",
        "MCTrack.fStartX", "MCTrack.fStartY", "MCTrack.fStartZ", "MCTrack.fPx", "MCTrack.fPy", "MCTrack.fPz", "MCTrack.fPoints", "MCTrack.fPdgCode", "MCTrack.fStartT", "MCTrack.fGeneratorFlags"
    };

    std::vector<std::string> digiParamNames = 
    {
        "STTHit.fRefIndex", "STTHit.fX", "STTHit.fY", "STTHit.fZ", "STTHit.fDetectorID", "STTHit.fTubeID",
        "STTHit.fDepCharge", "STTHit.fIsochrone",
    };

    std::cout << "Filling the digi trees..." << std::endl;
    dfDigi.Foreach(fillDigiTrees, digiParamNames);

    std::cout << "Filling the sim trees..." << std::endl;
    dfSim .Foreach(fillSimTrees , simParamNames);

    cells    ->Write();
    hits     ->Write();
    truth    ->Write();
    particles->Write();
    
    outputFile->Write();
    outputFile->Close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds." << std::endl;

    return 1;
}