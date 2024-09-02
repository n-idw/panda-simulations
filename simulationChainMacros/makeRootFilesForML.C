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
    // RDataFrame dfSim ("pndsim", prefix + "_sim.root" );
    RDataFrame dfDigi("pndsim", prefix + "_digi.root");

    auto dfTest = dfDigi.Range(10);

    // // Take the sim parameters from the data frame
    // auto mcTrackPx            = dfSim.Take<RVecD>("MCTrack.fPx"            );
    // auto mcTrackPy            = dfSim.Take<RVecD>("MCTrack.fPy"            );
    // auto mcTrackPz            = dfSim.Take<RVecD>("MCTrack.fPz"            );
    // auto mcTrackStartX        = dfSim.Take<RVecD>("MCTrack.fStartX"        );
    // auto mcTrackStartY        = dfSim.Take<RVecD>("MCTrack.fStartY"        );
    // auto mcTrackStartZ        = dfSim.Take<RVecD>("MCTrack.fStartZ"        );
    // auto mcTrackStartT        = dfSim.Take<RVecD>("MCTrack.fStartT"        );
    // auto mcTrackGeneratorFlag = dfSim.Take<RVecI>("MCTrack.fGeneratorFlags");
    // auto mcTrackPdgCode       = dfSim.Take<RVecI>("MCTrack.fPdgCode"       );
    // auto mcTrackPoints        = dfSim.Take<RVecI>("MCTrack.fPoints"        );

    // auto sttPointX       = dfSim.Take<RVecD>("STTPoint.fX"      );
    // auto sttPointY       = dfSim.Take<RVecD>("STTPoint.fY"      );
    // auto sttPointZ       = dfSim.Take<RVecD>("STTPoint.fZ"      );
    // auto sttPointT       = dfSim.Take<RVecD>("STTPoint.fTime"   );
    // auto sttPointPx      = dfSim.Take<RVecD>("STTPoint.fPx"     );
    // auto sttPointPy      = dfSim.Take<RVecD>("STTPoint.fPy"     );
    // auto sttPointPz      = dfSim.Take<RVecD>("STTPoint.fPz"     );
    // auto sttPointTrackID = dfSim.Take<RVecI>("STTPoint.fTrackID");
    
    // Take the digi parameters from the data frame
    // auto sttHitRefIndex   = dfDigi.Take<RVecI>("STTHit.fRefIndex"  );
    // auto sttHitDepCharge  = dfDigi.Take<RVecD>("STTHit.fDepCharge" );
    // auto sttHitTubeID     = dfDigi.Take<RVecI>("STTHit.fTubeID"    );
    // auto sttHitDetectorID = dfDigi.Take<RVecI>("STTHit.fDetectorID");
    // auto sttHitIsochrone  = dfDigi.Take<RVecD>("STTHit.fIsochrone" );
    // auto sttHitX          = dfDigi.Take<RVecD>("STTHit.fX"         );
    // auto sttHitY          = dfDigi.Take<RVecD>("STTHit.fY"         );
    // auto sttHitZ          = dfDigi.Take<RVecD>("STTHit.fZ"         );

    std::cout << std::endl;
    std::cout << "Write the data to ROOT file..." << std::endl;
    std::cout << "Output Directory: " << outputDir << "/mlData.root" << std::endl;

    TFile* outputFile = new TFile(Form("%s/mlData.root",outputDir.c_str()), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) 
    {
        std::cerr << "Failed to open destination file!" << std::endl;
        return 0;
    }

    TTree* cells     = new TTree("cells"    , "cells");
    TTree* hits      = new TTree("hits"     , "hits");
    TTree* truth     = new TTree("truth"    , "truth");
    TTree* particles = new TTree("particles", "particles");

    
    RVecD sttHitX, sttHitY, sttHitZ, depCharge, isochrone, energyloss;
    RVecI sttHitRefIndexVec, sttHitVolume_id, sttHitLayer_id, sttHitModule_id, skewed;
    
    TBranch* sttHitRefIndexBranch = hits->Branch("hit_id", &sttHitRefIndexVec);
    TBranch* sttHitXBranch        = hits->Branch("x", &sttHitX);
    TBranch* sttHitYBranch        = hits->Branch("y", &sttHitY);
    TBranch* sttHitZBranch        = hits->Branch("z", &sttHitZ);
    TBranch* volume_idBranch      = hits->Branch("volume_id", &sttHitVolume_id);
    TBranch* layer_idBranch       = hits->Branch("layer_id", &sttHitLayer_id);
    TBranch* module_idBranch      = hits->Branch("module_id", &sttHitModule_id);

    TBranch* sttHitRefIndexBranch = cells->Branch("hit_id", &sttHitRefIndexVec);
    TBranch* volume_idBranch      = cells->Branch("volume_id", &sttHitVolume_id);
    TBranch* layer_idBranch       = cells->Branch("layer_id", &sttHitLayer_id);
    TBranch* module_idBranch      = cells->Branch("module_id", &sttHitModule_id);
    TBranch* depChargeBranch      = cells->Branch("dep_charge", &depCharge);
    TBranch* isochroneBranch      = cells->Branch("isochrone", &isochrone);
    TBranch* energylossBranch     = cells->Branch("energy_loss", &energyloss);
    TBranch* skewedBranch         = cells->Branch("skewed", &skewed);

    auto fillHitsAndCells = [hits, cells, sttGeoData, &sttHitX, &sttHitY, &sttHitZ, &sttHitRefIndexVec, &sttHitVolume_id, &sttHitLayer_id, &sttHitModule_id](RVecD x, RVecD y, RVecD z, RVecI hit_id, RVecI volume_id, RVecI module_id)
    {
        sttHitX = x;
        sttHitY = y;
        sttHitZ = z;
        sttHitRefIndexVec = hit_id;
        sttHitVolume_id = volume_id;
        sttHitModule_id = module_id;

        for(RVecI::iterator moduleID = sttHitModule_id.begin(); moduleID != sttHitModule_id.end(); ++moduleID)
            sttHitLayer_id.push_back(stoi(sttGeoData.at(*moduleID-1).at(1)));

        hits->Fill();
        cells->Fill();
    };
    
    std::vector<std::string> hitParamNames = {"STTHit.fX", "STTHit.fY", "STTHit.fZ","STTHit.fRefIndex", "STTHit.fDetectorID", "STTHit.fTubeID"};

    dfTest.Foreach(fillHitsAndCells, hitParamNames);

    cells->Write();

    // Initiate the number of events
    // int nEvents = sttHitRefIndex->size();

    // for (Long64_t hitNum = 0; hitNum < nEvents; hitNum++)
    // {
    //     // Fill the branches
    // }
    
    outputFile->Write();
    outputFile->Close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds." << std::endl;

    return 1;
}