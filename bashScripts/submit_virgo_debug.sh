#!/bin/bash

#SBATCH --job-name=debug_mum5mup5
#SBATCH --output=%x_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=n_inde01@uni-muenster.de

#SBATCH --partition=debug
#SBATCH --time=0-00:05:00
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

# Activate the PandaRoot environment
source /lustre/panda/ninderw/PandaRoot/bin/config.sh -p

# Output directory for the simulation, cvs and log files
export OUTPUT_DIR=/lustre/panda/ninderw/simFiles

export CVS_GEN_FLAG="WithoutIdeal"	# With/Without IdealTrackFinder to Fill fParticles CSV
export NUM_GEN_EVENTS=10		# Number of events to be generated
export GENERATOR_NAME=DBoxGEN		# Specify the MC generator that should be used to generate the events
export BEAM_MOMENTUM=1.642		# Beam momentum in GeV/c

bash runSimChain_virgo.sh
