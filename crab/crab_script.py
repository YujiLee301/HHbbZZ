#!/usr/bin/env python3
import os
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis
from PhysicsTools.NanoAODTools.postprocessing.analysis.nanoAOD_skim.wvAnalysisModule import wvAnalysisModule

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *

jetmetCorrector = createJMECorrector(isMC=True, dataYear=2018, jesUncert="All", redojec=True, jetType = "AK4PFchs")
fatJetCorrector = createJMECorrector(isMC=True, dataYear=2018, jesUncert="All", redojec=True, jetType = "AK8PFPuppi")


p=PostProcessor(".",inputFiles(),"","keep_and_drop.txt",modules=[jetmetCorrector(),fatJetCorrector(),wvAnalysisModule()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())

p.run()

print("DONE")
