#!/usr/bin/env python3
import os
import sys
import argparse

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import createJMECorrector
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSFProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *

from PhysicsTools.NATModules.modules.electronSF import ElectronSF as electronSF_natlib
from PhysicsTools.NATModules.modules.eleScaleRes import eleScaleRes as eleScaleRes_natlib
from PhysicsTools.NATModules.modules.muonSF import MuonSF as muonSF_natlib
from PhysicsTools.NATModules.modules.muonScaleRes import muonScaleRes as muonScaleRes_natlib
from PhysicsTools.NATModules.modules.puWeightProducer import puWeightProducer as puWeightProducer_natlib

# Custom module imports
from H4Lmodule import *
from H4LCppModule import *
from JetSFMaker import *

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFile", default="", type=str, help="Input file name")
    parser.add_argument("-n", "--entriesToRun", default=100, type=int, help="Set  to 0 if need to run over all entries else put number of entries to run")
    parser.add_argument("-d", "--DownloadFileToLocalThenRun", default=True, type=bool, help="Download file to local then run")
    parser.add_argument("--NOsyst", default=False, action="store_true", help="Do not run systematics")
    parser.add_argument("--overwritePt", default=True, type=bool, help="overwright muon pt from muon scale and res corrections")
    return parser.parse_args()


def getListFromFile(filename):
    """Read file list from a text file."""
    with open(filename, "r") as file:
        return ["root://cms-xrd-global.cern.ch/" + line.strip() for line in file]


def main():
    args = parse_arguments()

    # Initial setup
    testfilelist = []
    modulesToRun = []
    isMC = True
    isFSR = False
    year = None
    data_tag = ""
    cfgFile = None
    jsonFileName = None
    sfFileName = None
    overwritePt = args.overwritePt

    entriesToRun = int(args.entriesToRun)
    DownloadFileToLocalThenRun = args.DownloadFileToLocalThenRun

    # Determine list of files to process
    if args.inputFile.endswith(".txt"):
        testfilelist = getListFromFile(args.inputFile)
    elif args.inputFile.endswith(".root"):
        testfilelist.append(args.inputFile)
    else:
        print("INFO: No input file specified. Using default file list.")
        testfilelist = getListFromFile("ExampleInputFileList.txt")
    print(("DEBUG: Input file list: {}".format(testfilelist)))
    if len(testfilelist) == 0:
        print("ERROR: No input files found. Exiting.")
        exit(1)

    """Determine the year and type (MC or Data) of input ROOT file:
    For data the string "/data/" is always there. So, we take this
    as handle to decide if the root file is MC or data.
    """
    first_file = testfilelist[0]
    isMC = "/data/" not in first_file

    print(first_file, "\n isMC = ", isMC)
    if "Summer22" in first_file or "Run2022" in first_file:
        """Summer22 and Run2022 for identification of 2022 MC and data respectiverly
        """
        year = 2022
        if "22EE" in first_file : 
            data_tag = "post_EE"
        else :
            data_tag = "pre_EE"
        cfgFile = "Input_2022.yml"
        jsonFileName = "golden_Json/Cert_Collisions2022_355100_362760_Golden.json"
        sfFileName = "DeepCSV_102XSF_V2.csv" # FIXME: Update for year 2022
        # json here: https://gitlab.cern.ch/cms-muonPOG/muonscarekit, from ref: https://muon-wiki.docs.cern.ch/code/ptcorr/
        
        # --- Electron SF ---
        if "pre_EE" in data_tag:
            eleSF = electronSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22/electron.json.gz")
            era = "2022Re-recoBCD"
        else:
            eleSF = electronSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22EE/electron.json.gz")
            era = "2022Re-recoE+PromptFG"
        set_name = "Electron-ID-SF"
        # --- Reco SF ---
        def reco_wp(pt):
            if pt < 20: return "RecoBelow20"
            elif pt <= 75: return "Reco20to75"
            else: return "RecoAbove75"
        eleSF.addCorrection(set_name, era, reco_wp, "sf",     "Reco_sf")
        eleSF.addCorrection(set_name, era, reco_wp, "sfup",   "Reco_sfUp")
        eleSF.addCorrection(set_name, era, reco_wp, "sfdown", "Reco_sfDown")
        # --- ID SF ---
        for wp in ["Loose", "Medium"]:
            eleSF.addCorrection(set_name, era, wp, "sf",     f"ID_{wp}_sf")
            eleSF.addCorrection(set_name, era, wp, "sfup",   f"ID_{wp}_sfUp")
            eleSF.addCorrection(set_name, era, wp, "sfdown", f"ID_{wp}_sfDown")
        modulesToRun.append(eleSF)
        
        # --- Muon SF ---
        if "post_EE" not in data_tag:
            muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22/muon_Z.json.gz")
        else:
            muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22EE/muon_Z.json.gz")
        muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "nominal", "MuonSF")
        muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "syst", "MuonSFsyst")
        modulesToRun.append(muSF)
        
        # --- Muon scale/resolution ---
        if "pre_EE" in data_tag:
            muonScale_json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22/muon_scalesmearing.json.gz"
        else:
            muonScale_json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22EE/muon_scalesmearing.json.gz"

        modulesToRun.append(muonScaleRes_natlib(muonScale_json, is_mc = isMC, overwritePt = overwritePt, minPt=3.))  
            
        # --- Electron scale/resolution ---
        EtDependent = True
        if "pre_EE" in data_tag:
            eleScale_json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22/electronSS_EtDependent.json.gz"
            scaleKey = "EGMScale_Compound_Ele_2022preEE"                
            smearKey = "EGMSmearAndSyst_ElePTsplit_2022preEE"
        else:
            eleScale_json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22EE/electronSS_EtDependent.json.gz"
            scaleKey = "EGMScale_Compound_Ele_2022postEE"
            smearKey = "EGMSmearAndSyst_ElePTsplit_2022postEE"

        modulesToRun.append(eleScaleRes_natlib(eleScale_json, scaleKey, smearKey, overwritePt))

    elif "Summer23" in first_file or "Run2023" in first_file:
        """Summer23 and Run2023 for identification of 2022 MC and data respectiverly
        """
        year = 2023
        if "23BPix" in first_file : 
            data_tag = "post_BPix"
        else : 
            data_tag = "pre_BPix"
        cfgFile = "Input_2023.yml"
        jsonFileName = "golden_Json/Cert_Collisions2022_355100_362760_Golden.json"
        sfFileName = "DeepCSV_102XSF_V2.csv" # FIXME: Update for year 2022
        
        # --- Electron SF ---
        if "pre_BPix" in data_tag:
            eleSF = electronSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2023_Summer23/electron.json.gz")
            era = "2023PromptC"
        else:
            eleSF = electronSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2023_Summer23BPix/electron.json.gz")
            era = "2023PromptD"
        set_name = "Electron-ID-SF"
        # --- Reco SF ---
        def reco_wp(pt):
            if pt < 20: 
                return "RecoBelow20"
            elif pt <= 75: 
                return "Reco20to75"
            else: 
                return "RecoAbove75"
        eleSF.addCorrection(set_name, era, reco_wp, "sf",     "Reco_sf")
        eleSF.addCorrection(set_name, era, reco_wp, "sfup",   "Reco_sfUp")
        eleSF.addCorrection(set_name, era, reco_wp, "sfdown", "Reco_sfDown")
        # --- ID SF ---
        for wp in ["Loose", "Medium"]:
            eleSF.addCorrection(set_name, era, wp, "sf",     f"ID_{wp}_sf")
            eleSF.addCorrection(set_name, era, wp, "sfup",   f"ID_{wp}_sfUp")
            eleSF.addCorrection(set_name, era, wp, "sfdown", f"ID_{wp}_sfDown")
        modulesToRun.append(eleSF)
        
        # --- Muon SF ---
        if "post_BPix" not in data_tag:
            muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2023_Summer23/muon_Z.json.gz")
        else:
            muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2023_Summer23BPix/muon_Z.json.gz")
        muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "nominal", "MuonSF")
        muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "syst", "MuonSFsyst")
        modulesToRun.append(muSF)
        
        # --- Muon scale/resolution ---
        if "pre_BPix" in data_tag:
            muon_json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2023_Summer23/muon_scalesmearing.json.gz"
        else:
            muon_json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2023_Summer23BPix/muon_scalesmearing.json.gz"
        modulesToRun.append(muonScaleRes_natlib(muon_json, is_mc = isMC, overwritePt = overwritePt, minPt=3.))
        
        # --- Electron scale/resolution
        EtDependent = True
        if "post_BPix" not in data_tag:
            json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2023_Summer23/electronSS_EtDependent.json.gz"
            scaleKey = "EGMScale_Compound_Ele_2023preBPIX"
            smearKey = "EGMSmearAndSyst_ElePTsplit_2023preBPIX"
        else:
            json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2023_Summer23BPix/electronSS_EtDependent.json.gz"
            scaleKey = "EGMScale_Compound_Ele_2023postBPIX"
            smearKey = "EGMSmearAndSyst_ElePTsplit_2023postBPIX"
        modulesToRun.append(eleScaleRes_natlib(json_path, scaleKey, smearKey, overwritePt))

    elif "UL18" in first_file or "UL2018" in first_file:
        """UL2018 for identification of 2018 UL data and UL18 for identification of 2018 UL MC
        """
        # ---Electron SF ---
        eleSF = electronSF_natlib(json_path)
        set_name = "UL-Electron-ID-SF"
        eleSF.addCorrection(set_name, "2018", "RecoAbove20", "sf", "RecoAbove20_sf")
        eleSF.addCorrection(set_name, "2018", "RecoBelow20", "sf", "RecoBelow20_sf")
        eleSF.addCorrection(set_name, "2018", "Medium", "sf", "ID_sf")
        modulesToRun.append(eleSF)
        # --- Muon SF ---
        muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2018_UL/muon_Z.json.gz")
        muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "nominal", "MuonSF")
        muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "syst", "MuonSFsyst")
        modulesToRun.append(muSF)
        # --- Muon/Electron scale&resolution ---
        year = 2018
        cfgFile = "Input_2018.yml"
        jsonFileName = "golden_Json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        modulesToRun.extend([muonScaleRes2018()])
        
        json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2018_UL/electronSS.json.gz"
        scaleKey = "EGMScale_Compound_Ele_2018"
        smearKey = "EGMSmearAndSyst_Ele_2018"
        modulesToRun.append(eleScaleRes_natlib(json_path, scaleKey, smearKey, isMC, True))
        
    elif "UL17" in first_file or "UL2017" in first_file:
        
        # ---Electron SF ---
        eleSF = electronSF_natlib(json_path)
        set_name = "UL-Electron-ID-SF"
        eleSF.addCorrection(set_name, "2017", "RecoAbove20", "sf", "RecoAbove20_sf")
        eleSF.addCorrection(set_name, "2017", "RecoBelow20", "sf", "RecoBelow20_sf")
        eleSF.addCorrection(set_name, "2017", "Medium", "sf", "ID_sf")
        modulesToRun.append(eleSF)
        # --- Muon SF ---
        muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2017_UL/muon_Z.json.gz")
        muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "nominal", "MuonSF")
        muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "syst", "MuonSFsyst")
        modulesToRun.append(muSF)
        # --- Muon/Electron scale&resolution ---
        year = 2017
        cfgFile = "Input_2017.yml"
        jsonFileName="golden_Json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        modulesToRun.extend([muonScaleRes2017()])
        
        json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2017_UL/electronSS.json.gz"
        scaleKey = "EGMScale_Compound_Ele_2017"
        smearKey = "EGMSmearAndSyst_Ele_2017"
        modulesToRun.append(eleScaleRes_natlib(json_path, scaleKey, smearKey, isMC, True))
        
    elif "UL16" in first_file or "UL2016" in first_file:
        
        # ---Electron SF ---
        eleSF = electronSF_natlib(json_path)
        set_name = "UL-Electron-ID-SF"
        tag = "2016preVFP" if "preVFP" in first_file or "APV" in first_file else "2016postVFP"
        eleSF.addCorrection(set_name, tag, "RecoAbove20", "sf", "RecoAbove20_sf")
        eleSF.addCorrection(set_name, tag, "RecoBelow20", "sf", "RecoBelow20_sf")
        eleSF.addCorrection(set_name, tag, "Medium", "sf", "ID_sf")
        modulesToRun.append(eleSF)
        # --- Muon SF ---
        if "APV" in first_file or "preVFP" in first_file:
            muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2016preVFP_UL/muon_Z.json.gz")
        else:
            muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2016postVFP_UL/muon_Z.json.gz")
        muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "nominal", "MuonSF")
        muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "syst", "MuonSFsyst")
        modulesToRun.append(muSF)
        # --- Muon/Electron scale&resolution ---
        year = 2016
        jsonFileName = "golden_Json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        modulesToRun.extend([muonScaleRes2016()])
        
        if "APV" in first_file or "preVFP" in first_file:
            # Electron energy scale and smearing corrections
            json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2016preVFP_UL/electronSS.json.gz"
            scaleKey = "EGMScale_Compound_Ele_2016preVFP"
            smearKey = "EGMSmearAndSyst_Ele_2016preVFP"
        else:
            # Electron energy scale and smearing corrections
            json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2016postVFP_UL/electronSS.json.gz"
            scaleKey = "EGMScale_Compound_Ele_2016postVFP"
            smearKey = "EGMSmearAndSyst_Ele_2016postVFP"
        modulesToRun.append(eleScaleRes_natlib(json_path, scaleKey, smearKey, isMC, True))
        
    modulesToRun.append(HZZAnalysisCppProducer(year,cfgFile, isMC, isFSR))
        
    print(("Input json file: {}".format(jsonFileName)))
    print(("Input cfg file: {}".format(cfgFile)))
    print(("isMC: {}".format(isMC)))
    print(("isFSR: {}".format(isFSR)))

    if isMC:
        if (not args.NOsyst):
            # FIXME: JES not used properly
            #jetmetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK4PFchs")
            #fatJetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK8PFPuppi")
            # btagSF = lambda: btagSFProducer("UL"+str(year), algo="deepjet",selectedWPs=['L','M','T','shape_corr'], sfFileName=sfFileName)
            #btagSF = lambda: btagSFProducer(era = "UL"+str(year), algo = "deepcsv")
            puidSF = lambda: JetSFMaker("%s" % year)
            #modulesToRun.extend([jetmetCorrector(), fatJetCorrector()])#, puidSF()
            # # modulesToRun.extend([jetmetCorrector(), fatJetCorrector(), btagSF(), puidSF()])

        # PU reweight
        if year == 2018: modulesToRun.extend([puAutoWeight_2018()])
        if year == 2017: modulesToRun.extend([puAutoWeight_2017()])
        if year == 2016: modulesToRun.extend([puAutoWeight_2016()])

        if year == 2022 :
            if "post_EE" not in data_tag :
                json = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/data/puWeights/puWeights_2022_Summer22.json.gz" % os.environ['CMSSW_BASE']
                key = "Collisions2022_355100_357900_eraBCD_GoldenJson"
            else :
                json = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/data/puWeights/puWeights_2022_Summer22EE.json.gz" % os.environ['CMSSW_BASE']
                key = "Collisions2022_359022_362760_eraEFG_GoldenJson"
            modulesToRun.insert(0, puWeightProducer_natlib(json, key))

        if year == 2023 :
            if "pre_BPix" in data_tag :
                json = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/data/puWeights/puWeights_2023_Summer23.json.gz" % os.environ['CMSSW_BASE']
                key = "Collisions2023_366403_369802_eraBC_GoldenJson"
            else :
                json = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/data/puWeights/puWeights_2023_Summer23BPix.json.gz" % os.environ['CMSSW_BASE']
                key = "Collisions2023_369803_370790_eraD_GoldenJson"
            modulesToRun.insert(0, puWeightProducer_natlib(json, key))


        # INFO: Keep the `fwkJobReport=False` to trigger `haddnano.py`
        #            otherwise the output file will have larger size then expected. Reference: https://github.com/cms-nanoAOD/nanoAOD-tools/issues/249
        p=PostProcessor(
            ".",
            testfilelist, 
            cut = "(Sum$(Muon_pt>3) + Sum$(Electron_pt>5) >= 4) && (nJet>=2)",
            branchsel = None,
            modules = modulesToRun, 
            provenance = True,
            fwkJobReport = True,
            haddFileName = "skimmed_nano.root", 
            maxEntries = entriesToRun, 
            prefetch = DownloadFileToLocalThenRun, 
            outputbranchsel = "keep_and_drop.txt"
        )
    else:
        #if (not args.NOsyst):
            # FIXME: JES not used properly
            #jetmetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK4PFchs")
            #fatJetCorrector = createJMECorrector(isMC=isMC, dataYear=year, jesUncert="All", jetType = "AK8PFPuppi")
            #modulesToRun.extend([jetmetCorrector(), fatJetCorrector()])

        p=PostProcessor(".",testfilelist, None, None, modules = modulesToRun, provenance=True, fwkJobReport=True,haddFileName="skimmed_nano.root", jsonInput=jsonFileName, maxEntries=entriesToRun, prefetch=DownloadFileToLocalThenRun, outputbranchsel="keep_and_drop_data.txt")

    p.run()


if __name__ == "__main__":
    main()
