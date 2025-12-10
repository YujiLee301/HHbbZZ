#!/usr/bin/env python3
import os
import sys
import argparse

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import createJMECorrector
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSFProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *

# Corrections configuration
from corrections_config import get_corrections_modules, get_pu_weight_module

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
        
        # Get corrections from config
        modulesToRun.extend(get_corrections_modules(year, data_tag, first_file, isMC, overwritePt))

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
        
        # Get corrections from config
        modulesToRun.extend(get_corrections_modules(year, data_tag, first_file, isMC, overwritePt))

    elif "UL18" in first_file or "UL2018" in first_file:
        """UL2018 for identification of 2018 UL data and UL18 for identification of 2018 UL MC
        """
        year = 2018
        cfgFile = "Input_2018.yml"
        jsonFileName = "golden_Json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        
        # Get corrections from config
        modulesToRun.extend(get_corrections_modules(year, data_tag, first_file, isMC, overwritePt))
        
    elif "UL17" in first_file or "UL2017" in first_file:
        year = 2017
        cfgFile = "Input_2017.yml"
        jsonFileName="golden_Json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        
        # Get corrections from config
        modulesToRun.extend(get_corrections_modules(year, data_tag, first_file, isMC, overwritePt))
        
    elif "UL16" in first_file or "UL2016" in first_file:
        year = 2016
        jsonFileName = "golden_Json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
        sfFileName = "DeepCSV_102XSF_V2.csv"
        
        # Get corrections from config
        modulesToRun.extend(get_corrections_modules(year, data_tag, first_file, isMC, overwritePt))
        
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

        # PU weight for 2022 and 2023
        pu_weight_module = get_pu_weight_module(year, data_tag)
        if pu_weight_module is not None:
            modulesToRun.insert(0, pu_weight_module)


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
