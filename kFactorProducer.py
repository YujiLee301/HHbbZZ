# kFactorProducer.py
import os
import sys
import ROOT
import numpy as np
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

class kFactorProducer(Module):
    """
    Compute k-factors for ZZ production using GenAnalysis for generator information.
    qqZZ: NNLO/NLO
    ggZZ: NNLO/LO
    """
    
    def __init__(self, year, sample_path="", kfactor_dir=None):
        self.year = year
        self.sample_path = sample_path
        
        if kfactor_dir is None:
            cmssw_base = os.environ.get("CMSSW_BASE", "")
            kfactor_dir = os.path.join(
                cmssw_base,
                "src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/data/kFactor"
            )
        
        self.kfactor_dir = kfactor_dir
        
        print(f"\n{'='*60}")
        print(f"[kFactorProducer] Initializing...")
        print(f"[kFactorProducer] Year: {year}")
        print(f"[kFactorProducer] Sample path: {sample_path}")
        print(f"{'='*60}\n")
        
        print(f"[kFactorProducer] kfactor_dir: {self.kfactor_dir}")
        print(f"[kFactorProducer] exists? {os.path.exists(self.kfactor_dir)}")
        
        full_path_lower = sample_path.lower().strip()

        # ggZZ continuum
        self.is_ggzz = "gluglutocontinto2z" in full_path_lower

        # qqZZ background:
        self.is_qqzz = ("zzto4l" in full_path_lower) and ("glugluhtozzto4l" not in full_path_lower)

        self.apply_kfactor = self.is_ggzz or self.is_qqzz
        
        print(f"[kFactorProducer] is_ggZZ: {self.is_ggzz}")
        print(f"[kFactorProducer] is_qqZZ: {self.is_qqzz}")
        print(f"[kFactorProducer] apply_kfactor: {self.apply_kfactor}")
        print(f"{'='*60}\n")
        
        # Load GenAnalysis
        self._load_genanalysis()
        
        # Initialize k-factor data
        self.ggzz_nnlo_file = None
        self.ggzz_nlo_file = None
        self.spkfactor_ggzz_nnlo = {}
        self.spkfactor_ggzz_nlo = {}
        self.ewk_table = None
        
        # Load k-factor files (only if needed)
        if self.apply_kfactor:
            self._load_kfactor_files()
        else:
            print(f"[kFactorProducer] Sample does not require k-factors, skipping file loading")
        
        print(f"[kFactorProducer] Initialization complete\n")
    
    def _load_genanalysis(self):
        """Load GenAnalysis library."""
        current_dir = os.path.dirname(os.path.abspath(__file__))
        lib_path = os.path.join(current_dir, "src", "GenAnalysis_cc.so")
        
        if not os.path.exists(lib_path):
            src_path = os.path.join(current_dir, "src", "GenAnalysis.cc")
            if os.path.exists(src_path):
                ROOT.gROOT.ProcessLine(f".L {src_path}+")
            else:
                print(f"[kFactorProducer] ERROR: GenAnalysis source not found")
                sys.exit(1)
        else:
            ROOT.gSystem.Load(lib_path)
        
        self.gen_analyzer = ROOT.GenAnalysis()
    
    def _load_kfactor_files(self):
        """Load ROOT files and initialize splines/histograms."""
        # Load ggZZ k-factors
        if self.is_ggzz:
            self._load_ggzz_kfactors()
        
        # Load qqZZ k-factors (EWK table)
        if self.is_qqzz:
            self._load_qqzz_kfactors()
    
    def _load_ggzz_kfactors(self):
        """Load ggZZ NNLO and NLO k-factor splines."""
        nnlo_path = os.path.join(self.kfactor_dir, 
                             "Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root")
        nlo_path = os.path.join(self.kfactor_dir,
                            "Kfactor_Collected_ggHZZ_2l2l_NLO_NNPDF_NarrowWidth_13TeV.root")
    
                
        print(f"[ggZZ] nnlo_path = {nnlo_path}, exists = {os.path.exists(nnlo_path)}")
        print(f"[ggZZ] nlo_path  = {nlo_path}, exists = {os.path.exists(nlo_path)}")
        
    
        variations = ["Nominal", "PDFScaleDn", "PDFScaleUp", 
                  "QCDScaleDn", "QCDScaleUp", "AsDn", "AsUp",
                  "PDFReplicaDn", "PDFReplicaUp"]
    
        if os.path.exists(nnlo_path):
            f = ROOT.TFile.Open(nnlo_path, "READ")
            if f and not f.IsZombie():
                for var in variations:
                    spline = f.Get(f"sp_kfactor_{var}")
                    if spline:
                        self.spkfactor_ggzz_nnlo[var] = spline.Clone(f"{var}_NNLO")
                f.Close()
    
        if os.path.exists(nlo_path):
            f = ROOT.TFile.Open(nlo_path, "READ")
            if f and not f.IsZombie():
                for var in variations:
                    spline = f.Get(f"sp_kfactor_{var}")
                    if spline:
                        self.spkfactor_ggzz_nlo[var] = spline.Clone(f"{var}_NLO")
                f.Close()
    
    def _load_qqzz_kfactors(self):
        """Load qqZZ EWK correction table."""
        ewk_path = os.path.join(self.kfactor_dir, "ZZ_EwkCorrections.dat")
        if os.path.exists(ewk_path):
            self.ewk_table = self._read_ewk_table(ewk_path)
    
    def _read_ewk_table(self, filepath):
        """Read EWK corrections from text file."""
        table = []
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split()
                    if len(parts) >= 5:
                        table.append([float(parts[0]), float(parts[1]), 
                                    float(parts[2]), float(parts[3]), float(parts[4])])
        except:
            pass
        return table
    
    def _find_ewk_correction(self, sqrt_s_hat, t_hat, quark_type):
        """Find EWK correction for given sqrt_s_hat, t_hat, and quark flavor."""
        if not self.ewk_table:
            return 1.0
        
        n_t_bins = 200
        best_s_idx = 0
        best_s_diff = 0.8E+04
        
        for i in range(0, len(self.ewk_table), n_t_bins):
            s_diff = abs(sqrt_s_hat - self.ewk_table[i][0])
            if s_diff < best_s_diff:
                best_s_diff = s_diff
                best_s_idx = i
            else:
                break
        
        if sqrt_s_hat > 0.8E+04:
            best_s_idx = 39800
        
        best_t_idx = best_s_idx
        best_t_diff = 1.0E+09
        
        max_t_in_block = self.ewk_table[best_s_idx + n_t_bins - 1][1]
        
        if t_hat > max_t_in_block:
            best_t_idx = best_s_idx + n_t_bins - 1
        else:
            for k in range(best_s_idx, min(best_s_idx + n_t_bins, len(self.ewk_table))):
                t_diff = abs(t_hat - self.ewk_table[k][1])
                if t_diff < best_t_diff:
                    best_t_diff = t_diff
                    best_t_idx = k
                else:
                    break
        
        ewk_uc = self.ewk_table[best_t_idx][2]
        ewk_ds = self.ewk_table[best_t_idx][3]
        ewk_b = self.ewk_table[best_t_idx][4]
        
        if quark_type in [1, 3]:   # d, s
            return 1.0 + ewk_ds
        elif quark_type in [2, 4]:  # u, c
            return 1.0 + ewk_uc
        elif quark_type == 5:        # b
            return 1.0 + ewk_b
        return 1.0
    
    def beginJob(self):
        pass
    
    def endJob(self):
        pass
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        """Declare output branches."""
        self.out = wrappedOutputTree
        
        if self.is_ggzz:
            self.out.branch("ggZZ_kf_NNLO_Nominal", "F")
            self.out.branch("ggZZ_kf_NNLO_PDFScaleUp", "F")
            self.out.branch("ggZZ_kf_NNLO_PDFScaleDn", "F")
            self.out.branch("ggZZ_kf_NNLO_QCDScaleUp", "F")
            self.out.branch("ggZZ_kf_NNLO_QCDScaleDn", "F")
            self.out.branch("ggZZ_kf_NNLO_AsUp", "F")
            self.out.branch("ggZZ_kf_NNLO_AsDn", "F")
            self.out.branch("ggZZ_kf_NNLO_PDFReplicaUp", "F")
            self.out.branch("ggZZ_kf_NNLO_PDFReplicaDn", "F")
            self.out.branch("ggZZ_kf_NLO_Nominal", "F")
        
        if self.is_qqzz:
            self.out.branch("qqZZ_kf_QCD_NLO", "F")
            self.out.branch("qqZZ_kf_QCD_NNLO", "F") 
            self.out.branch("qqZZ_kf_EWK", "F")
            self.out.branch("qqZZ_kf_QCD_NLO_dPhi", "F")
            self.out.branch("qqZZ_kf_QCD_NLO_Pt", "F")
            self.out.branch("qqZZ_kf_total", "F")
    
    # ========== k-factor calculation ==========
    def _kfactor_qqzz_qcd_M(self, gen_mass, final_state, order):
        """
        Mass-dependent QCD k-factor for qqZZ.
        Order: 1=NLO, 2=NNLO
        order=1 -> NLO/LO, order=2 -> NNLO/LO
        """
        xsec_table = self._get_xsec_table_qqzz()
    
        if final_state not in [1, 2]:
            return 1.0
        
        fs_idx = final_state - 1
        nbins = len(xsec_table[fs_idx])
        mass_points = [b[0] for b in xsec_table[fs_idx]]
        min_mass, max_mass = mass_points[0], mass_points[-1]
        m = max(min_mass, min(gen_mass, max_mass))
        
        cbin = -1
        for ix in range(nbins - 1):
            if xsec_table[fs_idx][ix][0] <= m < xsec_table[fs_idx][ix+1][0]:
                cbin = ix
                break
        
        if cbin < 0:
            cbin = nbins - 1 if m >= max_mass else 0
        
        xsec_lo = xsec_table[fs_idx][cbin][1]
        xsec_nlo = xsec_table[fs_idx][cbin][2]
        xsec_nnlo = xsec_table[fs_idx][cbin][3]
    
        if xsec_lo <= 0:
            return 1.0
    
        if order == 1:
            return xsec_nlo / xsec_lo   # NLO/LO
        elif order == 2:
            return xsec_nnlo / xsec_lo  # NNLO/LO
        return 1.0
    
    def _get_xsec_table_qqzz(self):
        """Return qqZZ cross section table."""
        return [
            # final_state = 1 (4e/4mu/4tau)
            [
                [0, 1.0958, 1.6231, 2.0064],
                [25, 1.7050, 2.5568, 3.0055],
                [50, 0.7215, 1.0767, 1.2602],
                [75, 9.9359, 12.4965, 12.8891],
                [100, 1.8898, 2.3230, 2.4458],
                [125, 0.8069, 1.0625, 1.1824],
                [150, 0.4201, 0.5770, 0.6541],
                [175, 0.5738, 0.7550, 0.8332],
                [200, 0.6426, 0.8436, 0.9284],
                [225, 0.4609, 0.6155, 0.6831],
                [250, 0.3251, 0.4406, 0.4938],
                [275, 0.2345, 0.3212, 0.3584],
                [300, 0.1731, 0.2392, 0.2724],
                [325, 0.1303, 0.1818, 0.2088],
                [350, 0.1002, 0.1404, 0.1609],
                [375, 0.0781, 0.1100, 0.1260],
                [400, 0.0619, 0.0874, 0.0995],
                [425, 0.0495, 0.0702, 0.0811],
                [450, 0.0400, 0.0570, 0.0648],
                [475, 0.0327, 0.0468, 0.0534],
            ],
            # final_state = 2 (2e2mu/2e2tau/2mu2tau)
            [
                [0, 5.8066, 8.3841, 10.4881],
                [25, 4.6286, 6.8487, 8.3869],
                [50, 1.6151, 2.3943, 2.8561],
                [75, 19.1061, 24.0053, 25.1089],
                [100, 3.8682, 4.7521, 5.1476],
                [125, 1.6456, 2.1684, 2.3851],
                [150, 0.8558, 1.1722, 1.3680],
                [175, 1.1630, 1.5267, 1.6855],
                [200, 1.2952, 1.6984, 1.8783],
                [225, 0.9262, 1.2358, 1.3679],
                [250, 0.6549, 0.8854, 0.9846],
                [275, 0.4710, 0.6440, 0.7311],
                [300, 0.3473, 0.4809, 0.5381],
                [325, 0.2620, 0.3631, 0.4135],
                [350, 0.2012, 0.2806, 0.3239],
                [375, 0.1573, 0.2210, 0.2593],
                [400, 0.1241, 0.1746, 0.2097],
                [425, 0.0993, 0.1399, 0.1663],
                [450, 0.0804, 0.1141, 0.1353],
                [475, 0.0654, 0.0938, 0.1048],
            ]
        ]
    
    def _kfactor_qqzz_qcd_dPhi(self, abs_dphi, final_state):
        """dPhi-dependent QCD k-factor for qqZZ."""
        if final_state == 1:
            bins = [
                (0.0, 0.1, 1.5158), (0.1, 0.2, 1.4963), (0.2, 0.3, 1.4955),
                (0.3, 0.4, 1.4833), (0.4, 0.5, 1.4656), (0.5, 0.6, 1.4915),
                (0.6, 0.7, 1.4412), (0.7, 0.8, 1.4408), (0.8, 0.9, 1.4143),
                (0.9, 1.0, 1.4225), (1.0, 1.1, 1.4010), (1.1, 1.2, 1.4085),
                (1.2, 1.3, 1.3812), (1.3, 1.4, 1.3706), (1.4, 1.5, 1.3473),
                (1.5, 1.6, 1.3401), (1.6, 1.7, 1.3127), (1.7, 1.8, 1.2901),
                (1.8, 1.9, 1.2553), (1.9, 2.0, 1.2545), (2.0, 2.1, 1.2240),
                (2.1, 2.2, 1.1788), (2.2, 2.3, 1.1626), (2.3, 2.4, 1.1054),
                (2.4, 2.5, 1.0747), (2.5, 2.6, 1.0219), (2.6, 2.7, 0.9463),
                (2.7, 2.8, 0.8575), (2.8, 2.9, 0.7166), (2.9, 3.1416, 1.1328),
            ]
        elif final_state == 2:
            bins = [
                (0.0, 0.1, 1.5138), (0.1, 0.2, 1.5417), (0.2, 0.3, 1.4978),
                (0.3, 0.4, 1.5350), (0.4, 0.5, 1.4782), (0.5, 0.6, 1.5043),
                (0.6, 0.7, 1.5206), (0.7, 0.8, 1.5070), (0.8, 0.9, 1.4942),
                (0.9, 1.0, 1.4505), (1.0, 1.1, 1.4608), (1.1, 1.2, 1.4716),
                (1.2, 1.3, 1.4677), (1.3, 1.4, 1.4224), (1.4, 1.5, 1.3972),
                (1.5, 1.6, 1.3756), (1.6, 1.7, 1.3919), (1.7, 1.8, 1.3686),
                (1.8, 1.9, 1.3179), (1.9, 2.0, 1.3140), (2.0, 2.1, 1.2746),
                (2.1, 2.2, 1.2423), (2.2, 2.3, 1.2447), (2.3, 2.4, 1.1463),
                (2.4, 2.5, 1.1078), (2.5, 2.6, 1.0421), (2.6, 2.7, 0.9736),
                (2.7, 2.8, 0.8722), (2.8, 2.9, 0.7345), (2.9, 3.1416, 1.1632),
            ]
        else:
            return 1.1
        
        for low, high, k in bins:
            if low < abs_dphi <= high:
                return k
        return 1.1
    
    def _kfactor_qqzz_qcd_Pt(self, gen_pt, final_state):
        """pT-dependent QCD k-factor for qqZZ."""
        if final_state == 1:
            bins = [
                (0.0, 5.0, 0.6416), (5.0, 10.0, 1.0999), (10.0, 15.0, 1.2939),
                (15.0, 20.0, 1.3786), (20.0, 25.0, 1.4243), (25.0, 30.0, 1.4504),
                (30.0, 35.0, 1.4702), (35.0, 40.0, 1.4883), (40.0, 45.0, 1.5057),
                (45.0, 50.0, 1.5021), (50.0, 55.0, 1.5092), (55.0, 60.0, 1.5246),
                (60.0, 65.0, 1.5240), (65.0, 70.0, 1.5242), (70.0, 75.0, 1.5542),
                (75.0, 80.0, 1.5254), (80.0, 85.0, 1.5790), (85.0, 90.0, 1.5303),
                (90.0, 95.0, 1.5615), (95.0, 100.0, 1.5447), (100.0, float('inf'), 1.5722),
            ]
        elif final_state == 2:
            bins = [
                (0.0, 5.0, 0.7436), (5.0, 10.0, 1.1479), (10.0, 15.0, 1.3382),
                (15.0, 20.0, 1.4142), (20.0, 25.0, 1.4551), (25.0, 30.0, 1.4757),
                (30.0, 35.0, 1.4905), (35.0, 40.0, 1.5062), (40.0, 45.0, 1.5033),
                (45.0, 50.0, 1.5219), (50.0, 55.0, 1.5204), (55.0, 60.0, 1.5398),
                (60.0, 65.0, 1.5349), (65.0, 70.0, 1.5177), (70.0, 75.0, 1.5449),
                (75.0, 80.0, 1.5776), (80.0, 85.0, 1.5508), (85.0, 90.0, 1.5708),
                (90.0, 95.0, 1.5616), (95.0, 100.0, 1.5418), (100.0, float('inf'), 1.5849),
            ]
        else:
            return 1.1
        
        for low, high, k in bins:
            if low < gen_pt <= high:
                return k
        return 1.1
    
    def _compute_ggzz_kfactor(self, gen_mass):
        """Compute ggZZ k-factors from splines."""
        result = {}
        mass = float(gen_mass)
        
        variations = ["Nominal", "PDFScaleUp", "PDFScaleDn",
                      "QCDScaleUp", "QCDScaleDn", "AsUp", "AsDn",
                      "PDFReplicaUp", "PDFReplicaDn"]
        
        # NNLO splines
        for var in variations:
            key = f"NNLO_{var}"
            result[key] = 1.0
            if var in self.spkfactor_ggzz_nnlo:
                spline = self.spkfactor_ggzz_nnlo[var]
                try:
                    xmin = spline.GetXmin()
                    xmax = spline.GetXmax()
                    m = max(min(mass, xmax), xmin)
                    kf_value = spline.Eval(m)
                    result[key] = kf_value

                except Exception as e:
                    print(f"  Error evaluating NNLO spline for {var}: {e}")
        
        # NLO nominal
        result["NLO_Nominal"] = 1.0
        if "Nominal" in self.spkfactor_ggzz_nlo:
            try:
                spline = self.spkfactor_ggzz_nlo["Nominal"]
                xmin = spline.GetXmin()
                xmax = spline.GetXmax()
                m = max(min(mass, xmax), xmin)
                result["NLO_Nominal"] = spline.Eval(m)
            except:
                pass
        
        return result
    
    def _compute_qqzz_kfactor(self, gen_mass, gen_pt, gen_dphi, final_state, 
                              sqrt_s_hat=None, t_hat=None, quark_type=2):
        """
        Compute qqZZ k-factors.
        """
        result = {
            "QCD_NLO": 1.0,      # NLO/LO
            "QCD_NNLO": 1.0,      # NNLO/NLO
            "EWK": 1.0,
            "QCD_NLO_dPhi": 1.0,
            "QCD_NLO_Pt": 1.0,
            "total": 1.0
        }
        
        k_nlo_lo = self._kfactor_qqzz_qcd_M(gen_mass, final_state, 1)   # NLO/LO
        k_nnlo_lo = self._kfactor_qqzz_qcd_M(gen_mass, final_state, 2)  # NNLO/LO
        
        if k_nlo_lo > 0:
            result["QCD_NLO"] = k_nlo_lo
            result["QCD_NNLO"] = k_nnlo_lo / k_nlo_lo  # NNLO/NLO
        else:
            result["QCD_NNLO"] = 1.0
        
        result["QCD_NLO_dPhi"] = self._kfactor_qqzz_qcd_dPhi(gen_dphi, final_state)
        result["QCD_NLO_Pt"] = self._kfactor_qqzz_qcd_Pt(gen_pt, final_state)
        
        if (
            sqrt_s_hat is not None
            and sqrt_s_hat > 0
            and t_hat is not None
            and t_hat != -1
        ):
            result["EWK"] = self._find_ewk_correction(
                sqrt_s_hat,
                t_hat,
                quark_type
            )
        result["total"] = result["QCD_NNLO"] * result["EWK"]
        
        return result
    
    def analyze(self, event):
        """Process event and compute k-factors using GenAnalysis for inputs."""
        if not self.apply_kfactor:
            return True
        
        # Reset GenAnalysis for this event
        self.gen_analyzer.Initialize()
        
        # Fill generator particle information
        try:
            genparts = Collection(event, "GenPart")
            n_gen_parts = len(genparts)
            
            if n_gen_parts == 0:
                self._fill_default_values()
                return True
            
            for gp in genparts:
                mother_idx = -1
                if hasattr(gp, 'genPartIdxMother') and gp.genPartIdxMother >= 0:
                    mother_idx = gp.genPartIdxMother
                
                self.gen_analyzer.SetGenParts(
                    gp.pt, gp.eta, gp.phi, gp.mass,
                    gp.pdgId, gp.status, mother_idx
                )
        except:
            self._fill_default_values()
            return True
        
        # Fill jet information
        try:
            if hasattr(event, 'nGenJet') and event.nGenJet > 0:
                genjets = Collection(event, "GenJet")
                for gj in genjets:
                    hadronFlavour = gj.hadronFlavour if hasattr(gj, 'hadronFlavour') else 0
                    self.gen_analyzer.SetGenJets(
                        gj.pt, gj.eta, gj.phi, gj.mass, hadronFlavour
                    )
        except:
            pass
        
        # Run GenAnalysis
        try:
            self.gen_analyzer.SetGenVariables()
        except:
            self._fill_default_values()
            return True
        
        # Get number of leptons from GenAnalysis
        n_gen_leptons = 0
        try:
            n_gen_leptons = len(self.gen_analyzer.GENlep_pt)
        except:
            n_gen_leptons = 0
        
        if n_gen_leptons < 4:
            self._fill_default_values()
            return True
        
        # Compute k-factors
        if self.is_ggzz:
            gen_mass = self.gen_analyzer.GENmassZZ
            
            if gen_mass <= 0 or np.isnan(gen_mass):
                self._fill_default_values()
                return True
            kf = self._compute_ggzz_kfactor(gen_mass)
            for name, value in kf.items():
                self.out.fillBranch(f"ggZZ_kf_{name}", value)
        
        elif self.is_qqzz:
            gen_mass = self.gen_analyzer.GENmassZZ
            gen_pt = self.gen_analyzer.GENpTZZ
            gen_dphi = self.gen_analyzer.GEN_dPhiZZ
            final_state = self.gen_analyzer.GEN_final_state
            quark_type = self.gen_analyzer.GEN_quark_type
            sqrt_s_hat = self.gen_analyzer.GEN_sqrt_s_hat
            t_hat = self.gen_analyzer.GEN_t_hat
            
            if (
                gen_mass <= 0
                or gen_pt < 0
                or gen_dphi < 0
                or final_state not in [1, 2]
                or np.isnan(gen_mass)
            ):
                self._fill_default_values()
                return True
            
            kf = self._compute_qqzz_kfactor(
                gen_mass, gen_pt, gen_dphi, final_state,
                sqrt_s_hat, t_hat, quark_type
            )
            
            self.out.fillBranch("qqZZ_kf_QCD_NLO", kf["QCD_NLO"])
            self.out.fillBranch("qqZZ_kf_QCD_NNLO", kf["QCD_NNLO"])
            self.out.fillBranch("qqZZ_kf_EWK", kf["EWK"])
            self.out.fillBranch("qqZZ_kf_QCD_NLO_dPhi", kf["QCD_NLO_dPhi"])
            self.out.fillBranch("qqZZ_kf_QCD_NLO_Pt", kf["QCD_NLO_Pt"])
            self.out.fillBranch("qqZZ_kf_total", kf["total"])
            
        
        return True
    
    def _fill_default_values(self):
        """Fill default values when no valid 4l system found."""
        try:
            if self.is_ggzz:
                for var in ["NNLO_Nominal", "NNLO_PDFScaleUp", "NNLO_PDFScaleDn",
                           "NNLO_QCDScaleUp", "NNLO_QCDScaleDn", "NNLO_AsUp",
                           "NNLO_AsDn", "NNLO_PDFReplicaUp", "NNLO_PDFReplicaDn",
                           "NLO_Nominal"]:
                    self.out.fillBranch(f"ggZZ_kf_{var}", 1.0)
            elif self.is_qqzz:
                for name in ["qqZZ_kf_QCD_NLO", "qqZZ_kf_QCD_NNLO", "qqZZ_kf_EWK",
                            "qqZZ_kf_QCD_NLO_dPhi", "qqZZ_kf_QCD_NLO_Pt", "qqZZ_kf_total"]:
                    self.out.fillBranch(name, 1.0)
        except:
            pass


def create_kfactor_producer(year, sample_path, kfactor_dir=None):
    """
    Factory function to create kFactorProducer with proper initialization.
    """
    return kFactorProducer(year, sample_path, kfactor_dir)
