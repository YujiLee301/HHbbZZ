#ifndef ZXCRSelection_h
#define ZXCRSelection_h

#include "H4LSelection.h"
#include <algorithm>
#include <vector>

namespace h4l {

struct ZXLepton {
    TLorentzVector bareP4, p4;
    int pdgId, sourceIndex, lepIndex, genIndex;
    bool tight;
};

struct ZXCandidate {
    // 0: none, 1: OS Z1+l 2P1F, 2: OS 2P2F, 3: OS 3P1F.
    int region = 0;
    std::vector<unsigned int> leptons;
    double z1Distance = 0., z2PtSum = 0.;
};

inline bool passesZXPairCuts(const std::vector<ZXLepton>& leptons,
                            const std::vector<unsigned int>& indices){
    for (std::size_t i = 0; i < indices.size(); ++i){
        for (std::size_t j = i + 1; j < indices.size(); ++j){
            const auto& first = leptons[indices[i]];
            const auto& second = leptons[indices[j]];
            if (first.bareP4.DeltaR(second.bareP4) <= 0.02) return false;
            if (first.pdgId * second.pdgId < 0 &&
                (first.bareP4 + second.bareP4).M() < 4.) return false;
        }
    }
    return true;
}

// Port of the OS selectors in yhbbzz/zxcr.py. Kinematic thresholds are
// supplied by the host framework; a tight Z1 is mandatory in every region.
inline ZXCandidate selectZXCandidate(const std::vector<ZXLepton>& leptons,
                                     double met, double zMass, double zMin,
                                     double zMax, double z1Min, double m4lMin){
    ZXCandidate best;
    struct Pair { unsigned int first, second; TLorentzVector p4; };
    std::vector<Pair> pairs;
    for (unsigned int i = 0; i < leptons.size(); ++i){
        for (unsigned int j = i + 1; j < leptons.size(); ++j){
            if (leptons[i].pdgId + leptons[j].pdgId != 0) continue;
            const auto p4 = leptons[i].p4 + leptons[j].p4;
            if (p4.M() > zMin && p4.M() < zMax) pairs.push_back({i, j, p4});
        }
    }
    bool hasTightZZ = false;
    for (const auto& z1 : pairs){
        if (!leptons[z1.first].tight || !leptons[z1.second].tight) continue;
        const double distance = std::fabs(z1.p4.M() - zMass);
        if (leptons.size() == 3){
            // OS fake-rate denominator: exactly three loose leptons, MET<25,
            // |mZ1-mZ|<7, and a failing third lepton.
            if (!std::isfinite(met) || met < 0. || met >= 25. || distance >= 7.) continue;
            if (std::max(leptons[z1.first].p4.Pt(), leptons[z1.second].p4.Pt()) <= 20. ||
                std::min(leptons[z1.first].p4.Pt(), leptons[z1.second].p4.Pt()) <= 10.) continue;
            const unsigned int probe = 3 - z1.first - z1.second;
            const std::vector<unsigned int> indices = {z1.first, z1.second, probe};
            if (leptons[probe].tight || !passesZXPairCuts(leptons, indices)) continue;
            if (!best.region || distance < best.z1Distance){
                best.region = 1;
                best.leptons = indices;
                best.z1Distance = distance;
            }
            continue;
        }
        if (z1.p4.M() <= z1Min) continue;
        for (const auto& z2 : pairs){
            if (z1.first == z2.first || z1.first == z2.second ||
                z1.second == z2.first || z1.second == z2.second) continue;
            std::vector<unsigned int> indices = {z1.first, z1.second, z2.first, z2.second};
            if (!passesZXPairCuts(leptons, indices)) continue;
            std::array<double, 4> pts;
            for (unsigned int i = 0; i < 4; ++i) pts[i] = leptons[indices[i]].p4.Pt();
            std::sort(pts.rbegin(), pts.rend());
            if (pts[0] <= 20. || pts[1] <= 10.) continue;
            if ((z1.p4 + z2.p4).M() <= m4lMin) continue;
            bool smartCut = true;
            if (std::abs(leptons[z1.first].pdgId) == std::abs(leptons[z2.first].pdgId)){
                // Inspect both unfiltered alternative pairings, including
                // pairs below zMin, as required by the reference smart cut.
                for (unsigned int k = 2; k < 4; ++k){
                    const auto& a = leptons[indices[0]];
                    const auto& b = leptons[indices[k]];
                    if (a.pdgId + b.pdgId != 0) continue;
                    const auto za = a.p4 + b.p4;
                    const auto zb = leptons[indices[1]].p4 + leptons[indices[5-k]].p4;
                    const auto ordered = orderZCandidates(0, 1, za, zb, zMass);
                    const auto& nearZ = ordered.first == 0 ? za : zb;
                    const auto& farZ = ordered.first == 0 ? zb : za;
                    if (std::fabs(nearZ.M() - zMass) < distance && farZ.M() < zMin) smartCut = false;
                }
            }
            if (!smartCut) continue;
            const int failed = int(!leptons[z2.first].tight) + int(!leptons[z2.second].tight);
            if (failed == 0){
                // Only the nominal mass-ordered tight pairing defines SR.
                if (distance <= std::fabs(z2.p4.M() - zMass)) hasTightZZ = true;
                continue;
            }
            const double ptSum = leptons[z2.first].p4.Pt() + leptons[z2.second].p4.Pt();
            if (best.region && !isBetterZZCandidate(distance, ptSum, best.z1Distance, best.z2PtSum)) continue;
            // The failing probe is L3; ties are resolved by decreasing pT.
            if ((leptons[indices[2]].tight && !leptons[indices[3]].tight) ||
                (leptons[indices[2]].tight == leptons[indices[3]].tight &&
                 leptons[indices[3]].p4.Pt() > leptons[indices[2]].p4.Pt())) std::swap(indices[2], indices[3]);
            best.region = failed == 1 ? 3 : 2;
            best.leptons = indices;
            best.z1Distance = distance;
            best.z2PtSum = ptSum;
        }
    }
    // Fully tight four-lepton candidates take priority even before jet cuts.
    return hasTightZZ ? ZXCandidate() : best;
}

} // namespace h4l
#endif
