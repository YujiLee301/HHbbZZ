#ifndef H4LSelection_h
#define H4LSelection_h

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

#include <TLorentzVector.h>

namespace h4l {

struct CandidateLepton {
    TLorentzVector bareP4;
    int charge;
};

inline bool passesGhostRemoval(
    const std::array<CandidateLepton, 4>& leptons,
    double minDeltaR = 0.02
){
    for (std::size_t first = 0; first < leptons.size(); ++first){
        for (std::size_t second = first + 1; second < leptons.size(); ++second){
            if (leptons[first].bareP4.DeltaR(leptons[second].bareP4) <= minDeltaR){
                return false;
            }
        }
    }
    return true;
}

inline bool passesOppositeSignPairMass(
    const std::array<CandidateLepton, 4>& leptons,
    double minMass = 4.0
){
    for (std::size_t first = 0; first < leptons.size(); ++first){
        for (std::size_t second = first + 1; second < leptons.size(); ++second){
            if (leptons[first].charge + leptons[second].charge != 0){
                continue;
            }
            if ((leptons[first].bareP4 + leptons[second].bareP4).M() < minMass){
                return false;
            }
        }
    }
    return true;
}

inline std::pair<unsigned int, unsigned int> orderZCandidates(
    unsigned int first,
    unsigned int second,
    const TLorentzVector& firstP4,
    const TLorentzVector& secondP4,
    double nominalZMass
){
    if (std::fabs(firstP4.M() - nominalZMass) <=
        std::fabs(secondP4.M() - nominalZMass)){
        return std::make_pair(first, second);
    }
    return std::make_pair(second, first);
}

inline bool isBetterZZCandidate(
    double candidateZ1Distance,
    double candidateZ2PtSum,
    double currentZ1Distance,
    double currentZ2PtSum
){
    return candidateZ1Distance < currentZ1Distance ||
        (candidateZ1Distance == currentZ1Distance &&
         candidateZ2PtSum > currentZ2PtSum);
}

inline bool passesZ1Mass(const TLorentzVector& z1, double minMass){
    return z1.M() > minMass;
}

inline bool passesFourLeptonMass(
    const TLorentzVector& z1,
    const TLorentzVector& z2,
    double minMass
){
    return (z1 + z2).M() > minMass;
}

}  // namespace h4l

#endif
