#include "../interface/H4LSelection.h"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

namespace {

h4l::CandidateLepton makeLepton(
    double pt,
    double eta,
    double phi,
    int charge
){
    h4l::CandidateLepton lepton;
    lepton.bareP4.SetPtEtaPhiM(pt, eta, phi, 0.000511);
    lepton.charge = charge;
    return lepton;
}

TLorentzVector makeParticle(double pt, double phi, double mass){
    TLorentzVector particle;
    particle.SetPtEtaPhiM(pt, 0.0, phi, mass);
    return particle;
}

void testGhostRemovalChecksAllLeptonPairs(){
    std::array<h4l::CandidateLepton, 4> leptons = {{
        makeLepton(25.0, 0.0, 0.0, 1),
        makeLepton(20.0, 0.0, 0.01, -1),
        makeLepton(18.0, 1.0, 1.0, 1),
        makeLepton(15.0, -1.0, -1.0, -1),
    }};
    assert(!h4l::passesGhostRemoval(leptons));

    leptons[1] = makeLepton(20.0, 0.0, 0.03, -1);
    assert(h4l::passesGhostRemoval(leptons));
}

void testOppositeSignMassUsesTheMatchingFourVectors(){
    std::array<h4l::CandidateLepton, 4> leptons = {{
        makeLepton(25.0, 0.0, 0.0, 1),
        makeLepton(20.0, 1.0, 2.0, -1),
        makeLepton(18.0, -1.0, -2.0, 1),
        makeLepton(15.0, 0.0, 0.10, -1),
    }};
    assert(!h4l::passesOppositeSignPairMass(leptons));

    leptons[3] = makeLepton(15.0, 0.0, 1.0, -1);
    assert(h4l::passesOppositeSignPairMass(leptons));
}

void testZ1IsClosestToTheNominalMass(){
    const TLorentzVector first = makeParticle(0.0, 0.0, 55.0);
    const TLorentzVector second = makeParticle(0.0, 0.0, 85.0);
    const std::pair<unsigned int, unsigned int> ordered =
        h4l::orderZCandidates(3, 7, first, second, 91.1876);
    assert(ordered.first == 7);
    assert(ordered.second == 3);
}

void testCandidateRankingUsesZ2PtForEqualZ1Distance(){
    assert(h4l::isBetterZZCandidate(2.0, 35.0, 3.0, 100.0));
    assert(h4l::isBetterZZCandidate(2.0, 45.0, 2.0, 40.0));
    assert(!h4l::isBetterZZCandidate(2.0, 35.0, 2.0, 40.0));
}

void testZ1MassUsesTheConfiguredThreshold(){
    const TLorentzVector boundary = makeParticle(0.0, 0.0, 40.0);
    const TLorentzVector passing = makeParticle(0.0, 0.0, 40.1);
    assert(!h4l::passesZ1Mass(boundary, 40.0));
    assert(h4l::passesZ1Mass(passing, 40.0));
}

void testM4lCutUsesTheFourBodyInvariantMass(){
    const TLorentzVector z1 = makeParticle(40.0, 0.0, 20.0);
    const TLorentzVector z2 = makeParticle(40.0, std::acos(-1.0), 20.0);
    assert(z1.M() + z2.M() < 70.0);
    assert((z1 + z2).M() > 70.0);
    assert(h4l::passesFourLeptonMass(z1, z2, 70.0));

    const TLorentzVector lowMassZ1 = makeParticle(0.0, 0.0, 20.0);
    const TLorentzVector lowMassZ2 = makeParticle(0.0, 0.0, 20.0);
    assert(!h4l::passesFourLeptonMass(lowMassZ1, lowMassZ2, 70.0));

    const TLorentzVector boundaryZ1 = makeParticle(0.0, 0.0, 35.0);
    const TLorentzVector boundaryZ2 = makeParticle(0.0, 0.0, 35.0);
    assert(!h4l::passesFourLeptonMass(boundaryZ1, boundaryZ2, 70.0));
}

}  // namespace

int main(){
    testGhostRemovalChecksAllLeptonPairs();
    testOppositeSignMassUsesTheMatchingFourVectors();
    testZ1IsClosestToTheNominalMass();
    testCandidateRankingUsesZ2PtForEqualZ1Distance();
    testZ1MassUsesTheConfiguredThreshold();
    testM4lCutUsesTheFourBodyInvariantMass();
    std::cout << "H4L candidate selection tests passed" << std::endl;
    return 0;
}
