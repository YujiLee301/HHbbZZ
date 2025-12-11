import ROOT


def ZZSelection(Electrons, Muons, Eid, muid):
    Zmass = 91.1876
    foundZZCandidate = False
    z1 = ROOT.TLorentzVector()
    z2 = ROOT.TLorentzVector()
    Elelist = []
    Elechg = []
    Mulist = []
    Muchg = []
    Muiso = []

    for x in Electrons:
        if(x.pdgId>0):
            Elechg.append(-1)
        else:
            Elechg.append(1)

        Ele = ROOT.TLorentzVector()
        Ele.SetPtEtaPhiM(x.pt,x.eta,x.phi,x.mass)
        Elelist.append(Ele)

    for x in Muons:
        if(x.pdgId>0):
            Muchg.append(-1)
        else:
            Muchg.append(1)

        Mu = ROOT.TLorentzVector()
        Mu.SetPtEtaPhiM(x.pt,x.eta,x.phi,x.mass)
        Mulist.append(Mu)
        Muiso.append(x.pfRelIso03_all)

    nTightEle = 0
    nTightMu = 0
    nTightEleChgSum = 0
    nTightMuChgSum = 0
    TightEleindex = []
    TightMuindex = []
    for a in range(len(Eid)):
        if(Eid[a]==1):
            nTightEle += 1
            TightEleindex.append(a)
            nTightEleChgSum += Elechg[a]
    
    for a in range(len(muid)):
        if((muid[a]==1)&Muiso[a]<0.35):
            nTightMu += 1
            TightMuindex.append(a)
            nTightMuChgSum += Muchg[a]
    
    if((nTightMu+nTightEle)<4):
        return foundZZCandidate, z1, z2
    
    if((abs(nTightEleChgSum)+abs(nTightMuChgSum))>(nTightMu+nTightEle-4)):
        return foundZZCandidate, z1, z2

    #Find Z candidates
    Zlist = []
    Zlep1index = []
    Zlep2index = []
    Zlep1pt = []
    Zlep2pt = []
    Zlep1eta = []
    Zlep2eta = []
    Zlep1phi = []
    Zlep2phi = []
    Zlep1mass = []
    Zlep2mass = []
    Zlep1chg = []
    Zlep2chg = []
    Zflavor = []


    if(len(TightEleindex)>1):
        for k in range(0,len(TightEleindex)):
            for j in range(1,len(TightEleindex)):
                if ((Elechg[TightEleindex[k]]+Elechg[TightEleindex[j]])==0):
                    Zcan = ROOT.TLorentzVector()
                    Zcan = Elelist[TightEleindex[k]] + Elelist[TightEleindex[j]]
                    if((Zcan.M()>12)&(Zcan.M()<120)):
                        Zlist.append(Zcan)
                        Zlep1index.append(TightEleindex[k])
                        Zlep2index.append(TightEleindex[j])
                        Zflavor.append('e')
                        Zlep1pt.append(Elelist[TightEleindex[k]].Pt())
                        Zlep2pt.append(Elelist[TightEleindex[j]].Pt())
                        Zlep1eta.append(Elelist[TightEleindex[k]].Eta())
                        Zlep2eta.append(Elelist[TightEleindex[j]].Eta())
                        Zlep1phi.append(Elelist[TightEleindex[k]].Phi())
                        Zlep2phi.append(Elelist[TightEleindex[j]].Phi())
                        Zlep1mass.append(Elelist[TightEleindex[k]].M())
                        Zlep2mass.append(Elelist[TightEleindex[j]].M())
                        Zlep1chg.append(Elechg[TightEleindex[k]])
                        Zlep2chg.append(Elechg[TightEleindex[j]])

    if(len(TightMuindex)>1):
        for k in range(0,len(TightMuindex)):
            for j in range(1,len(TightMuindex)):
                if ((Muchg[TightMuindex[k]]+Muchg[TightMuindex[j]])==0):
                    Zcan = ROOT.TLorentzVector()
                    Zcan = Mulist[TightMuindex[k]] + Mulist[TightMuindex[j]]
                    if((Zcan.M()>12)&(Zcan.M()<120)):
                        Zlist.append(Zcan)
                        Zlep1index.append(TightMuindex[k])
                        Zlep2index.append(TightMuindex[j])
                        Zflavor.append('mu')
                        Zlep1pt.append(Mulist[TightMuindex[k]].Pt())
                        Zlep2pt.append(Mulist[TightMuindex[j]].Pt())
                        Zlep1eta.append(Mulist[TightMuindex[k]].Eta())
                        Zlep2eta.append(Mulist[TightMuindex[j]].Eta())
                        Zlep1phi.append(Mulist[TightMuindex[k]].Phi())
                        Zlep2phi.append(Mulist[TightMuindex[j]].Phi())
                        Zlep1mass.append(Mulist[TightMuindex[k]].M())
                        Zlep2mass.append(Mulist[TightMuindex[j]].M())
                        Zlep1chg.append(Muchg[TightEleindex[k]])
                        Zlep2chg.append(Muchg[TightEleindex[j]])

    if(len(Zlist)<2):
        return foundZZCandidate, SelectedEleIndex, SelectedMuIndex
    
    #Find ZZ candidates
    Z1CanIndex = []
    Z2CanIndex = []
    for m in range(0,len(Zlist)):
        for n in range(1,len(Zlist)):
            if (Zflavor[m]==Zflavor[n]):
               if ((Zlep1index[m] == Zlep1index[n]) | (Zlep2index[m] == Zlep1index[n])): continue  #non-overlapping
               if ((Zlep1index[m] == Zlep2index[n]) | (Zlep2index[m] == Zlep2index[n])): continue
            if (Zlist[m].DeltaR(Zlist[n])<0.02): continue #ghost removal
            nPassPt20 = (Zlep1pt[m]>20) | (Zlep2pt[m]>20) | (Zlep1pt[n]>20) | (Zlep2pt[n]>20)
            nPassPt10 = 0
            if (Zlep1pt[m]>10): nPassPt10 += 1 
            if (Zlep2pt[m]>10): nPassPt10 += 1 
            if (Zlep1pt[n]>10): nPassPt10 += 1 
            if (Zlep2pt[n]>10): nPassPt10 += 1 
            if (nPassPt10 < 2): continue
            if (nPassPt20 == False): continue #lep Pt requirements

            if ((Zlep1chg[m]+Zlep1chg[n])==0):
                lepA = ROOT.TLorentzVector()
                lepB = ROOT.TLorentzVector()
                lepAB = ROOT.TLorentzVector()
                lepA.SetPtEtaPhiM(Zlep1pt[m],Zlep1eta[m],Zlep1phi[m],Zlep1mass[m])
                lepB.SetPtEtaPhiM(Zlep1pt[n],Zlep1eta[n],Zlep1phi[n],Zlep1mass[n])
                lepAB = lepA + lepB
                if(lepAB.M()<4): continue  #QCD suppression

            if ((Zlep1chg[m]+Zlep2chg[n])==0):
                lepA = ROOT.TLorentzVector()
                lepB = ROOT.TLorentzVector()
                lepAB = ROOT.TLorentzVector()
                lepA.SetPtEtaPhiM(Zlep1pt[m],Zlep1eta[m],Zlep1phi[m],Zlep1mass[m])
                lepB.SetPtEtaPhiM(Zlep2pt[n],Zlep2eta[n],Zlep2phi[n],Zlep2mass[n])
                lepAB = lepA + lepB
                if(lepAB.M()<4): continue  #QCD suppression

            if ((Zlep2chg[m]+Zlep1chg[n])==0):
                lepA = ROOT.TLorentzVector()
                lepB = ROOT.TLorentzVector()
                lepAB = ROOT.TLorentzVector()
                lepA.SetPtEtaPhiM(Zlep2pt[m],Zlep2eta[m],Zlep2phi[m],Zlep2mass[m])
                lepB.SetPtEtaPhiM(Zlep1pt[n],Zlep1eta[n],Zlep1phi[n],Zlep1mass[n])
                lepAB = lepA + lepB
                if(lepAB.M()<4): continue  #QCD suppression

            if ((Zlep1chg[m]+Zlep1chg[n])==0):
                lepA = ROOT.TLorentzVector()
                lepB = ROOT.TLorentzVector()
                lepAB = ROOT.TLorentzVector()
                lepA.SetPtEtaPhiM(Zlep2pt[m],Zlep2eta[m],Zlep2phi[m],Zlep2mass[m])
                lepB.SetPtEtaPhiM(Zlep1pt[n],Zlep1eta[n],Zlep1phi[n],Zlep1mass[n])
                lepAB = lepA + lepB
                if(lepAB.M()<4): continue  #QCD suppression

            if (Zlist[m].M()<40 & Zlist[n].M()<40) : continue #Z1 mass

            Z1 = ROOT.TLorentzVector()
            Z2 = ROOT.TLorentzVector()
            if(Zlist[m].M()>Zlist[n].M()):
                Z1 = Zlist[m]
                Z2 = Zlist[n]
            else:
                Z1 = Zlist[n]
                Z2 = Zlist[m]

            passSmartCut = True #Perform Smart cut
            if (Zflavor[m]==Zflavor[n]):
                Za = ROOT.TLorentzVector()
                Zb = ROOT.TLorentzVector()
                lepM1 = ROOT.TLorentzVector()
                lepM2 = ROOT.TLorentzVector()
                lepN1 = ROOT.TLorentzVector()
                lepN2 = ROOT.TLorentzVector()
                lepM1.SetPtEtaPhiM(Zlep1pt[m],Zlep1eta[m],Zlep1phi[m],Zlep1mass[m])
                lepM2.SetPtEtaPhiM(Zlep2pt[m],Zlep2eta[m],Zlep2phi[m],Zlep2mass[m])
                lepM1.SetPtEtaPhiM(Zlep1pt[n],Zlep1eta[n],Zlep1phi[n],Zlep1mass[n])
                lepM2.SetPtEtaPhiM(Zlep2pt[n],Zlep2eta[n],Zlep2phi[n],Zlep2mass[n])
                lepM1chg = Zlep1chg[m]
                lepM2chg = Zlep2chg[m]
                lepN1chg = Zlep1chg[n]
                lepN2chg = Zlep2chg[n]
                if(lepM1chg == lepN1chg):
                    Za = lepM1 + lepN2
                    Zb = lepN1 + lepM2
                else:
                    Za = lepM1 + lepN1
                    Zb = lepN2 + lepM2
                if (abs(Za.M()-Zmass)<abs(Zb.M()-Zmass)):
                    if ( abs(Za.M()-Zmass)<abs(Z1.M()-Zmass) & Zb.M()<12 ):
                        passSmartCut=False
                else:
                    if ( abs(Zb.M()-Zmass)<abs(Z1.M()-Zmass) & Za.M()<12 ):
                        passSmartCut=False

            if (passSmartCut==False): continue
            if (Z1.M()+Z2.M()<70): continue
            foundZZCandidate = True
            if(Zlist[m].M()>Zlist[n].M()):
                Z1CanIndex.append(m)
                Z2CanIndex.append(n)
            else:
                Z1CanIndex.append(n)
                Z2CanIndex.append(m)
    if(foundZZCandidate == False): return foundZZCandidate, z1, z2
    Z1index = Z1CanIndex[0]
    Z2index = Z2CanIndex[0]
    if(len(Z1CanIndex)>1):
        for i in range(len(Z1CanIndex)):
            if(abs(Zlist[Z1CanIndex[i]].M()-Zmass)<abs(Zlist[Z1index].M())):
                Z1index = Z1CanIndex[i]
                Z2index = Z2CanIndex[i]
    
    z1 = Zlist[Z1index]
    z2 = Zlist[Z2index]       
        


    return foundZZCandidate, z1, z2