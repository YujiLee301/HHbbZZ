from array import array
import math
import uuid


def h1(ROOT, edges, title=""):
    hist = ROOT.TH1D("h_" + uuid.uuid4().hex, title, len(edges) - 1, array("d", edges))
    hist.SetDirectory(0)
    hist.Sumw2()
    return hist


def h2(ROOT, xedges, yedges, title=""):
    hist = ROOT.TH2D("h_" + uuid.uuid4().hex, title,
                     len(xedges) - 1, array("d", xedges), len(yedges) - 1, array("d", yedges))
    hist.SetDirectory(0)
    hist.Sumw2()
    return hist


def clone(hist):
    copied = hist.Clone("h_" + uuid.uuid4().hex)
    copied.SetDirectory(0)
    return copied


def fill(hist, value, weight, fold=True):
    if fold:
        axis = hist.GetXaxis()
        if value < axis.GetXmin():
            value = axis.GetBinCenter(1)
        elif value >= axis.GetXmax():
            value = axis.GetBinCenter(hist.GetNbinsX())
    hist.Fill(value, weight)


def total(hist):
    return {"yield": sum(hist.GetBinContent(i) for i in range(hist.GetNbinsX() + 2)),
            "stat_error": math.sqrt(sum(hist.GetBinError(i) ** 2 for i in range(hist.GetNbinsX() + 2)))}
