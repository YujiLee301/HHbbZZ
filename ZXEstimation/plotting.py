"""PyROOT CR histograms and signed MC stacks with data/MC ratio panels."""
import uuid

from histograms import h1, clone, fill
from root_io import EventReader, write_object


def poisson_errors(ROOT, count):
    alpha = 1. - 0.682689492137
    lower = ROOT.Math.gamma_quantile(alpha / 2., count, 1.) if count > 0 else 0.
    upper = ROOT.Math.gamma_quantile_c(alpha / 2., count + 1., 1.)
    return count - lower, upper - count


def draw(ROOT, cfg, data, mc, spec, title, luminosity, destination):
    uid = uuid.uuid4().hex
    canvas = ROOT.TCanvas("canvas_" + uid, "", 800, 800)
    top = ROOT.TPad("top_" + uid, "", 0., 0.30, 1., 1.)
    bottom = ROOT.TPad("bottom_" + uid, "", 0., 0., 1., 0.30)
    for pad in (top, bottom):
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.04)
        pad.Draw()
    top.SetBottomMargin(0.025)
    bottom.SetTopMargin(0.025)
    bottom.SetBottomMargin(0.33)
    top.cd()
    summed = clone(data)
    summed.Reset()
    for hist in mc.values():
        summed.Add(hist)
    frame = clone(data)
    frame.Reset()
    frame.SetTitle("")
    frame.GetXaxis().SetLabelSize(0.)
    frame.GetYaxis().SetTitle("Events / bin")
    frame.GetYaxis().SetTitleOffset(1.45)
    positive, negative = [], []
    for i in range(1, frame.GetNbinsX() + 1):
        positive.append(sum(max(0., h.GetBinContent(i)) for h in mc.values()))
        negative.append(sum(min(0., h.GetBinContent(i)) for h in mc.values()))
    upper = max([1.] + positive + [data.GetBinContent(i) + poisson_errors(ROOT, data.GetBinContent(i))[1]
                                  for i in range(1, data.GetNbinsX() + 1)] +
                [summed.GetBinContent(i) + summed.GetBinError(i) for i in range(1, summed.GetNbinsX() + 1)])
    lower = min([0.] + negative + [summed.GetBinContent(i) - summed.GetBinError(i)
                                  for i in range(1, summed.GetNbinsX() + 1)])
    frame.SetMinimum(lower * 1.2)
    frame.SetMaximum(upper * 1.48)
    frame.Draw("AXIS")
    boxes, legend_objects = [], []
    # Separate positive and negative stacks preserve signed NLO contributions.
    # THStack can display negative bins incorrectly; never clip their yields.
    ypos, yneg = [0.] * frame.GetNbinsX(), [0.] * frame.GetNbinsX()
    for process, settings in cfg["processes"].items():
        hist = mc[process]
        color = ROOT.TColor.GetColor(settings["color"])
        hist.SetFillColor(color)
        hist.SetLineColor(ROOT.kBlack)
        legend_objects.append(hist)
        for i in range(1, hist.GetNbinsX() + 1):
            value = hist.GetBinContent(i)
            running = ypos if value >= 0 else yneg
            start, stop = running[i - 1], running[i - 1] + value
            box = ROOT.TBox(hist.GetXaxis().GetBinLowEdge(i), min(start, stop),
                            hist.GetXaxis().GetBinUpEdge(i), max(start, stop))
            box.SetFillColor(color)
            box.SetLineColor(ROOT.kBlack)
            box.SetLineWidth(1)
            box.Draw("SAME")
            boxes.append(box)
            running[i - 1] = stop
    band = clone(summed)
    band.SetFillStyle(3354)
    band.SetFillColor(ROOT.kGray + 2)
    band.SetMarkerSize(0)
    band.SetLineColor(ROOT.kGray + 2)
    band.Draw("E2 SAME")
    points = ROOT.TGraphAsymmErrors()
    for i in range(1, data.GetNbinsX() + 1):
        y = data.GetBinContent(i)
        lo, hi = poisson_errors(ROOT, y)
        points.SetPoint(i - 1, data.GetBinCenter(i), y)
        points.SetPointError(i - 1, 0., 0., lo, hi)
    points.SetMarkerStyle(20)
    points.SetMarkerSize(0.85)
    points.SetLineColor(ROOT.kBlack)
    points.Draw("PZ SAME")
    legend = ROOT.TLegend(0.57, 0.64, 0.94, 0.91)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetNColumns(2)
    legend.AddEntry(points, "Data", "pe")
    for process, settings in cfg["processes"].items():
        legend.AddEntry(mc[process], settings["label"], "f")
    legend.AddEntry(band, "MC stat.", "f")
    legend.Draw()
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.035)
    text.DrawLatex(0.14, 0.95, cfg["plotting"]["label"])
    text.SetTextAlign(31)
    text.DrawLatex(0.96, 0.95, f"{luminosity:g} fb^{{-1}}")
    text.SetTextAlign(11)
    text.SetTextSize(0.033)
    text.DrawLatex(0.17, 0.87, title)
    top.RedrawAxis()
    bottom.cd()
    ratio_frame = clone(frame)
    ratio_frame.SetMinimum(float(cfg["plotting"]["ratio_range"][0]))
    ratio_frame.SetMaximum(float(cfg["plotting"]["ratio_range"][1]))
    ratio_frame.GetYaxis().SetTitle("Data / MC")
    ratio_frame.GetYaxis().SetNdivisions(505)
    ratio_frame.GetYaxis().SetTitleSize(0.10)
    ratio_frame.GetYaxis().SetTitleOffset(0.58)
    ratio_frame.GetYaxis().SetLabelSize(0.09)
    ratio_frame.GetXaxis().SetTitle(spec["label"])
    ratio_frame.GetXaxis().SetTitleSize(0.12)
    ratio_frame.GetXaxis().SetLabelSize(0.10)
    ratio_frame.Draw("AXIS")
    ratio, ratio_band = ROOT.TGraphAsymmErrors(), ROOT.TGraphAsymmErrors()
    omitted = []
    point = 0
    for i in range(1, data.GetNbinsX() + 1):
        denominator = summed.GetBinContent(i)
        if denominator <= 0:
            omitted.append(i)
            continue
        value = data.GetBinContent(i)
        lo, hi = poisson_errors(ROOT, value)
        ratio.SetPoint(point, data.GetBinCenter(i), value / denominator)
        ratio.SetPointError(point, 0., 0., lo / denominator, hi / denominator)
        ratio_band.SetPoint(point, data.GetBinCenter(i), 1.)
        err = summed.GetBinError(i) / denominator
        ratio_band.SetPointError(point, data.GetBinWidth(i)/2., data.GetBinWidth(i)/2., err, err)
        point += 1
    ratio_band.SetFillColor(ROOT.kGray + 2)
    ratio_band.SetFillStyle(3354)
    ratio_band.Draw("2 SAME")
    unity = ROOT.TLine(frame.GetXaxis().GetXmin(), 1., frame.GetXaxis().GetXmax(), 1.)
    unity.SetLineStyle(2)
    unity.Draw("SAME")
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(0.85)
    ratio.Draw("PZ SAME")
    bottom.RedrawAxis()
    destination.parent.mkdir(parents=True, exist_ok=True)
    for extension in cfg["plotting"]["formats"]:
        canvas.SaveAs(str(destination.with_suffix("." + extension)))
    canvas.Close()
    return summed, omitted


def control_plots(ROOT, cfg, output):
    regions = cfg["plotting"]["regions"]
    variables = {name: spec for name, spec in cfg["variables"].items() if set(spec["regions"]) & set(regions)}
    eras = list(cfg["eras"]) + ["combined"]
    histograms = {}
    for era in eras:
        for region in regions:
            for channel in cfg["channels"]:
                for variable, spec in variables.items():
                    if region not in spec["regions"]:
                        continue
                    key = era, region, channel, variable
                    histograms[key] = {p: h1(ROOT, spec["edges"]) for p in ["data"] + list(cfg["processes"])}
    reader = EventReader(ROOT, cfg)
    yields = {}
    for event in reader.events(regions, variables=variables):
        sample = event["sample"]
        process = "data" if sample["kind"] == "data" else sample["process"]
        for era in (sample["era"], "combined"):
            for channel in event["channels"]:
                group = f"{era}/{event['region']}/{channel}/{process}"
                bucket = yields.setdefault(group, {"entries": 0, "sumw": 0., "sumw2": 0.})
                bucket["entries"] += 1
                bucket["sumw"] += event["weight"]
                bucket["sumw2"] += event["weight"]**2
                for variable, value in event["values"].items():
                    key = era, event["region"], channel, variable
                    fill(histograms[key][process], value, event["weight"], variables[variable].get("fold", True))
    root_file = ROOT.TFile.Open(str(output), "RECREATE")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Cannot create {output}")
    audit = {"input": dict(reader.audit), "yields": yields, "undefined_ratio_bins": {}}
    try:
        for key, histograms_by_process in histograms.items():
            era, region, channel, variable = key
            data = histograms_by_process["data"]
            mc = {p: histograms_by_process[p] for p in cfg["processes"]}
            path = "/".join(key)
            lumi = (sum(float(v["lumi_fb"]) for v in cfg["eras"].values()) if era == "combined"
                    else float(cfg["eras"][era]["lumi_fb"]))
            summed, omitted = draw(ROOT, cfg, data, mc, variables[variable],
                                    f"{era}, {region}, {channel}", lumi, cfg["_output"] / "plots" / path)
            for process, hist in histograms_by_process.items():
                write_object(root_file, path, hist, process)
            write_object(root_file, path, summed, "total_mc")
            if omitted:
                audit["undefined_ratio_bins"][path] = omitted
    finally:
        root_file.Close()
    return audit


def rate_plots(ROOT, cfg, rates):
    for (era, flavor), (hist, valid) in rates.maps.items():
        canvas = ROOT.TCanvas("rates_" + uuid.uuid4().hex, "", 800, 650)
        canvas.SetRightMargin(0.16)
        canvas.SetLeftMargin(0.12)
        hist.SetTitle(f"{era}, |PDG ID| = {flavor};Probe p_{{T}} [GeV];Probe |#eta|;Fake rate")
        hist.SetMinimum(0.)
        hist.SetMaximum(1.)
        ROOT.gStyle.SetPaintTextFormat(".3f")
        hist.Draw("COLZ TEXT E")
        boxes, labels = [], []
        for ix in range(1, hist.GetNbinsX() + 1):
            for iy in range(1, hist.GetNbinsY() + 1):
                if valid.GetBinContent(ix, iy) == 1:
                    continue
                box = ROOT.TBox(hist.GetXaxis().GetBinLowEdge(ix), hist.GetYaxis().GetBinLowEdge(iy),
                                hist.GetXaxis().GetBinUpEdge(ix), hist.GetYaxis().GetBinUpEdge(iy))
                box.SetFillColor(ROOT.kGray)
                box.Draw("SAME")
                label = ROOT.TLatex(hist.GetXaxis().GetBinCenter(ix), hist.GetYaxis().GetBinCenter(iy), "invalid")
                label.SetTextAlign(22)
                label.SetTextSize(0.02)
                label.Draw("SAME")
                boxes.append(box)
                labels.append(label)
        destination = cfg["_output"] / "plots" / "fake_rates" / era
        destination.mkdir(parents=True, exist_ok=True)
        for extension in cfg["plotting"]["formats"]:
            canvas.SaveAs(str(destination / f"fr_{flavor}.{extension}"))
        canvas.Close()
