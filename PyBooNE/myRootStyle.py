#!/usr/local/bin/python

import ROOT as rt
__author__ = 'miceli'

def myRootStyle():
    mysty = rt.TStyle("mysty","My Style")

    # Turn off some borders
    mysty.SetCanvasBorderMode(0);
    mysty.SetFrameBorderMode(0);
    mysty.SetPadBorderMode(0);
    mysty.SetDrawBorder(0);
    mysty.SetCanvasBorderSize(0);
    mysty.SetFrameBorderSize(0);
    mysty.SetPadBorderSize(1);
    mysty.SetTitleBorderSize(0);

    # Say it in black and white!
    mysty.SetAxisColor(1, "xyz");
    mysty.SetCanvasColor(0);
    mysty.SetFrameFillColor(0);
    mysty.SetFrameLineColor(1);
    mysty.SetHistFillColor(0);
    mysty.SetHistLineColor(1);
    #mysty.SetPadColor(1);
    mysty.SetPadColor(rt.kWhite);
    mysty.SetStatColor(0);
    mysty.SetStatTextColor(1);
    mysty.SetTitleColor(1);
    mysty.SetTitleTextColor(1);
    mysty.SetLabelColor(1,"xyz");
    # Show functions in red...
    mysty.SetFuncColor(2);

    # Set the size of the default canvas
    # 600x500 looks almost square
    mysty.SetCanvasDefH(500);
    mysty.SetCanvasDefW(600);
    mysty.SetCanvasDefX(10);
    mysty.SetCanvasDefY(10);

    # Fonts:  I use Helvetica, upright, normal
    #         I sort of wish they had something like "HIGZ portable" of PAW
    style_label_font=42;
    mysty.SetLabelFont(style_label_font,"xyz");
    mysty.SetLabelSize(0.045,"xyz");
    mysty.SetLabelOffset(0.015,"xyz");
    mysty.SetStatFont(style_label_font);
    mysty.SetTitleFont(style_label_font,"xyz"); # axis titles
    mysty.SetTitleFont(style_label_font,"h"); # histogram title
    mysty.SetTitleSize(0.08,"xyz");#mysty.SetTitleSize(0.05,"xyz"); # axis titles
    mysty.SetTitleSize(0.09,"h")#mysty.SetTitleSize(0.05,"h"); # histogram title

    mysty.SetTitleOffset(1.1,"x");
    mysty.SetTitleOffset(1.2,"y");
    mysty.SetStripDecimals(rt.kFALSE); # if we have 1.5 do not set 1.0 -> 1
    mysty.SetTitleX(0.12); # spot where histogram title goes
    mysty.SetTitleY(0.95);
    mysty.SetTitleW(0.78); # width computed so that title is centered
    rt.TGaxis.SetMaxDigits(2); # restrict the number of digits in labels

    # Set Line Widths
    mysty.SetFrameLineWidth(1);
    mysty.SetFuncWidth(2);
    mysty.SetHistLineWidth(2);

    # Set all fill styles to be empty and line styles to be solid
    mysty.SetFrameFillStyle(0);
    mysty.SetHistFillStyle(1001);
    mysty.SetFrameLineStyle(0);
    mysty.SetHistLineStyle(0);
    mysty.SetTitleStyle(0);
    mysty.SetFuncStyle(1);

    # Set margins -- I like to shift the plot a little up and to the
    # right to make more room for axis labels
    mysty.SetPadTopMargin(0.2);
    mysty.SetPadBottomMargin(0.2);
    mysty.SetPadLeftMargin(0.2);
    mysty.SetPadRightMargin(0.1);

    # Set Data/Stat/... and other options
    mysty.SetOptDate(0);
    mysty.SetDateX(0.01);
    mysty.SetDateY(0.01);
    mysty.SetOptFile(0);
    mysty.SetOptFit(111);
    mysty.SetOptLogx(0);
    mysty.SetOptLogy(0);
    mysty.SetOptLogz(0);
    mysty.SetOptStat(1110);# no histogram title
    mysty.SetStatFormat("6.4f");
    mysty.SetFitFormat("6.4f");
    mysty.SetStatStyle(0); # hollow
    #mysty.SetStatStyle(1001); // filled
    mysty.SetStatBorderSize(0);
    mysty.SetStatW(0.25);
    mysty.SetStatH(0.125);
    #mysty.SetStatX(0.90);
    #mysty.SetStatY(0.90);
    mysty.SetStatX(1.0-mysty.GetPadRightMargin()-0.02);
    mysty.SetStatY(1.0-mysty.GetPadTopMargin()-0.02);
    mysty.SetOptTitle(1);
    # Set tick marks and turn off grids
    #mysty.SetNdivisions(1005,"xyz");
    mysty.SetNdivisions(510,"xyz");
    mysty.SetPadTickX(1);
    mysty.SetPadTickY(1);
    mysty.SetTickLength(0.02,"xyz");
    mysty.SetPadGridX(0);
    mysty.SetPadGridY(0);

    # no supressed zeroes!
    mysty.SetHistMinimumZero(rt.kTRUE);


    # Set paper size for life in the US
    #mysty.SetPaperSize(rt.TStyle.kUSLetter);#20,24
    #mysty.SetPaperSize(2000,2400)
    # or europe
    #mysty.SetPaperSize(TStyle::kA4);

    # use a pretty palette for color plots
    # mysty.SetPalette(1,0);

    #mysty.SetLabelColor(1,"xyz");
    return mysty