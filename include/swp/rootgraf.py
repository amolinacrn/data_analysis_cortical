import numpy as np
import ROOT 


def funcion_normal(x,par):
    xx=x[0]
    f = par[2]*par[0]*np.exp(-np.power(xx-par[1],2)/(2*np.power(par[0],2)))
    return f


def graf_root():
    c=ROOT.TCanvas("cV5","migrafico",300,500,550,470)

    data=np.genfromtxt("cfswp.txt")#"MC_pesada.txt")

    contador = len(data)
    #print(data[0])

    H_x = np.array(data,dtype = float)
    H_y = np.ones(contador)



    h1 = ROOT.TH1F("h1","h1 title",80,0.55,0.72);
    h1.FillN(contador, H_x, H_y)

    h = h1.Clone("h");

    h.Scale(1./h.Integral())
    #h.GetBinContent(i)
    #h.GetBinWidth(i)
    #h.SetFillColor(47) 
    h.SetMarkerStyle(33)
    #h.SetMarkerColorAlpha(1)
    h.SetLineColorAlpha(1,1)
    h.SetMarkerColor(1) 

    f3 = ROOT.TF1("fnormal",funcion_normal,0.55,0.72,3)
    f3.SetNpx(1000)
    f3.SetLineWidth(2)
    f3.SetLineColor(2)
    f3.SetParameters(0.66,0.24,0.2)
    #h.Draw("HIST")
    h.Fit(f3,"","",0,1)

    c.Update()

    input("Prompt: ")