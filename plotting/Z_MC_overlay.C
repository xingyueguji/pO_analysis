void Z_MC_overlay()
{
    // 8 versions of Signal Z MC
    TFile *f1[8];

    for (int i = 0; i < 8; i++)
    {
        f1[i] = new TFile(Form("../skim/rootfile/WToMuNu_pO_PFMet_MC_%d_hist.root", i + 1), "READ");
        if (f1[i] == nullptr)
        {
            std::cout << "Cannot found root file!" << std::endl;
            break;
        }
    }

    TH1D *h_MET[8];
    TH1D *h_MT[8];

    for (int i = 0; i < 8; i++)
    {
        h_MET[i] = (TH1D *)f1[i]->Get("hMet_Final");
        h_MT[i] = (TH1D *)f1[i]->Get("hMt_Final");
    }


}