#include <TFile.h>
#include <TGraphErrors.h>
#include <TString.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <sstream>

// --------------------------------------------------
void plotRpOtheory(const char *fileList = "filelist_theory.txt",
                   const char *outRoot = "./RpO_rootfile/RpO_FB_graphs.root")
{
    TFile *fout = TFile::Open(outRoot, "RECREATE");
    if (!fout || fout->IsZombie())
    {
        std::cerr << "ERROR: cannot create output ROOT file " << outRoot << "\n";
        return;
    }

    std::ifstream finList(fileList);
    if (!finList.is_open())
    {
        std::cerr << "ERROR: cannot open file list " << fileList << "\n";
        return;
    }

    std::string txtFile;
    int nGraphs = 0;

    while (std::getline(finList, txtFile))
    {
        if (txtFile.empty())
            continue;

        std::ifstream fin(txtFile);
        if (!fin.is_open())
        {
            std::cerr << "WARN: cannot open " << txtFile << ", skipping\n";
            continue;
        }

        // --- store input as y -> (value, error)
        std::map<double, std::pair<double, double>> data;

        std::string line;
        while (std::getline(fin, line))
        {
            if (line.empty() || line[0] == '#')
                continue;

            double y, ratio, err;
            std::istringstream ss(line);
            if (!(ss >> y >> ratio >> err))
                continue;

            data[y] = {ratio, err};
        }
        fin.close();

        if (data.empty())
        {
            std::cerr << "WARN: no valid data in " << txtFile << "\n";
            continue;
        }

        // --- build FB ratio
        std::vector<double> vx, vy, vex, vey;

        for (const auto &kv : data)
        {
            double y = kv.first;
            if (y <= 0)
                continue; // only forward

            auto itB = data.find(-y);
            if (itB == data.end())
                continue; // no backward partner

            double Rf = kv.second.first;
            double ef = kv.second.second;
            double Rb = itB->second.first;
            double eb = itB->second.second;

            if (Rb <= 0 || Rf <= 0)
                continue;

            double Rfb = Rf / Rb;
            double eRfb = Rfb * std::sqrt(
                                    (ef / Rf) * (ef / Rf) +
                                    (eb / Rb) * (eb / Rb));

            vx.push_back(y);
            vy.push_back(Rfb);
            vex.push_back(0.0);
            vey.push_back(eRfb);
        }

        if (vx.empty())
        {
            std::cerr << "WARN: no FB points for " << txtFile << "\n";
            continue;
        }

        // --- graph name
        TString gname = txtFile;
        if (gname.BeginsWith("./"))
            gname.Remove(0, 2);
        gname.ReplaceAll("/", "_");
        gname.ReplaceAll(".txt", "");
        gname += "_FB";

        TGraphErrors *gr = new TGraphErrors(
            vx.size(),
            vx.data(),
            vy.data(),
            vex.data(),
            vey.data());

        gr->SetName(gname);
        gr->SetTitle(gname);

        fout->cd();
        gr->Write();
        ++nGraphs;
    }

    finList.close();
    fout->Close();

    std::cout << "Done.\n"
              << "  Wrote " << nGraphs << " FB TGraphErrors to " << outRoot << "\n";
}