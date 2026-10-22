// SANDBOX ONLY (scratchpad copy of the repo, 2026-10-20): the electron W-skim
// SF replaced by the pO-derived, CHARGE-DEPENDENT scale factors of the Z->ee
// tag-and-probe split by probe charge (correction/tnp_charge_ele.C, its CSV
// copied to skim/sf/tnp_charge_ele.csv): per (charge, eta_SC cell) the
// cut-and-count SF of the reco link (SC-energy probe four-vector) and of the
// WP90 ID+iso link, cells |eta_SC| < 1.6 / fwd eta_SC <= -1.6 / bwd eta_SC >= +1.6.
// The trigger factor stays the pO MB one (one inclusive number, already
// applied by electron_sf.h). Structure as electron_sf.h::W: RECO x [iso ?
// IDISO : 1] x TRIG of the leading electron.
#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace pOChargeSF
{
const double kEtaSplit = 1.6; // = tnp_charge_ele.C's pOTnpQ::kEtaSplit
const char *const kCellName[6] = {"em_lo", "em_fwd", "em_bwd", "ep_lo", "ep_fwd", "ep_bwd"};
inline int CellOf(int q, double etaSC) // = pOTnpQ::CellOf
{
  const int base = (q < 0) ? 0 : 3;
  if (std::fabs(etaSC) < kEtaSplit) return base;
  return base + (etaSC < 0.0 ? 1 : 2);
}

// one CSV line -> fields (double quotes protect commas)
inline std::vector<std::string> SplitCsv(const std::string &line)
{
  std::vector<std::string> out;
  std::string cur;
  bool q = false;
  for (char c : line)
  {
    if (c == '"') q = !q;
    else if (c == ',' && !q) { out.push_back(cur); cur.clear(); }
    else cur += c;
  }
  out.push_back(cur);
  return out;
}

struct Table
{
  double reco[6] = {0, 0, 0, 0, 0, 0}, idiso[6] = {0, 0, 0, 0, 0, 0};
  // diagnostics over the nominal W selection: sum w_gen x {old, new} SF per cell
  double sumOld[6] = {0, 0, 0, 0, 0, 0}, sumNew[6] = {0, 0, 0, 0, 0, 0}, sumW[6] = {0, 0, 0, 0, 0, 0};

  bool Load(const std::string &path, std::ostream &log)
  {
    std::ifstream in(path.c_str());
    if (!in) { log << "[CHARGESF] ERR cannot open " << path << "\n"; return false; }
    std::string line;
    std::getline(in, line);
    const std::vector<std::string> head = SplitCsv(line);
    int iLink = -1, iBin = -1, iSf = -1;
    for (size_t k = 0; k < head.size(); ++k)
    {
      if (head[k] == "link") iLink = (int)k;
      if (head[k] == "bin") iBin = (int)k;
      if (head[k] == "sf_cc") iSf = (int)k;
    }
    if (iLink < 0 || iBin < 0 || iSf < 0) { log << "[CHARGESF] ERR " << path << ": no link/bin/sf_cc column\n"; return false; }
    int found = 0;
    while (std::getline(in, line))
    {
      const std::vector<std::string> f = SplitCsv(line);
      if ((int)f.size() <= iSf) continue;
      for (int c = 0; c < 6; ++c)
        if (f[iBin] == kCellName[c])
        {
          const double v = std::stod(f[iSf]);
          if (f[iLink] == "reco") { reco[c] = v; ++found; }
          if (f[iLink] == "idiso") { idiso[c] = v; ++found; }
        }
    }
    if (found != 12) { log << "[CHARGESF] ERR " << path << ": " << found << " of 12 cells found\n"; return false; }
    log << "[CHARGESF] pO charge-dependent electron SF (cut-and-count, reco SC-energy link) from " << path << "\n";
    for (int c = 0; c < 6; ++c)
      log << Form("[CHARGESF]   %-6s reco %.4f  idiso %.4f  -> reco x idiso %.4f\n", kCellName[c], reco[c], idiso[c],
                  reco[c] * idiso[c]);
    return true;
  }

  // the replacement of electron_sf.h::W's nominal factor; trig = its trigger factor
  double W(int q, double scEta, bool passIso, double trig) const
  {
    const int c = CellOf(q, scEta);
    return reco[c] * (passIso ? idiso[c] : 1.0) * trig;
  }
  void Count(int q, double scEta, double wGen, double sfOld, double sfNew)
  {
    const int c = CellOf(q, scEta);
    sumW[c] += wGen;
    sumOld[c] += wGen * sfOld;
    sumNew[c] += wGen * sfNew;
  }
  void Print(std::ostream &log) const
  {
    log << "[CHARGESF] <SF> over the nominal W selection (iso-pass, pT > 25), EGM -> pO charge-dependent:\n";
    for (int c = 0; c < 6; ++c)
      if (sumW[c] > 0)
        log << Form("[CHARGESF]   %-6s <SF_EGM> %.4f  <SF_pO> %.4f  ratio %.4f\n", kCellName[c], sumOld[c] / sumW[c],
                    sumNew[c] / sumW[c], sumNew[c] / sumOld[c]);
  }
};
} // namespace pOChargeSF
