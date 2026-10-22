// correction/abcd_ygroups.h
//
// The RAPIDITY GROUPS of the rapidity-resolved in-fit ABCD (QCD_MODE=abcdy,
// 2026-10-19, user decision) -- SINGLE SOURCE, included by
// correction/qcd_abcd.C (the per-group ABCD counts abcd_counts_*_<bin>_<g>),
// plotting/mtandmet.C (the qcd_abcdy SR template, the per-group CR dirs and
// the group sidecar combine_input_W_leppt_mt40_abcdy.txt that the Combine
// fork reads) and the diagnostics (correction/qcd_ysplit_pergroup.C).
//
// Three groups of the 12 rapidity bins, y = -eta_lab (proton-going > 0, the
// skim's convention), the SAME bin INDICES for the lab and the FB binning:
//   fwd = y10, y11        (group 1, most forward)
//   bwd = y0,  y1         (group 2, most backward)
//   cen = y2 .. y9        (group 3)
// Each group gets its own transfer factor from its own control regions,
// A_g = B_g * C40_g / D_g; the split INSIDE a group, the lepton-pT shape and
// the electron fake-factor weighting stay inclusive (the qcd_abcd template).
#ifndef PO_ABCD_YGROUPS_H
#define PO_ABCD_YGROUPS_H

#include <string>
#include <vector>

namespace pOAbcdY
{
inline constexpr int kNGroups = 3;
inline constexpr int kNYBins = 12; // = pOSkim::kNY (checked by the users)

// group names, in the user's order (group 1, 2, 3)
inline const char *const kGroupName[kNGroups] = {"fwd", "bwd", "cen"};

// group index (0 = fwd, 1 = bwd, 2 = cen) of rapidity bin iy (lab and fb alike)
inline int GroupOf(int iy)
{
  if (iy >= 10) return 0; // fwd: y10, y11
  if (iy <= 1) return 1;  // bwd: y0, y1
  return 2;               // cen: y2 .. y9
}

// the bin indices of group g, ascending
inline std::vector<int> BinsOf(int g)
{
  std::vector<int> v;
  for (int iy = 0; iy < kNYBins; ++iy)
    if (GroupOf(iy) == g) v.push_back(iy);
  return v;
}
} // namespace pOAbcdY

#endif
