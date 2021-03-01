#ifndef PLOT_H
#define PLOT_H

#include <TCanvas.h>
#include <TColor.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <fmt/core.h>
#include <vector>

namespace plot {

  const int kMpBlue   = TColor::GetColor(0x1f, 0x77, 0xb4);
  const int kMpOrange = TColor::GetColor(0xff, 0x7f, 0x0e);
  const int kMpGreen  = TColor::GetColor(0x2c, 0xa0, 0x2c);
  const int kMpRed    = TColor::GetColor(0xd6, 0x27, 0x28);
  const int kMpPurple = TColor::GetColor(0x94, 0x67, 0xbd);
  const int kMpBrown  = TColor::GetColor(0x8c, 0x56, 0x4b);
  const int kMpPink   = TColor::GetColor(0xe3, 0x77, 0xc2);
  const int kMpGrey   = TColor::GetColor(0x7f, 0x7f, 0x7f);
  const int kMpMoss   = TColor::GetColor(0xbc, 0xbd, 0x22);
  const int kMpCyan   = TColor::GetColor(0x17, 0xbe, 0xcf);

  const std::vector<int> kPalette = {kMpBlue,  kMpOrange, kMpGreen, kMpRed,  kMpPurple,
                                     kMpBrown, kMpPink,   kMpGrey,  kMpMoss, kMpCyan};

  void draw_label(int ebeam, int pbeam, const std::string_view detector)
  {
    auto t = new TPaveText(.15, 0.800, .7, .925, "NB NDC");
    t->SetFillColorAlpha(kWhite, 0.4);
    t->SetTextFont(43);
    t->SetTextSize(25);
    t->AddText(fmt::format("#bf{{{} }}SIMULATION", detector).c_str());
    t->AddText(fmt::format("{} GeV on {} GeV", ebeam, pbeam).c_str());
    t->SetTextAlign(12);
    t->Draw();
  }

} // namespace plot

#endif
