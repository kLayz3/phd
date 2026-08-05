#include "monad/monad.hxx"

#include <csignal>

#include "TFOOTCalCont.h"
#include "TFOOTCalProc.h"
#include "TFOOTMapCont.h"
#include "TFRSCalCont.h"
#include "TFRSCalProc.h"
#include "TFRSMapCont.h"
#include "util/CLI.h"

using namespace std::literals;
using namespace mnd;

int main(int argc, char *argv[]) {
  using namespace indicators;
  signal(SIGINT, sig_callback_handler);
  signal(SIGSEGV, sig_callback_handler);

  CLI::App app{
      "This program will analyse the mapped ROOT file and perform the clustering of the FOOT data + calibrating FRS data.\n\
Always remember: PHYSICS IS FUN <(^.^)>"};

  std::string fileName, outFile;
  std::string setupFile = PROG_PATH "/params/frs_setup.json";
  std::string footSetupFile = PROG_PATH "/params/foot_setup.json";
  u64 maxEvents = -1;

  add_logged_option<DisplayDefault::No>(app, "-f,--file", fileName,
                                        "Input ROOT file")
      ->required()
      ->expected(1)
      ->check(CLI::ExistingFile);

  add_logged_option<DisplayDefault::No>(
      app, "-o,--output", outFile,
      "Specify output file name. Default same as the input file with \'_cal\' "
      "suffix.")
      ->expected(0, 1);

  add_logged_option<DisplayDefault::No>(
      app, "-m,--max-events", maxEvents,
      "Specify total number of events. Default: all events in the input ROOT "
      "file.")
      ->check(CLI::PositiveNumber);

  add_logged_option(app, "-s,--setup", setupFile,
                    "Specify FRS JSON setup file name.")
      ->check(CLI::ExistingFile);

  add_logged_option(app, "-p,--foot-setup", footSetupFile,
                    "Specify FOOT JSON setup file name.")
      ->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);
  if (outFile.empty())
    outFile = fileName.substr(0, fileName.find('.')) + "_cal.root";

  srand(time(NULL));
  std::vector<TimePoint> tv;

  /* Set up the containers. */
  TFOOTMapCont mfoot[N_FOOT]{}; // input map container.
  for (int i = 0; i < N_FOOT; ++i) {
    mfoot[i].Init({{"FOOT_ID"s, std::to_string(::static_detectors[i])}});
    mfoot[i].Setup();
  }

  TFRSMapCont mfrs{};
  mfrs.Setup();

  TFRSCalCont cfrs{};
  cfrs.Init({{"Setup", setupFile}});
  cfrs.Setup();

  TFOOTCalCont cfoot[N_FOOT]; // output container.
  cfoot[0].SetRegisterBox(true);

  for (int i = 0; i < N_FOOT; ++i) {
    cfoot[i].Init(
        {{"ID"s, std::to_string(mfoot[i].FOOT_N)}, {"Setup"s, footSetupFile}});
    cfoot[i].Setup();
  }

  /* Set up the process pool. */
  auto pool = TAnalysisProcess<>(fileName, outFile, "h103")
                  .emplace_process<TFOOTCalProc>(cfoot[0], mfoot[0])
                  .emplace_process<TFOOTCalProc>(cfoot[1], mfoot[1])
                  .emplace_process<TFOOTCalProc>(cfoot[2], mfoot[2])
                  .emplace_process<TFOOTCalProc>(cfoot[3], mfoot[3])
                  .emplace_process<TFOOTCalProc>(cfoot[4], mfoot[4])
                  .emplace_process<TFOOTCalProc>(cfoot[5], mfoot[5])
                  .emplace_process<TFOOTCalProc>(cfoot[6], mfoot[6])
                  .emplace_process<TFOOTCalProc>(cfoot[7], mfoot[7])
                  .emplace_process<TFRSCalProc>(cfrs, mfrs)
#ifdef MND_DEBUG_ENABLED
                  .MakePool<1>(512);
#else
                  .MakePool<8>(4092);
#endif

  ProgressBar bar{option::BarWidth{50},
                  option::Start{"["},
                  option::Fill{"="},
                  option::Lead{">"},
                  option::Remainder{" "},
                  option::End{"]"},
                  option::PostfixText{"Clustering & Calibration (per event)"},
                  option::ForegroundColor{Color::yellow},
                  option::ShowPercentage{true},
                  option::ShowElapsedTime{true},
                  option::ShowRemainingTime{true},
                  option::FontStyles{std::vector{FontStyle::bold}}};

  tv.emplace_back(TimePoint("start"));

  pool.Start(bar, maxEvents);
  pool.Collect();

  tv.emplace_back(TimePoint("end"));

  PrintElapsed<kSECOND>(std::move(tv));
}
