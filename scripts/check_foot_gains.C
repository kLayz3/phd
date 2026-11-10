// Run from the repository root:
//   root -l -b -q 'scripts/check_foot_gains.C()'
// Optional arguments: calibration directory, setup JSON, verbose differences.
//   root -l -b -q 'scripts/check_foot_gains.C("caldir", "params/foot_setup.json", true)'
// Returns 0 if every file matches all eight detectors, 1 for mismatches/missing
// data, and 2 for configuration errors. Opens ROOT files read-only.

R__ADD_INCLUDE_PATH(includes)
R__LOAD_LIBRARY(libStructures.so)

#include "../includes/util/JSONParser.h"
#include <TFile.h>
#include <TKey.h>
#include <TStreamerInfo.h>
#include <algorithm>
#include <array>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

int check_foot_gains(const char* caldir = "caldir",
                     const char* setup_file = "params/foot_setup.json",
                     bool verbose = false) {
    using json = nlohmann::json;
    namespace fs = std::filesystem;
    try {
        const auto setup = ParseJSON(setup_file);
        std::array<std::string, 8> names{};
        std::array<json, 8> expected;
        for (auto it = setup.begin(); it != setup.end(); ++it) {
            if (!it.value().is_object() || !it.value().contains("gain")) continue;
            const int n = it.value().at("N").get<int>();
            if (n < 0 || n >= 8 || !names[n].empty())
                throw std::runtime_error("Invalid or duplicate detector N in " + it.key());
            names[n] = it.key();
            expected[n] = it.value().at("gain");
        }
        for (int n = 0; n < 8; ++n)
            if (names[n].empty())
                throw std::runtime_error("Missing gain configuration for detector N=" + std::to_string(n));

        std::vector<fs::path> files;
        for (const auto& entry : fs::recursive_directory_iterator(caldir))
            if (entry.is_regular_file() && entry.path().extension() == ".root")
                files.push_back(entry.path());
        std::sort(files.begin(), files.end());
        if (files.empty()) throw std::runtime_error("No ROOT files found in " + std::string(caldir));

        unsigned matched_files = 0, matched = 0, different = 0, missing = 0, errors = 0;
        for (const auto& path : files) {
            std::cout << path.string() << '\n';
            std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
            if (!file || file->IsZombie()) {
                std::cout << "  ERROR: cannot open file\n";
                ++errors;
                continue;
            }

            // Inspect the on-file schema before deserializing. Older FOOTParam
            // layouts share a class version with newer ones but have no gain;
            // reading those with the current dictionary can misinterpret data.
            bool found_schema = false, has_gain = true;
            for (auto* obj : *file->GetStreamerInfoList()) {
                auto* info = dynamic_cast<TStreamerInfo*>(obj);
                if (!info || std::string(info->GetName()) != "FOOTParam") continue;
                found_schema = true;
                if (!info->GetElements()->FindObject("gain")) has_gain = false;
            }

            unsigned file_matches = 0;
            for (int n = 0; n < 8; ++n) {
                const std::string key = "FOOT" + std::to_string(n) + "_setup";
                const std::string label = "  " + names[n] + " (" + key + "): ";
                auto* stored_key = file->GetKey(key.c_str());
                if (!stored_key) {
                    std::cout << label << "MISSING object\n";
                    ++missing;
                    continue;
                }
                if (std::string(stored_key->GetClassName()) != "FOOTParam" || !found_schema) {
                    std::cout << label << "ERROR: unexpected class or missing schema\n";
                    ++errors;
                    continue;
                }
                if (!has_gain) {
                    std::cout << label << "NO STORED GAIN FIELD (older schema)\n";
                    ++missing;
                    continue;
                }
                auto* param = file->Get<FOOTParam>(key.c_str());
                if (!param || param->N != n) {
                    std::cout << label << "ERROR: unreadable object or wrong detector N\n";
                    ++errors;
                    continue;
                }

                // Copy doubles directly: TBufferJSON's default decimal precision
                // would round coefficients and spoil an exact comparison.
                json actual = {{"fit", json::array()}, {"nominal_value", json::array()}};
                for (const auto& asic : param->gain.fit) {
                    json polys = json::array();
                    for (const auto& poly : asic.multi_poly)
                        polys.push_back({{"Z", poly.Z}, {"pol", poly.pol}});
                    actual["fit"].push_back({{"multi_poly", polys}});
                }
                for (const auto& nominal : param->gain.nominal_value)
                    actual["nominal_value"].push_back({{"Z", nominal.Z}, {"val", nominal.val}});

                const auto delta = json::diff(actual, expected[n]);
                if (delta.empty()) {
                    ++matched;
                    ++file_matches;
                    if (verbose) std::cout << label << "MATCH\n";
                } else {
                    ++different;
                    std::cout << label << "DIFFERENT (" << delta.size() << " JSON changes)\n";
                    if (verbose) std::cout << "    ROOT -> setup JSON: " << delta.dump() << '\n';
                    else std::cout << "    First change: " << delta.front().dump() << '\n';
                }
            }
            if (file_matches == 8) {
                ++matched_files;
                std::cout << "  MATCH: all 8 detectors\n";
            } else {
                std::cout << "  " << file_matches << "/8 detectors match\n";
            }
        }
        std::cout << "\nFiles: " << files.size() << ", all-eight matches: " << matched_files
                  << "\nDetectors: " << matched << " match, " << different << " different, "
                  << missing << " missing object/gain, " << errors << " errors\n";
        return different || missing || errors ? 1 : 0;
    } catch (const std::exception& e) {
        std::cerr << "check_foot_gains: " << e.what() << '\n';
        return 2;
    }
}
