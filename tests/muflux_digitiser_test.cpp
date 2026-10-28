// SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
//
// SPDX-License-Identifier: LGPL-3.0-or-later

// muflux_digitiser_test.cpp — wire table, distance to the wire and the
// digitisation rules of the 2018 muon-flux spectrometer

#include "muflux/muflux_digitiser.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace mf = SHiP::shannon::muflux;

int failures = 0;

void check(bool ok, std::string const& what) {
    if (!ok) {
        std::fprintf(stderr, "FAILED: %s\n", what.c_str());
        ++failures;
    }
}

bool near(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) <= tol;
}

// No smearing, and a fixed draw for the efficiencies.
struct FixedRng {
    double u = 0.5;
    double uniform53(double lo, double hi) { return lo + u * (hi - lo); }
    double gaussian(double mean, double /*sigma*/) { return mean; }
};

std::filesystem::path write_table(std::string const& name, std::string const& text) {
    auto const path = std::filesystem::temp_directory_path() / name;
    std::ofstream{path} << text;
    return path;
}

bool throws(std::string const& text) {
    auto const path = write_table("muflux_bad_table.txt", text);
    try {
        (void)mf::read_wire_table(path.string());
    } catch (std::runtime_error const&) {
        return true;
    }
    return false;
}

mf::HitView tube_hit(std::int32_t id, mf::Vec3 pos, mf::Vec3 mom, double time = 0,
                     std::int32_t track = 1) {
    return {.detector_id = mf::kDriftTubeDetectorId,
            .channel = id,
            .track_id = track,
            .pdg_code = 13,
            .position = pos,
            .momentum = mom,
            .time = time};
}

}  // namespace

int main() {
    // Two lines of the legacy table: an x-view tube (wire along y) and a
    // u-view tube (wire at 60 degrees to the vertical).
    auto const table =
        write_table("muflux_wires.txt",
                    "# detID  x1 y1 z1  x2 y2 z2  r\n"
                    "10002001 -21.77 -56.895 19.32 -21.77 52.905 19.32 1.815\n"
                    "11002001 42.5197 -50.0046 50.5382 -52.7610 4.5631 50.5382 1.815\n");
    auto const wires = mf::read_wire_table(table.string());
    check(wires.size() == 2, "wire table: two wires");
    check(near(wires.at(10002001).end1[0], -217.7) && near(wires.at(10002001).end2[1], 529.05) &&
              near(wires.at(10002001).gas_radius, 18.15),
          "wire table: cm -> mm");

    check(throws(""), "an empty table");
    check(throws("10002001 1 2 3\n"), "a short line");
    check(throws("1 0 0 0 0 1 0 1\n1 0 0 0 0 1 0 1\n"), "a duplicate id");
    check(throws("1 0 0 0 0 0 0 1\n"), "a zero-length wire");

    auto const& x = wires.at(10002001);
    // A track along z, 5 mm beside the wire in x.
    check(near(mf::distance_to_wire({-212.7, 0, 100}, {0, 0, 1}, x), 5.0, 1e-9),
          "distance: track along z");
    // The same track tilted in y (along the wire): still 5 mm.
    check(near(mf::distance_to_wire({-212.7, 0, 100}, {0, 0.3, 1}, x), 5.0, 1e-9),
          "distance: tilt along the wire changes nothing");
    // Tilted in x: the track passes the wire plane z = 193.2 at x = -212.7 + 0.1 * 93.2.
    {
        double const x_at = -212.7 + 0.1 * (193.2 - 100.);
        double const expected = std::abs(x_at - (-217.7)) / std::sqrt(1 + 0.01);
        check(near(mf::distance_to_wire({-212.7, 0, 100}, {0.1, 0, 1}, x), expected, 1e-9),
              "distance: tilt across the wire");
    }
    check(near(mf::distance_to_wire({-212.7, 0, 193.2}, {0, 1, 0}, x), 5.0, 1e-9),
          "distance: track parallel to the wire");

    // The stereo tube: a track along z through a point 3 mm from the wire,
    // measured perpendicular to the wire within the plane z = const.
    {
        auto const& u = wires.at(11002001);
        mf::Vec3 const w = mf::sub(u.end2, u.end1);
        double const wn = mf::norm(w);
        mf::Vec3 const perp{-w[1] / wn, w[0] / wn, 0};  // in-plane normal to the wire
        mf::Vec3 const mid{(u.end1[0] + u.end2[0]) / 2, (u.end1[1] + u.end2[1]) / 2, 0};
        mf::Vec3 const p{mid[0] + 3 * perp[0], mid[1] + 3 * perp[1], 400};
        check(near(mf::distance_to_wire(p, {0, 0, 1}, u), 3.0, 1e-9), "distance: stereo tube");
    }

    // Digitisation rules, first for entry points (sd_mode 'crossing').
    mf::DigitiserConfig config;
    config.step_segments = false;
    FixedRng rng;
    mf::Counters counters;
    std::vector<mf::HitView> hits{
        tube_hit(10002001, {-212.7, 0, 180}, {0, 0, 10}, 1.0, 7),    // 5 mm, early
        tube_hit(10002001, {-219.7, 0, 180}, {0, 0, 10}, 400.0, 8),  // 2 mm, but 400 ns later
        tube_hit(11002001, {0, 0, 500}, {0, 0, 10}, 0.0, 7),
        {.detector_id = mf::kOtherDetectorId,
         .channel = 10050,
         .track_id = 7,
         .pdg_code = 13,
         .position = {0, 0, 8800},
         .momentum = {0, 0, 10},
         .time = 0},
        {.detector_id = mf::kOtherDetectorId,
         .channel = 10050,
         .track_id = 9,
         .pdg_code = 11,
         .position = {0, 0, 8800},
         .momentum = {0, 0, 1},
         .time = 0},
        {.detector_id = mf::kOtherDetectorId,
         .channel = 11100,
         .track_id = 7,
         .pdg_code = 13,
         .position = {0, 0, 8801},
         .momentum = {0, 0, 10},
         .time = 0},
        {.detector_id = mf::kOtherDetectorId,
         .channel = 0,
         .track_id = 7,
         .pdg_code = 13,
         .position = {0, 0, 115},
         .momentum = {0, 0, 10},
         .time = 0},
        {.detector_id = mf::kOtherDetectorId,
         .channel = 0,
         .track_id = 7,
         .pdg_code = 13,
         .position = {0, 0, 1363},
         .momentum = {0, 0, 10},
         .time = 0},
    };
    auto const d = mf::digitise(hits, wires, config, rng, counters);
    check(d.dt_id.size() == 2 && d.dt_id[0] == 10002001 && d.dt_id[1] == 11002001,
          "one digit per tube, ordered by id");
    check(d.dt_id.size() == 2 && d.dt_track[0] == 7 && near(d.dt_drift[0], 5.0, 1e-5),
          "the earliest hit of a tube is kept");
    check(d.dt_id.size() == 2 && near(d.dt_drift_time[0], 5.0 * 72, 1e-3),
          "drift time = distance / v_drift");
    check(counters.dt_hits == 3 && counters.dt_later == 1, "counters: later hit");
    check(
        d.rpc_id.size() == 2 && d.rpc_id[0] == 10050 && d.rpc_track[0] == 7 && d.rpc_id[1] == 11100,
        "a strip fires once, with the track of its first hit");
    check(counters.rpc_hits == 3 && counters.rpc_shared == 1, "counters: shared strip");
    check(d.n_sa == 1 && d.n_sb == 1, "SA and SB told apart by z");

    // The earliest signal (hit time + drift time) wins, not the first hit.
    {
        mf::Counters k;
        auto const d1 = mf::digitise({tube_hit(10002001, {-212.7, 0, 180}, {0, 0, 1}, 1.0, 7),
                                      tube_hit(10002001, {-219.7, 0, 180}, {0, 0, 1}, 100.0, 8)},
                                     wires, config, rng, k);
        check(d1.dt_id.size() == 1 && d1.dt_track[0] == 8 && near(d1.dt_drift[0], 2.0, 1e-5),
              "a later hit with a shorter drift time replaces the first");
    }

    // Charges from PDG codes, and hits of neutral particles dropped.
    check(mf::three_charge(13) == -3 && mf::three_charge(-13) == 3 && mf::three_charge(14) == 0 &&
              mf::three_charge(22) == 0 && mf::three_charge(211) == 3 &&
              mf::three_charge(-321) == -3 && mf::three_charge(130) == 0 &&
              mf::three_charge(310) == 0 && mf::three_charge(2212) == 3 &&
              mf::three_charge(2112) == 0 && mf::three_charge(3122) == 0 &&
              mf::three_charge(-3222) == -3 && mf::three_charge(421) == 0 &&
              mf::three_charge(411) == 3 && mf::three_charge(521) == 3 &&
              mf::three_charge(1000020040) == 6,
          "three_charge");
    {
        mf::Counters k;
        auto nu = tube_hit(10002001, {-212.7, 0, 180}, {0, 0, 10}, 0.0, 9);
        nu.pdg_code = 14;
        auto const dn = mf::digitise({nu, tube_hit(10002001, {-219.7, 0, 180}, {0, 0, 10}, 5.0, 7)},
                                     wires, config, rng, k);
        check(dn.dt_id.size() == 1 && dn.dt_track[0] == 7 && k.dt_neutral == 1 && k.dt_hits == 1,
              "a neutrino does not take the tube of the muon behind it");
    }

    // Dead channel, and a distance beyond the gas radius.
    {
        mf::DigitiserConfig c2;
        c2.step_segments = false;
        c2.dead_channels = {11002001};
        mf::Counters k;
        auto const d2 = mf::digitise({tube_hit(11002001, {0, 0, 500}, {0, 0, 1}),
                                      tube_hit(10002001, {-217.7 + 30, 0, 180}, {0, 0, 1})},
                                     wires, c2, rng, k);
        check(d2.dt_id.size() == 1 && d2.dt_id[0] == 10002001 && k.dt_dead == 1,
              "a dead channel gives no digit");
        check(k.dt_clipped == 1 && d2.dt_id.size() == 1 && near(d2.dt_true_dist[0], 18.15, 1e-4),
              "a distance beyond the gas radius is set to it");
    }

    // Efficiencies: the fixed draw of 0.5 passes 0.6 and fails 0.4.
    {
        mf::DigitiserConfig c3;
        c3.step_segments = false;
        mf::Counters k;
        c3.dt_efficiency = 0.6;
        c3.rpc_efficiency = 0.4;
        auto const d3 = mf::digitise(hits, wires, c3, rng, k);
        check(d3.dt_id.size() == 2 && d3.rpc_id.empty() && k.rpc_inefficient == 3 &&
                  k.dt_inefficient == 0,
              "efficiencies");
    }

    // Steps (sd_mode 'scoring'): the distance of the step segment, and the
    // closest ionisation of all the hits of a tube.
    {
        auto const& x = wires.at(10002001);  // wire at x = -217.7 mm, z = 193.2 mm, along y
        check(near(mf::distance_segment_to_wire({-212.7, 0, 180}, {0, 0, 1}, 36.3, x), 5.0, 1e-9),
              "segment through the closest point");
        check(near(mf::distance_segment_to_wire({-212.7, 0, 180}, {0, 0, 1}, 3.0, x),
                   std::hypot(5.0, 10.2), 1e-9),
              "segment ending before the closest point");
        check(near(mf::distance_segment_to_wire({-212.7, 0, 200}, {0, 0, 1}, 0.0, x),
                   std::hypot(5.0, 6.8), 1e-9),
              "zero-length step: distance of its point");

        mf::DigitiserConfig cs;  // step_segments = true by default
        mf::Counters k;
        // a muon crossing the tube 5 mm from the wire, in two steps, and a delta
        // electron produced on its path pointing at the wire but stopping after
        // 1 mm: the muon's 5 mm wins
        auto mu1 = tube_hit(10002001, {-212.7, 0, 176}, {0, 0, 10}, 0.0, 7);
        mu1.path_length = 15;
        auto mu2 = tube_hit(10002001, {-212.7, 0, 191}, {0, 0, 10}, 0.05, 7);
        mu2.path_length = 20;
        auto delta = tube_hit(10002001, {-212.7, 0, 185}, {-1, 0, 0}, 0.03, 12);
        delta.pdg_code = 11;
        delta.path_length = 1;
        auto const ds = mf::digitise({mu1, delta, mu2}, wires, cs, rng, k);
        check(ds.dt_id.size() == 1 && ds.dt_track[0] == 7 && near(ds.dt_true_dist[0], 5.0, 1e-5),
              "steps: the muon's closest approach, not the delta's direction");
        // a delta that really ionises closer to the wire wins
        delta.path_length = 4;
        delta.position = {-212.7, 0, 193.2};
        auto const dc = mf::digitise({mu1, delta, mu2}, wires, cs, rng, k);
        check(dc.dt_id.size() == 1 && dc.dt_track[0] == 12 && near(dc.dt_true_dist[0], 1.0, 1e-5),
              "steps: a closer ionisation wins");
    }

    // An unknown channel stops the job.
    {
        bool threw = false;
        mf::Counters k;
        try {
            (void)mf::digitise({tube_hit(99999999, {0, 0, 0}, {0, 0, 1})}, wires, config, rng, k);
        } catch (std::runtime_error const&) {
            threw = true;
        }
        check(threw, "a channel missing from the wire table");
    }

    if (failures == 0) {
        std::puts("muflux_digitiser: all checks passed");
    }
    return failures == 0 ? 0 : 1;
}
