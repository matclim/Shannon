// SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
//
// SPDX-License-Identifier: LGPL-3.0-or-later

// muflux_digitise.cpp — Phlex module plugin
//
// Digitises the drift tubes and RPC strips of the 2018 muon-flux setup from
// the sim_hits of aegir's workflows/muflux2018_replay.jsonnet, one event per
// input entry (see muflux_digitiser.hpp for the model). A counter-based RNG is
// seeded per event from (seed, entry), so results do not depend on the number
// of threads.
//
// Configuration:
//   layer           layer of the input events (required)
//   wire_table      dt_wires.txt exported from the legacy geometry (required)
//   seed            RNG seed (default 0)
//   sigma_spatial   drift-distance resolution [mm] (default 0.27)
//   drift_time_per_mm  [ns/mm] (default 72)
//   dt_efficiency, rpc_efficiency  per hit (default 1: off, as in the legacy simulation)
//   dead_channels   drift-tube ids without a hit (default: the legacy simulation's list)
//   step_segments   true (default): hits are Geant4 steps with an energy deposit (aegir
//                   sd_mode 'scoring'); false: entry points of tracks (sd_mode 'crossing')

#include "muflux/muflux_digitiser.hpp"
#include "phlex/core/product_selector.hpp"
#include "phlex/model/data_cell_index.hpp"
#include "phlex/module.hpp"

#include <SHiP/SimHit.hpp>
#include <SHiP/random/philox_rng.hpp>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <spdlog/spdlog.h>
#include <stdexcept>
#include <string>
#include <vector>

using namespace phlex;

namespace {

namespace mf = SHiP::shannon::muflux;

constexpr std::uint32_t muflux_digitise_stream = 0x3F1A2018;

// Adds up the counters of all events and prints them when the job ends.
class Summary {
   public:
    void add(mf::Counters const& c) {
        std::scoped_lock lock{mutex_};
        total_ += c;
        ++events_;
    }
    ~Summary() {
        try {
            spdlog::info(
                "muflux_digitise: {} events; drift tubes: {} hits of charged particles ({} of "
                "neutral ones dropped), {} in dead channels, {} inefficient, {} later in the same "
                "tube, {} with the distance set to the gas radius; RPC: {} strip hits ({} neutral "
                "dropped), {} inefficient, {} on a strip already fired",
                events_, total_.dt_hits, total_.dt_neutral, total_.dt_dead, total_.dt_inefficient,
                total_.dt_later, total_.dt_clipped, total_.rpc_hits, total_.rpc_neutral,
                total_.rpc_inefficient, total_.rpc_shared);
        } catch (...) {
        }
    }

   private:
    std::mutex mutex_;
    mf::Counters total_;
    std::uint64_t events_ = 0;
};

}  // namespace

PHLEX_REGISTER_ALGORITHMS(m, config) {
    auto const layer = phlex::experimental::identifier{config.get<std::string>("layer")};
    auto const seed = static_cast<std::uint32_t>(config.get<int>("seed", 0));
    auto const wires = std::make_shared<std::map<std::int32_t, mf::Wire> const>(
        mf::read_wire_table(config.get<std::string>("wire_table")));

    mf::DigitiserConfig model;
    model.sigma_spatial_mm = config.get<double>("sigma_spatial", double{model.sigma_spatial_mm});
    model.drift_time_per_mm =
        config.get<double>("drift_time_per_mm", double{model.drift_time_per_mm});
    model.dt_efficiency = config.get<double>("dt_efficiency", double{model.dt_efficiency});
    model.rpc_efficiency = config.get<double>("rpc_efficiency", double{model.rpc_efficiency});
    model.step_segments = config.get<bool>("step_segments", bool{model.step_segments});
    if (auto const dead = config.get_if_present<std::vector<int>>("dead_channels")) {
        model.dead_channels = std::set<std::int32_t>(dead->begin(), dead->end());
    }
    for (auto const& id : model.dead_channels) {
        if (!wires->contains(id)) {
            throw std::runtime_error("muflux_digitise: dead channel " + std::to_string(id) +
                                     " is not in the wire table");
        }
    }
    spdlog::info(
        "muflux_digitise: {} wires, sigma {} mm, {} ns/mm, efficiencies {} (tubes) and {} (RPC), "
        "{} dead channels, seed {}",
        wires->size(), model.sigma_spatial_mm, model.drift_time_per_mm, model.dt_efficiency,
        model.rpc_efficiency, model.dead_channels.size(), seed);

    auto const summary = std::make_shared<Summary>();

    m.transform(
         "muflux_digitise",
         [seed, wires, model, summary](data_cell_index const& id,
                                       std::vector<::SHiP::SimHit> const& sim_hits) {
             std::vector<mf::HitView> hits;
             hits.reserve(sim_hits.size());
             for (auto const& h : sim_hits) {
                 hits.push_back({.detector_id = h.detector_id,
                                 .channel = h.geometry_node_id,
                                 .track_id = h.track_id,
                                 .pdg_code = h.pdg_code,
                                 .position = {h.position[0], h.position[1], h.position[2]},
                                 .momentum = {h.momentum[0], h.momentum[1], h.momentum[2]},
                                 .time = h.time,
                                 .path_length = h.path_length});
             }
             SHiP::random::PhiloxRng rng{seed, muflux_digitise_stream, id.number()};
             mf::Counters counters;
             auto digits = mf::digitise(hits, *wires, model, rng, counters);
             summary->add(counters);
             return digits;
         },
         concurrency::unlimited)
        .input_family(
            product_selector{.creator = "rntuple_source", .layer = layer, .suffix = "id"},
            product_selector{.creator = "rntuple_source", .layer = layer, .suffix = "sim_hits"})
        .output_product_suffixes("muflux_digits");
}
