// SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
//
// SPDX-License-Identifier: LGPL-3.0-or-later

// muflux_digits_output.cpp — Phlex module plugin
//
// Writes the digits of the 2018 muon-flux setup (muflux_digitise) to an
// RNTuple named `muflux_digits`, one entry per event. Every field is a
// fundamental type or a std::vector of one, so the file is read with plain
// ROOT, PyROOT or uproot, without dictionaries. The legacy data digits are
// converted to the same layout, so reconstruction treats data and simulation
// alike.
//
// Fields:
//   entry             entry of the event in the simulation input
//   dt_id             drift-tube channel id (legacy numbering)
//   dt_drift          drift distance [mm] (smeared)
//   dt_drift_time     drift time [ns]
//   dt_true_dist      true distance of closest approach [mm] (simulation)
//   dt_track, dt_pdg  Geant4 track id and PDG code of the hit (simulation)
//   rpc_id            RPC strip id (legacy numbering)
//   rpc_track         track id of the strip's first hit (simulation)
//   n_sa, n_sb        tracks entering the trigger scintillators (simulation)
//
// Events finish in any order under concurrency, so entries are not sorted by
// `entry`.
//
// Configuration: rntuple_file (default muflux_digits.root), layer (required),
// creator (default muflux_digitise).

#include "muflux/muflux_digitiser.hpp"
#include "phlex/core/product_selector.hpp"
#include "phlex/model/data_cell_index.hpp"
#include "phlex/module.hpp"

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>

#include <cstdint>
#include <memory>
#include <spdlog/spdlog.h>
#include <string>
#include <vector>

using namespace phlex;

namespace {

namespace mf = SHiP::shannon::muflux;

class DigitsWriter {
   public:
    explicit DigitsWriter(std::string const& filename) {
        auto model = ROOT::RNTupleModel::Create();
        entry_ = model->MakeField<std::uint64_t>("entry");
        dt_id_ = model->MakeField<std::vector<std::int32_t>>("dt_id");
        dt_drift_ = model->MakeField<std::vector<float>>("dt_drift");
        dt_drift_time_ = model->MakeField<std::vector<float>>("dt_drift_time");
        dt_true_dist_ = model->MakeField<std::vector<float>>("dt_true_dist");
        dt_track_ = model->MakeField<std::vector<std::int32_t>>("dt_track");
        dt_pdg_ = model->MakeField<std::vector<std::int32_t>>("dt_pdg");
        rpc_id_ = model->MakeField<std::vector<std::int32_t>>("rpc_id");
        rpc_track_ = model->MakeField<std::vector<std::int32_t>>("rpc_track");
        n_sa_ = model->MakeField<std::int32_t>("n_sa");
        n_sb_ = model->MakeField<std::int32_t>("n_sb");
        writer_ = ROOT::RNTupleWriter::Recreate(std::move(model), "muflux_digits", filename);
        spdlog::info("muflux_digits_output: writing '{}'", filename);
    }

    // Called serially (concurrency::serial): the writer is not thread-safe.
    void write(data_cell_index const& id, mf::Digits const& d) {
        *entry_ = id.number();
        *dt_id_ = d.dt_id;
        *dt_drift_ = d.dt_drift;
        *dt_drift_time_ = d.dt_drift_time;
        *dt_true_dist_ = d.dt_true_dist;
        *dt_track_ = d.dt_track;
        *dt_pdg_ = d.dt_pdg;
        *rpc_id_ = d.rpc_id;
        *rpc_track_ = d.rpc_track;
        *n_sa_ = d.n_sa;
        *n_sb_ = d.n_sb;
        writer_->Fill();
    }

   private:
    std::shared_ptr<std::uint64_t> entry_;
    std::shared_ptr<std::vector<std::int32_t>> dt_id_, dt_track_, dt_pdg_, rpc_id_, rpc_track_;
    std::shared_ptr<std::vector<float>> dt_drift_, dt_drift_time_, dt_true_dist_;
    std::shared_ptr<std::int32_t> n_sa_, n_sb_;
    std::unique_ptr<ROOT::RNTupleWriter> writer_;  // last: written out before the fields go
};

}  // namespace

PHLEX_REGISTER_ALGORITHMS(m, config) {
    auto const rntuple_file =
        config.get<std::string>("rntuple_file", std::string{"muflux_digits.root"});
    auto const layer = phlex::experimental::identifier{config.get<std::string>("layer")};
    auto const creator = phlex::experimental::identifier{
        config.get<std::string>("creator", std::string{"muflux_digitise"})};

    auto writer = m.make<DigitsWriter>(rntuple_file);
    writer.observe("write_muflux_digits", &DigitsWriter::write, concurrency::serial)
        .input_family(
            product_selector{.creator = "rntuple_source", .layer = layer, .suffix = "id"},
            product_selector{.creator = creator, .layer = layer, .suffix = "muflux_digits"});
}
