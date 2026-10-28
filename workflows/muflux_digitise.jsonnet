// SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
//
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// 2018 muon-flux test beam: digitise the drift tubes and RPC strips of
// aegir's workflows/muflux2018_replay.jsonnet output, one event per entry,
// into the muflux_digits RNTuple (src/muflux/muflux_digits_output.cpp).
//
// Top-level arguments (--tla-str / --tla-code):
//   infile      aegir replay output (required)
//   events      entries to digitise, at most those of infile (required)
//   wire_table  dt_wires.txt exported from the legacy geometry (required)
//   out         output file (default 'muflux_digits.root')
//   seed        RNG seed (default 0)
//   step_segments  true (default) for aegir sd_mode 'scoring' input, false for 'crossing'
function(infile, events, wire_table, out='muflux_digits.root', seed=0, step_segments=true)
  {
    stage: 'digitisation',
    driver: {
      cpp: 'generate_layers',
      layers: {
        event: { parent: 'job', total: events },
      },
    },
    sources: {
      rntuple_source: {
        cpp: 'read_sim_file',
        input_file: infile,
        ntuple_name: 'events',
        layer: 'event',
      },
    },
    modules: {
      muflux_digitise: {
        cpp: 'muflux_digitise',
        layer: 'event',
        wire_table: wire_table,
        seed: seed,
        step_segments: step_segments,
      },
      output: {
        cpp: 'muflux_digits_output',
        layer: 'event',
        rntuple_file: out,
      },
    },
  }
