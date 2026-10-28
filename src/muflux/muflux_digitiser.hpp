// SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
//
// SPDX-License-Identifier: LGPL-3.0-or-later

// muflux_digitiser.hpp — digitisation of the 2018 muon-flux spectrometer
//
// The 2018 muon-flux test beam (SPS H4) had four drift-tube stations (T1, T2
// with x-u and v-x views, T3, T4 with x views; 576 tubes) and an RPC muon
// tagger (5 stations, horizontal and vertical strips). aegir's
// workflows/muflux2018_replay.jsonnet records one SimHit where a track enters
// a tube's gas or an RPC strip, with the legacy channel id as
// geometry_node_id:
//   drift tube: station*1e7 + view*1e6 + plane*1e5 + layer*1e4 + 2000 + tube
//   RPC strip:  station*1e4 + direction*1e3 + strip (direction 0: H, 1: V)
//
// The model follows the legacy FairShip digitisation of this setup (FairShip
// branch muflux: python/MufluxDigi.py, charmdet/MufluxSpectrometerHit.cxx):
//   - the signal of a tube comes from the ionisation closest to its wire:
//     each hit is a Geant4 step with an energy deposit (aegir's sd_mode
//     'scoring': start point, momentum direction, step length), the distance
//     is that of the step segment to the wire, and of all the hits of a tube
//     (the track and its delta electrons, other tracks) the earliest signal,
//     hit time + distance / v_drift, is kept. With step_segments = false a
//     hit is the entry point of a track (sd_mode 'crossing') and the distance
//     that of its straight line to the wire;
//   - drift distance = |Gauss(that distance, sigma_spatial)|, drift time =
//     drift distance / v_drift;
//   - the tubes in dead_channels give no hit;
//   - an RPC strip fires once per event, whatever the number of tracks;
//   - only charged particles make hits. aegir's sd_mode 'crossing' records
//     every particle entering a volume, neutrinos and photons included; a
//     neutrino from a pi/K decay flies along its muon through T1/T2 and would
//     otherwise take its tubes.
// Inefficiencies are off by default, as in the legacy simulation.
//
// Kept free of phlex and ROOT so that it can be tested on its own.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace SHiP::shannon::muflux {

using Vec3 = std::array<double, 3>;

inline Vec3 sub(Vec3 const& a, Vec3 const& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
inline double dot(Vec3 const& a, Vec3 const& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline Vec3 cross(Vec3 const& a, Vec3 const& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
inline double norm(Vec3 const& a) {
    return std::sqrt(dot(a, a));
}

/// One drift tube: the two ends of its wire and the radius of its gas, in mm.
struct Wire {
    Vec3 end1{};
    Vec3 end2{};
    double gas_radius = 0;
};

/// Drift-tube channel id -> wire, read from the text table exported from the
/// legacy geometry (dt_wires.txt): one line per tube,
///   id x1 y1 z1 x2 y2 z2 r_gas      (lengths in cm)
/// with '#' starting a comment line.
inline std::map<std::int32_t, Wire> read_wire_table(std::string const& path) {
    std::ifstream in{path};
    if (!in) {
        throw std::runtime_error("muflux: cannot open wire table '" + path + "'");
    }
    std::map<std::int32_t, Wire> wires;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        auto const first = line.find_first_not_of(" \t");
        if (first == std::string::npos || line[first] == '#') {
            continue;
        }
        std::istringstream fields{line};
        std::int64_t id = 0;
        std::array<double, 7> v{};
        fields >> id;
        for (auto& x : v) {
            fields >> x;
        }
        if (!fields) {
            throw std::runtime_error("muflux: wire table '" + path + "', line " +
                                     std::to_string(line_number) + ": expected 8 numbers");
        }
        constexpr double cm = 10.0;  // the table is in cm, the hits in mm
        Wire const w{.end1 = {v[0] * cm, v[1] * cm, v[2] * cm},
                     .end2 = {v[3] * cm, v[4] * cm, v[5] * cm},
                     .gas_radius = v[6] * cm};
        if (!(norm(sub(w.end2, w.end1)) > 0) || !(w.gas_radius > 0)) {
            throw std::runtime_error("muflux: wire table '" + path + "', line " +
                                     std::to_string(line_number) + ": zero-length wire or radius");
        }
        if (!wires.emplace(static_cast<std::int32_t>(id), w).second) {
            throw std::runtime_error("muflux: wire table '" + path + "', line " +
                                     std::to_string(line_number) + ": duplicate id " +
                                     std::to_string(id));
        }
    }
    if (wires.empty()) {
        throw std::runtime_error("muflux: wire table '" + path + "' has no wires");
    }
    return wires;
}

/// Distance of closest approach between the straight track through `point`
/// along `direction` and the wire (both taken as infinite lines), in mm.
inline double distance_to_wire(Vec3 const& point, Vec3 const& direction, Wire const& wire) {
    Vec3 const w = sub(wire.end2, wire.end1);
    Vec3 const d = sub(point, wire.end1);
    Vec3 const n = cross(direction, w);
    double const n_norm = norm(n);
    double const w_norm = norm(w);
    double const dir_norm = norm(direction);
    // A track (nearly) parallel to the wire: its distance from the wire line.
    if (dir_norm == 0 || n_norm <= 1e-12 * dir_norm * w_norm) {
        return norm(cross(d, w)) / w_norm;
    }
    return std::abs(dot(d, n)) / n_norm;
}

/// Distance between the segment from `point` along `direction` over
/// `length` and the wire (an infinite line), in mm. The distance to the line
/// is convex along the segment, so the closest point of the segment is the
/// line's closest point clamped to the segment.
inline double distance_segment_to_wire(Vec3 const& point, Vec3 const& direction, double length,
                                       Wire const& wire) {
    Vec3 const w = sub(wire.end2, wire.end1);
    double const wn = norm(w);
    Vec3 const wu{w[0] / wn, w[1] / wn, w[2] / wn};
    double const dn = norm(direction);
    auto const distance_of = [&](Vec3 const& q) { return norm(cross(sub(q, wire.end1), wu)); };
    if (dn == 0 || length <= 0) {
        return distance_of(point);
    }
    Vec3 const u{direction[0] / dn, direction[1] / dn, direction[2] / dn};
    // closest point of the infinite track line to the wire line
    Vec3 const r = sub(point, wire.end1);
    double const b = dot(u, wu);
    double const denom = 1 - b * b;
    double s = 0;
    if (denom > 1e-12) {
        s = (b * dot(r, wu) - dot(r, u)) / denom;
    }
    s = std::clamp(s, 0.0, length);
    return distance_of({point[0] + s * u[0], point[1] + s * u[1], point[2] + s * u[2]});
}

/// Three times the electric charge of a particle from its PDG code: quarks,
/// leptons, gauge bosons, mesons, baryons (quark content) and nuclei.
inline int three_charge(std::int32_t pdg) {
    constexpr int quark[7] = {0, -1, 2, -1, 2, -1, 2};  // d u s c b t
    int const a = pdg < 0 ? -pdg : pdg;
    int const sign = pdg < 0 ? -1 : 1;
    int q = 0;
    if (a >= 1000000000) {  // nucleus 10LZZZAAAI
        q = 3 * (a / 10000 % 1000);
    } else if (a < 100) {
        if (a >= 1 && a <= 6) {
            q = quark[a];
        } else if (a == 11 || a == 13 || a == 15 || a == 17) {
            q = -3;
        } else if (a == 24 || a == 37) {
            q = 3;
        }
    } else {
        int const n1 = a / 1000 % 10, n2 = a / 100 % 10, n3 = a / 10 % 10;
        if (n1 == 0) {  // meson n2 n3
            if (n2 <= 6 && n3 <= 6) {
                q = (n2 == 3 || n2 == 5) ? quark[n3] - quark[n2] : quark[n2] - quark[n3];
            }
        } else if (n1 <= 6 && n2 <= 6 && n3 <= 6) {  // baryon n1 n2 n3
            q = quark[n1] + quark[n2] + quark[n3];
        }
    }
    return sign * q;
}

/// Model parameters, defaults as in the legacy FairShip simulation
/// (geometry/charm-geometry_config.py, python/MufluxDigi.py).
struct DigitiserConfig {
    double sigma_spatial_mm = 0.27;  ///< 0.027 cm
    double drift_time_per_mm = 72.;  ///< ns/mm: v_drift = 1/(72 ns/mm)
    double dt_efficiency = 1.0;      ///< per tube hit
    double rpc_efficiency = 1.0;     ///< per strip hit
    bool step_segments = true;       ///< hits are steps (sd_mode 'scoring'), else entry points
    std::set<std::int32_t> dead_channels{10112001, 20112003, 30002041, 30012026,
                                         30102025, 30112013, 30112018, 40012014};
};

/// The digits of one event; the vectors of each detector are parallel.
struct Digits {
    // drift tubes, ordered by channel id
    std::vector<std::int32_t> dt_id;
    std::vector<float> dt_drift;         ///< smeared drift distance [mm]
    std::vector<float> dt_drift_time;    ///< drift distance / v_drift [ns]
    std::vector<float> dt_true_dist;     ///< true distance of closest approach [mm]
    std::vector<std::int32_t> dt_track;  ///< Geant4 track id of the hit kept
    std::vector<std::int32_t> dt_pdg;
    // RPC strips, ordered by channel id
    std::vector<std::int32_t> rpc_id;
    std::vector<std::int32_t> rpc_track;  ///< track id of the first hit of the strip
    // trigger scintillators: charged tracks crossing them
    std::int32_t n_sa = 0;
    std::int32_t n_sb = 0;
};

/// Event and job counters, for the summary at the end of a job.
struct Counters {
    std::uint64_t dt_hits = 0;         ///< simulated drift-tube hits of charged particles
    std::uint64_t dt_neutral = 0;      ///< of neutral particles, dropped
    std::uint64_t dt_dead = 0;         ///< tubes hit in a dead channel
    std::uint64_t dt_inefficient = 0;  ///< tubes dropped by dt_efficiency
    std::uint64_t dt_later = 0;        ///< further hits in a tube, with a later signal
    std::uint64_t dt_clipped = 0;      ///< distance beyond the gas radius, set to it
    std::uint64_t rpc_hits = 0;
    std::uint64_t rpc_neutral = 0;
    std::uint64_t rpc_inefficient = 0;
    std::uint64_t rpc_shared = 0;  ///< further hits on a strip already fired

    Counters& operator+=(Counters const& o) {
        dt_hits += o.dt_hits;
        dt_neutral += o.dt_neutral;
        rpc_neutral += o.rpc_neutral;
        dt_dead += o.dt_dead;
        dt_inefficient += o.dt_inefficient;
        dt_later += o.dt_later;
        dt_clipped += o.dt_clipped;
        rpc_hits += o.rpc_hits;
        rpc_inefficient += o.rpc_inefficient;
        rpc_shared += o.rpc_shared;
        return *this;
    }
};

/// What the digitiser needs from a simulated hit, so that the model can be
/// tested without the data model.
struct HitView {
    std::int32_t detector_id;
    std::int32_t channel;  ///< SimHit::geometry_node_id
    std::int32_t track_id;
    std::int32_t pdg_code;
    Vec3 position;           ///< mm
    Vec3 momentum;           ///< GeV/c
    double time;             ///< ns
    double path_length = 0;  ///< step length [mm] (sd_mode 'scoring')
};

constexpr std::int32_t kDriftTubeDetectorId = 2;  // SHiP::detector_id::StrawTubes
constexpr std::int32_t kOtherDetectorId = 99;     // aegir's test_fixture: RPC strips, SA, SB
constexpr double kSaSbBoundaryMm = 700.;          // SA at z = 115 mm, SB at 1363 mm

/// Digitise the hits of one event. `rng` provides uniform53(lo, hi) and
/// gaussian(mean, sigma).
template <typename Rng>
Digits digitise(std::vector<HitView> const& hits, std::map<std::int32_t, Wire> const& wires,
                DigitiserConfig const& config, Rng& rng, Counters& counters) {
    struct Signal {
        double time;  // hit time + drift time of its closest ionisation
        double dist;  // true distance of that ionisation to the wire
        std::int32_t track;
        std::int32_t pdg;
    };
    std::map<std::int32_t, Signal> tubes;         // earliest signal per tube
    std::map<std::int32_t, std::int32_t> strips;  // id -> track of the first hit
    std::set<std::int32_t> sa_tracks, sb_tracks, dead;
    Digits digits;

    for (auto const& h : hits) {
        bool const neutral = three_charge(h.pdg_code) == 0;
        if (h.detector_id == kDriftTubeDetectorId) {
            if (neutral) {
                ++counters.dt_neutral;
                continue;
            }
            ++counters.dt_hits;
            auto const wire = wires.find(h.channel);
            if (wire == wires.end()) {
                throw std::runtime_error(
                    "muflux: drift-tube hit with channel " + std::to_string(h.channel) +
                    " that is not in the wire table (is copy_number_depth set in aegir?)");
            }
            if (config.dead_channels.contains(h.channel)) {
                dead.insert(h.channel);
                continue;
            }
            double dist =
                config.step_segments
                    ? distance_segment_to_wire(h.position, h.momentum, h.path_length, wire->second)
                    : distance_to_wire(h.position, h.momentum, wire->second);
            if (dist > wire->second.gas_radius) {
                ++counters.dt_clipped;
                dist = wire->second.gas_radius;
            }
            Signal const sig{.time = h.time + dist * config.drift_time_per_mm,
                             .dist = dist,
                             .track = h.track_id,
                             .pdg = h.pdg_code};
            auto [it, inserted] = tubes.try_emplace(h.channel, sig);
            if (!inserted) {
                ++counters.dt_later;
                if (sig.time < it->second.time) {
                    it->second = sig;
                }
            }
        } else if (h.detector_id == kOtherDetectorId) {
            if (neutral) {
                if (h.channel > 0) {
                    ++counters.rpc_neutral;
                }
                continue;
            }
            if (h.channel > 0) {
                ++counters.rpc_hits;
                if (config.rpc_efficiency < 1.0 &&
                    rng.uniform53(0.0, 1.0) >= config.rpc_efficiency) {
                    ++counters.rpc_inefficient;
                    continue;
                }
                if (!strips.try_emplace(h.channel, h.track_id).second) {
                    ++counters.rpc_shared;
                }
            } else {
                (h.position[2] < kSaSbBoundaryMm ? sa_tracks : sb_tracks).insert(h.track_id);
            }
        }
    }

    // inefficiency and smearing per tube, in channel order
    counters.dt_dead += dead.size();
    for (auto const& [id, sig] : tubes) {
        if (config.dt_efficiency < 1.0 && rng.uniform53(0.0, 1.0) >= config.dt_efficiency) {
            ++counters.dt_inefficient;
            continue;
        }
        double const drift = std::abs(rng.gaussian(sig.dist, config.sigma_spatial_mm));
        digits.dt_id.push_back(id);
        digits.dt_drift.push_back(static_cast<float>(drift));
        digits.dt_drift_time.push_back(static_cast<float>(drift * config.drift_time_per_mm));
        digits.dt_true_dist.push_back(static_cast<float>(sig.dist));
        digits.dt_track.push_back(sig.track);
        digits.dt_pdg.push_back(sig.pdg);
    }
    digits.n_sa = static_cast<std::int32_t>(sa_tracks.size());
    digits.n_sb = static_cast<std::int32_t>(sb_tracks.size());
    for (auto const& [id, track] : strips) {
        digits.rpc_id.push_back(id);
        digits.rpc_track.push_back(track);
    }
    return digits;
}

}  // namespace SHiP::shannon::muflux
