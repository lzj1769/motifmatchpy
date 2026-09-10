// nanobind bindings for the vendored MOODS core.
//
// This replaces the SWIG interface files (core/*.i) that upstream MOODS uses.
// The exposed surface mirrors the upstream MOODS.tools / MOODS.scan /
// MOODS.parsers / MOODS.misc modules one-to-one, so anything written against
// MOODS keeps working; the friendlier layer lives in Python on top of this.

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>

#include <sstream>
#include <stdexcept>

#include "moods.h"
#include "match_types.h"
#include "moods_misc.h"
#include "moods_parsers.h"
#include "moods_scan.h"
#include "moods_tools.h"
#include "scanner.h"

namespace nb = nanobind;
using namespace nb::literals;

using MOODS::match;
using MOODS::match_with_variant;
using MOODS::variant;
using MOODS::scan::Scanner;

namespace {

// Upstream returns an empty matrix for anything it fails to parse. Turn that
// into a real exception so callers get a diagnosable failure instead of a
// silently empty result.
score_matrix check_parsed(score_matrix mat, const std::string &filename,
                          const char *what) {
    if (mat.empty()) {
        throw std::runtime_error("could not parse " + std::string(what) +
                                 " matrix file: " + filename);
    }
    return mat;
}

std::string repr_double(double v) {
    std::ostringstream os;
    os.precision(12);
    os << v;
    return os.str();
}

}  // namespace

NB_MODULE(_core, m) {
    m.doc() =
        "Low-level bindings to the MOODS C++ core (vendored from "
        "github.com/jhkorhonen/MOODS).";
    m.attr("MOODS_VERSION") = "1.9.4.1";

    // ------------------------------------------------------------------
    // match types
    // ------------------------------------------------------------------
    nb::class_<match>(m, "Match",
                      "A motif occurrence: 0-based start position and score.")
        .def(nb::init<>())
        .def("__init__",
             [](match *self, size_t pos, double score) {
                 new (self) match{pos, score};
             },
             "pos"_a, "score"_a)
        .def_rw("pos", &match::pos, "0-based start position in the sequence.")
        .def_rw("score", &match::score, "Score of the occurrence.")
        .def("__iter__",
             [](const match &s) {
                 return nb::iter(nb::make_tuple(s.pos, s.score));
             })
        .def("__len__", [](const match &) { return 2; })
        .def("__eq__",
             [](const match &a, const match &b) {
                 return a.pos == b.pos && a.score == b.score;
             },
             nb::is_operator())
        .def("__repr__", [](const match &s) {
            return "Match(pos=" + std::to_string(s.pos) +
                   ", score=" + repr_double(s.score) + ")";
        });

    nb::class_<variant>(
        m, "Variant",
        "A sequence variant: [start_pos, end_pos) replaced by modified_seq.\n"
        "end == start + 1 is a substitution, end == start an insertion, and an\n"
        "empty modified_seq a deletion.")
        .def(nb::init<>())
        .def(nb::init<size_t, size_t, std::string>(), "start_pos"_a,
             "end_pos"_a, "modified_seq"_a)
        .def_rw("start_pos", &variant::start_pos)
        .def_rw("end_pos", &variant::end_pos)
        .def_rw("modified_seq", &variant::modified_seq)
        .def("__lt__", [](const variant &a,
                          const variant &b) { return a < b; }, nb::is_operator())
        .def("__eq__",
             [](const variant &a, const variant &b) {
                 return a.start_pos == b.start_pos && a.end_pos == b.end_pos &&
                        a.modified_seq == b.modified_seq;
             },
             nb::is_operator())
        .def("__repr__", [](const variant &v) {
            return "Variant(start_pos=" + std::to_string(v.start_pos) +
                   ", end_pos=" + std::to_string(v.end_pos) +
                   ", modified_seq='" + v.modified_seq + "')";
        });

    nb::class_<match_with_variant>(
        m, "MatchWithVariant",
        "A motif occurrence that depends on one or more applied variants.\n"
        "`variants` holds indices into the variant list passed to the scanner.")
        .def(nb::init<>())
        .def_rw("pos", &match_with_variant::pos)
        .def_rw("score", &match_with_variant::score)
        .def_rw("variants", &match_with_variant::variants)
        .def("__repr__", [](const match_with_variant &s) {
            std::string v;
            for (size_t i = 0; i < s.variants.size(); ++i) {
                if (i) v += ", ";
                v += std::to_string(s.variants[i]);
            }
            return "MatchWithVariant(pos=" + std::to_string(s.pos) +
                   ", score=" + repr_double(s.score) + ", variants=[" + v + "])";
        });

    // ------------------------------------------------------------------
    // tools
    // ------------------------------------------------------------------
    nb::module_ tools = m.def_submodule("tools", "Matrix and background tools.");
    namespace T = MOODS::tools;

    tools.attr("DEFAULT_DP_PRECISION") = T::DEFAULT_DP_PRECISION;

    tools.def("flat_bg", &T::flat_bg, "alphabet_size"_a,
              "Uniform background distribution over `alphabet_size` symbols.");
    tools.def("bg_from_sequence_dna", &T::bg_from_sequence_dna, "seq"_a, "ps"_a,
              nb::call_guard<nb::gil_scoped_release>(),
              "Estimate an ACGT background distribution from a DNA sequence.");
    tools.def("snp_variants", &T::snp_variants, "seq"_a,
              nb::call_guard<nb::gil_scoped_release>(),
              "Variants for every IUPAC ambiguity code in the sequence.");

    tools.def("reverse_complement",
              static_cast<score_matrix (*)(const score_matrix &)>(
                  &T::reverse_complement),
              "mat"_a, "Reverse complement of a 0-order matrix.");
    tools.def("reverse_complement",
              static_cast<score_matrix (*)(const score_matrix &, size_t)>(
                  &T::reverse_complement),
              "mat"_a, "a"_a,
              "Reverse complement of a matrix over an alphabet of size `a`\n"
              "(handles both 0-order and high-order matrices).");

    tools.def("log_odds",
              static_cast<score_matrix (*)(const score_matrix &,
                                           const std::vector<double> &, double)>(
                  &T::log_odds),
              "mat"_a, "bg"_a, "ps"_a,
              "Log-odds transform of a count/frequency matrix (natural log).");
    tools.def("log_odds",
              static_cast<score_matrix (*)(const score_matrix &,
                                           const std::vector<double> &, double,
                                           double)>(&T::log_odds),
              "mat"_a, "bg"_a, "ps"_a, "log_base"_a,
              "Log-odds transform of a count/frequency matrix.");
    tools.def("log_odds",
              static_cast<score_matrix (*)(
                  const score_matrix &, const std::vector<std::vector<double>> &,
                  const std::vector<double> &, double, size_t)>(&T::log_odds),
              "mat"_a, "low_order_terms"_a, "bg"_a, "ps"_a, "a"_a,
              "High-order log-odds transform (natural log).");
    tools.def("log_odds",
              static_cast<score_matrix (*)(
                  const score_matrix &, const std::vector<std::vector<double>> &,
                  const std::vector<double> &, double, size_t, double)>(
                  &T::log_odds),
              "mat"_a, "low_order_terms"_a, "bg"_a, "ps"_a, "a"_a, "log_base"_a,
              "High-order log-odds transform.");

    tools.def("threshold_from_p",
              static_cast<double (*)(const score_matrix &,
                                     const std::vector<double> &, const double &)>(
                  &T::threshold_from_p),
              "pssm"_a, "bg"_a, "p"_a,
              nb::call_guard<nb::gil_scoped_release>(),
              "Score threshold corresponding to p-value `p`.");
    tools.def("threshold_from_p",
              static_cast<double (*)(const score_matrix &,
                                     const std::vector<double> &, const double &,
                                     size_t)>(&T::threshold_from_p),
              "pssm"_a, "bg"_a, "p"_a, "a"_a,
              nb::call_guard<nb::gil_scoped_release>(),
              "Score threshold for a high-order matrix over an alphabet of size `a`.");
    tools.def("threshold_from_p_with_precision",
              static_cast<double (*)(const score_matrix &,
                                     const std::vector<double> &, const double &,
                                     double)>(&T::threshold_from_p_with_precision),
              "pssm"_a, "bg"_a, "p"_a, "precision"_a,
              nb::call_guard<nb::gil_scoped_release>(),
              "Score threshold from p-value at an explicit DP precision.");
    tools.def("threshold_from_p_with_precision",
              static_cast<double (*)(const score_matrix &,
                                     const std::vector<double> &, const double &,
                                     double, size_t)>(
                  &T::threshold_from_p_with_precision),
              "pssm"_a, "bg"_a, "p"_a, "precision"_a, "a"_a,
              nb::call_guard<nb::gil_scoped_release>(),
              "High-order score threshold from p-value at an explicit DP precision.");

    tools.def("max_score",
              static_cast<double (*)(const score_matrix &)>(&T::max_score),
              "mat"_a, "Highest score the matrix can produce.");
    tools.def("max_score",
              static_cast<double (*)(const score_matrix &, size_t)>(&T::max_score),
              "mat"_a, "a"_a, "Highest score of a high-order matrix.");
    tools.def("min_score",
              static_cast<double (*)(const score_matrix &)>(&T::min_score),
              "mat"_a, "Lowest score the matrix can produce.");
    tools.def("min_score",
              static_cast<double (*)(const score_matrix &, size_t)>(&T::min_score),
              "mat"_a, "a"_a, "Lowest score of a high-order matrix.");
    tools.def("min_delta", &T::min_delta, "mat"_a,
              "Smallest non-zero score difference within a matrix column.");

    // ------------------------------------------------------------------
    // parsers
    // ------------------------------------------------------------------
    nb::module_ parsers = m.def_submodule("parsers", "Matrix file parsers.");
    namespace P = MOODS::parsers;

    parsers.def(
        "pfm",
        [](const std::string &filename) {
            return check_parsed(P::pfm(filename), filename, "pfm");
        },
        "filename"_a, "Read a JASPAR-style count matrix as-is.");
    parsers.def(
        "pfm_to_log_odds",
        [](const std::string &filename, const std::vector<double> &bg,
           double pseudocount, double log_base) {
            return check_parsed(
                P::pfm_to_log_odds(filename, bg, pseudocount, log_base),
                filename, "pfm");
        },
        "filename"_a, "bg"_a, "pseudocount"_a, "log_base"_a = -1.0,
        "Read a count matrix and log-odds transform it.");
    parsers.def(
        "adm_1o_terms",
        [](const std::string &filename, size_t a) {
            return check_parsed(P::adm_1o_terms(filename, a), filename, "adm");
        },
        "filename"_a, "a"_a = 4,
        "First-order terms of an adjacent dinucleotide model file.");
    parsers.def(
        "adm_0o_terms",
        [](const std::string &filename, size_t a) {
            return check_parsed(P::adm_0o_terms(filename, a), filename, "adm");
        },
        "filename"_a, "a"_a = 4,
        "Zero-order terms of an adjacent dinucleotide model file.");
    parsers.def(
        "adm_to_log_odds",
        [](const std::string &filename, const std::vector<double> &bg,
           double pseudocount, size_t a, double log_base) {
            return check_parsed(
                P::adm_to_log_odds(filename, bg, pseudocount, a, log_base),
                filename, "adm");
        },
        "filename"_a, "bg"_a, "pseudocount"_a, "a"_a = 4, "log_base"_a = -1.0,
        "Read an adjacent dinucleotide model and log-odds transform it.");

    // ------------------------------------------------------------------
    // scan
    // ------------------------------------------------------------------
    nb::module_ scan = m.def_submodule("scan", "Scanning entry points.");
    namespace S = MOODS::scan;

    scan.def("scan_dna", &S::scan_dna, "seq"_a, "matrices"_a, "bg"_a,
             "thresholds"_a, "window_size"_a = 7,
             nb::call_guard<nb::gil_scoped_release>(),
             "Scan a DNA sequence with the given matrices and thresholds.");
    scan.def("scan", &S::scan, "seq"_a, "matrices"_a, "bg"_a, "thresholds"_a,
             "window_size"_a, "alphabet"_a,
             nb::call_guard<nb::gil_scoped_release>(),
             "Scan a sequence over a custom alphabet.");
    scan.def("scan_best_hits_dna", &S::scan_best_hits_dna, "seq"_a, "matrices"_a,
             "target"_a, "iterations"_a = 10, "MULT"_a = 2, "LIMIT_MULT"_a = 10,
             "window_size"_a = 7, nb::call_guard<nb::gil_scoped_release>(),
             "Find approximately `target` best-scoring hits per matrix.");
    scan.def("naive_scan_dna",
             static_cast<std::vector<match> (*)(const std::string &, score_matrix,
                                                double)>(&S::naive_scan_dna),
             "seq"_a, "matrix"_a, "threshold"_a,
             nb::call_guard<nb::gil_scoped_release>(),
             "Brute-force scan; reference implementation for testing.");
    scan.def("naive_scan_dna",
             static_cast<std::vector<match> (*)(const std::string &, score_matrix,
                                                double, size_t)>(
                 &S::naive_scan_dna),
             "seq"_a, "matrix"_a, "threshold"_a, "a"_a,
             nb::call_guard<nb::gil_scoped_release>(),
             "Brute-force scan for high-order matrices.");

    nb::class_<Scanner>(
        scan, "Scanner",
        "Persistent scanner: pays the matrix preprocessing cost once and can\n"
        "then scan many sequences.")
        .def(nb::init<unsigned int>(), "window_size"_a)
        .def(nb::init<unsigned int, const std::vector<std::string> &>(),
             "window_size"_a, "alphabet"_a)
        .def("set_motifs", &Scanner::set_motifs, "matrices"_a, "bg"_a,
             "thresholds"_a, nb::call_guard<nb::gil_scoped_release>(),
             "Preprocess matrices for scanning. `bg` only affects search\n"
             "optimisation, not the reported scores.")
        .def("scan", &Scanner::scan, "seq"_a,
             nb::call_guard<nb::gil_scoped_release>(),
             "All hits above threshold, per matrix.")
        .def("scan_max_hits", &Scanner::scan_max_hits, "seq"_a, "max_hits"_a,
             nb::call_guard<nb::gil_scoped_release>(),
             "Hits per matrix, giving up once `max_hits` is exceeded.")
        .def("counts_max_hits", &Scanner::counts_max_hits, "seq"_a, "max_hits"_a,
             nb::call_guard<nb::gil_scoped_release>(),
             "Hit counts per matrix, giving up once `max_hits` is exceeded.")
        .def("variant_matches", &Scanner::variant_matches, "seq"_a, "variants"_a,
             "max_depth"_a = 0, nb::call_guard<nb::gil_scoped_release>(),
             "Hits that only exist once one or more variants are applied.")
        .def("size", &Scanner::size, "Number of motifs currently loaded.")
        .def("__len__", &Scanner::size);

    // ------------------------------------------------------------------
    // misc
    // ------------------------------------------------------------------
    nb::module_ misc = m.def_submodule("misc", "Bit-twiddling helpers.");
    namespace M = MOODS::misc;

    misc.def("shift", &M::shift, "a"_a, "Bits needed per symbol (ceil(log2 a)).");
    misc.def("mask", &M::mask, "a"_a, "Bit mask covering one symbol.");
    misc.def("q_gram_size", &M::q_gram_size, "rows"_a, "a"_a,
             "q-gram length implied by a matrix with `rows` rows.");
    misc.def("rc_tuple", &M::rc_tuple, "CODE"_a, "a"_a, "q"_a,
             "Reverse complement of an encoded q-gram.");
    misc.def("preprocess_seq", &M::preprocess_seq, "s"_a, "a"_a,
             "alphabet_map"_a, nb::call_guard<nb::gil_scoped_release>(),
             "Bounds of the maximal scannable (valid-symbol) intervals.");
}
