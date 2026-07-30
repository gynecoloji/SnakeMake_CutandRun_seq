import sys

sys.path.insert(0, "workflow/scripts")
import consensus_peaks as cp


def _write(p, lines):
    p.write_text("\n".join("\t".join(map(str, r)) for r in lines) + "\n")


def test_load_narrowpeak_uses_summit_offset(tmp_path):
    f = tmp_path / "a_peaks.narrowPeak"
    # chrom start end name score strand signal p q summitoffset
    _write(f, [["chr1", 100, 300, "p1", 50, ".", 5.0, 9.0, 7.5, 40]])
    peaks = cp.load_peaks(str(f), "a")
    assert peaks[0].summit == 140          # start + offset
    assert peaks[0].score == 7.5           # col 9 (-log10 q)


def test_load_broadpeak_uses_midpoint(tmp_path):
    f = tmp_path / "b_peaks.broadPeak"
    # broadPeak has 9 columns, no summit offset
    _write(f, [["chr1", 100, 300, "p1", 50, ".", 5.0, 9.0, 7.5]])
    peaks = cp.load_peaks(str(f), "b")
    assert peaks[0].summit == 200          # midpoint (100+300)//2
    assert peaks[0].score == 7.5


# ── mode-aware consensus geometry ─────────────────────────────────────────
KEEP = r"^chr([1-9]|1[0-9]|2[0-2]|X)$"


def _bl(tmp_path):
    p = tmp_path / "bl.bed"
    p.write_text("")            # empty blacklist
    return str(p)


def _consensus(tmp_path, fname, rows, mode, method="single"):
    f = tmp_path / fname
    _write(f, rows)
    return cp.build_consensus(
        groups={"G": ["s1"]}, group_method={"G": method},
        narrowpeak_paths={"s1": str(f)}, idr_paths={},
        min_reps=2, width=500, keep_regex=KEEP, blacklist_path=_bl(tmp_path),
        group_mode={"G": mode})


def test_broad_group_keeps_variable_width(tmp_path):
    # two non-overlapping broad domains of different widths -> widths PRESERVED
    cons = _consensus(tmp_path, "b_peaks.broadPeak", [
        ["chr1", 1000, 6000, "p1", 50, ".", 5.0, 9.0, 7.5],    # 5000 bp
        ["chr1", 20000, 20800, "p2", 40, ".", 4.0, 8.0, 6.0],  # 800 bp
    ], mode="broad")
    assert sorted(w["end"] - w["start"] for w in cons) == [800, 5000]  # not 500/500


def test_broad_group_merges_overlapping(tmp_path):
    # overlapping broad domains -> merged into one variable-width region
    cons = _consensus(tmp_path, "b_peaks.broadPeak", [
        ["chr1", 1000, 3000, "p1", 50, ".", 5.0, 9.0, 7.5],
        ["chr1", 2500, 4000, "p2", 40, ".", 4.0, 8.0, 6.0],   # overlaps p1
    ], mode="broad")
    assert len(cons) == 1
    assert (cons[0]["start"], cons[0]["end"]) == (1000, 4000)


def test_narrow_group_still_fixed_width(tmp_path):
    # narrow peak -> fixed 500 bp window centered on the summit (unchanged)
    cons = _consensus(tmp_path, "n_peaks.narrowPeak", [
        ["chr1", 1000, 1300, "p1", 50, ".", 5.0, 9.0, 7.5, 100],  # summit = 1100
    ], mode="narrow")
    assert len(cons) == 1
    assert cons[0]["end"] - cons[0]["start"] == 500
    assert cons[0]["start"] == 1100 - 250          # summit-centered


def test_mixed_narrow_fixed_and_broad_variable(tmp_path):
    # one narrow group (fixed) + one broad group (variable) in one consensus
    nf = tmp_path / "n_peaks.narrowPeak"
    _write(nf, [["chr1", 1000, 1300, "p", 50, ".", 5.0, 9.0, 7.5, 100]])   # summit 1100
    bf = tmp_path / "b_peaks.broadPeak"
    _write(bf, [["chr2", 5000, 12000, "q", 40, ".", 4.0, 8.0, 6.0]])       # 7000 bp
    cons = cp.build_consensus(
        groups={"N": ["n1"], "B": ["b1"]},
        group_method={"N": "single", "B": "single"},
        narrowpeak_paths={"n1": str(nf), "b1": str(bf)}, idr_paths={},
        min_reps=2, width=500, keep_regex=KEEP, blacklist_path=_bl(tmp_path),
        group_mode={"N": "narrow", "B": "broad"})
    widths = {w["chrom"]: w["end"] - w["start"] for w in cons}
    assert widths["chr1"] == 500      # narrow -> fixed
    assert widths["chr2"] == 7000     # broad -> variable
