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
