import sys

sys.path.insert(0, "workflow/scripts")
import seacr_consensus as sc


def test_seacr_bed_load(tmp_path):
    f = tmp_path / "s.bed"
    f.write_text("chr1\t100\t300\t12.5\t9.0\tchr1:180-200\n")
    ivs = sc.load_seacr_bed(str(f))
    assert ivs == [("chr1", 100, 300)]


def test_reproducible_union_majority():
    # condition X, 3 reps, min_reps=2: peak present in >=2 reps kept
    groups = {"X": ["a", "b", "c"]}
    method = {"X": "majority"}
    peaks = {
        "a": [("chr1", 100, 200)],
        "b": [("chr1", 150, 250)],           # overlaps a
        "c": [("chr2", 500, 600)],           # alone
    }
    out = sc.reproducible_union(groups, method, peaks, min_reps=2)
    # chr1 region reproducible (a,b overlap); chr2 not (only c)
    assert any(c == "chr1" for c, s, e in out)
    assert not any(c == "chr2" for c, s, e in out)


def test_reproducible_union_single_passthrough():
    groups = {"Y": ["a"]}
    method = {"Y": "single"}
    peaks = {"a": [("chr3", 10, 20)]}
    out = sc.reproducible_union(groups, method, peaks, min_reps=2)
    assert ("chr3", 10, 20) in out


def test_merge_and_filter():
    ivs = [("chr1", 100, 200), ("chr1", 150, 300), ("chrM", 0, 50)]
    merged = sc.merge_and_filter(ivs, keep_regex=r"^chr([1-9]|1[0-9]|2[0-2]|X)$",
                                 blacklist={})
    assert ("chr1", 100, 300) in merged
    assert not any(c == "chrM" for c, s, e in merged)
