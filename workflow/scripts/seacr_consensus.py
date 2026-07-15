#!/usr/bin/env python3
"""Overlap-based reproducible variable-width union of SEACR peaks → BED + SAF.

SEACR emits enrichment domains (no summit / p-value), so reproducibility is
overlap-based, not IDR/fixed-width. Per condition: keep a peak if it overlaps
peaks in >=K replicates (K = min_reps for >=3 reps; 2-of-2 for 2 reps;
passthrough for 1 rep). Union across conditions, merge overlapping, drop
blacklist / off-target chroms.
"""
import re
from pathlib import Path
from collections import defaultdict


def load_seacr_bed(path):
    """SEACR .bed: chrom start end total_signal max_signal max_region. Return
    [(chrom, start, end)]."""
    out = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        out.append((f[0], int(f[1]), int(f[2])))
    return out


def _overlaps(a, intervals):
    c, s, e = a
    for c2, s2, e2 in intervals:
        if c2 == c and s2 < e and s < e2:
            return True
    return False


def reproducible_union(groups, method, peaks_by_sample, min_reps):
    """Return list of (chrom,start,end) reproducible peaks unioned across
    conditions (not yet merged)."""
    kept = []
    for g, members in groups.items():
        m = method[g]
        rep_lists = [peaks_by_sample.get(s, []) for s in members]
        if m == "single":
            for s in members:
                kept.extend(peaks_by_sample.get(s, []))
            continue
        k = min_reps if m == "majority" else 2
        for i, s in enumerate(members):
            others = [rep_lists[j] for j in range(len(rep_lists)) if j != i]
            for p in peaks_by_sample.get(s, []):
                cover = 1 + sum(1 for o in others if _overlaps(p, o))
                if cover >= k:
                    kept.append(p)
    return kept


def merge_and_filter(intervals, keep_regex, blacklist):
    """Sort, chrom/blacklist filter, then merge overlapping into non-overlapping
    (chrom,start,end)."""
    keep_re = re.compile(keep_regex)
    bl = blacklist or {}

    def _bl_hit(c, s, e):
        for bs, be in bl.get(c, []):
            if bs < e and s < be:
                return True
        return False

    ivs = sorted((c, s, e) for (c, s, e) in intervals
                 if keep_re.fullmatch(c) and not _bl_hit(c, s, e))
    merged = []
    for c, s, e in ivs:
        if merged and merged[-1][0] == c and s <= merged[-1][2]:
            pc, ps, pe = merged[-1]
            merged[-1] = (pc, ps, max(pe, e))
        else:
            merged.append((c, s, e))
    return merged


def _load_blacklist(path):
    by = defaultdict(list)
    if not path:
        return by
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        by[f[0]].append((int(f[1]), int(f[2])))
    return by


def write_bed(consensus, path):
    lines = [f"{c}\t{s}\t{e}\tseacr_consensus_peak_{i}\t.\t."
             for i, (c, s, e) in enumerate(consensus, 1)]
    Path(path).write_text("\n".join(lines) + ("\n" if lines else ""))


def write_saf(consensus, path):
    lines = ["GeneID\tChr\tStart\tEnd\tStrand"]
    for i, (c, s, e) in enumerate(consensus, 1):
        lines.append(f"seacr_consensus_peak_{i}\t{c}\t{s + 1}\t{e}\t.")
    Path(path).write_text("\n".join(lines) + "\n")


if "snakemake" in globals():  # pragma: no cover
    sm = snakemake  # noqa: F821
    _groups = dict(sm.params.groups)
    _method = dict(sm.params.group_method)
    _paths = dict(sm.params.seacr_paths)  # sample -> bed path
    _peaks = {s: load_seacr_bed(p) for s, p in _paths.items()}
    _union = reproducible_union(_groups, _method, _peaks, int(sm.params.min_reps))
    _blacklist = _load_blacklist(str(sm.input.blacklist))
    _consensus = merge_and_filter(_union, sm.params.keep_regex, _blacklist)
    write_bed(_consensus, str(sm.output.bed))
    write_saf(_consensus, str(sm.output.saf))
