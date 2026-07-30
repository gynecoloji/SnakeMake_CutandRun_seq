#!/usr/bin/env python3
"""Reproducible, MODE-AWARE non-overlapping consensus peak set.

Replicate-count-adaptive reproducibility selects the surviving peaks per group:
  * >=3 reps -> majority vote (summit covered by >= min_reps replicates)
  * ==2 reps -> IDR peaks (precomputed by the reproducible_idr rule)
  * ==1 rep  -> passthrough (flagged; not reproducibility-filtered)
The consensus GEOMETRY then follows the group's peak_mode:
  * narrow -> Corces et al. 2018 fixed-width summit windows + SPM-ranked iterative
              overlap removal (point-source model; MACS2 gives a real summit).
  * broad  -> the surviving peaks' ORIGINAL variable-width intervals, merged into a
              non-overlapping union (broad marks are domains; a fixed window around
              the midpoint would discard their extent).
An all-narrow run is byte-identical to the fixed-width-only behavior.
"""
import re
import bisect
from pathlib import Path
from collections import defaultdict


class Peak:
    __slots__ = ("chrom", "start", "end", "score", "summit", "sample", "spm")

    def __init__(self, chrom, start, end, score, summit, sample):
        self.chrom = chrom
        self.start = start      # original peak start (0-based)
        self.end = end          # original peak end
        self.score = score      # MACS2 -log10(q) (narrowPeak col 9)
        self.summit = summit    # absolute summit = start + col10 offset
        self.sample = sample    # sample id (or group name for IDR peaks)
        self.spm = 0.0


def load_peaks(path, sample):
    """Parse a MACS2 narrowPeak (10-col) or broadPeak (9-col) file into Peak
    objects tagged with `sample`. narrowPeak: summit = start + col10 offset;
    broadPeak: summit = peak midpoint. score = col9 (-log10 q) for both."""
    peaks = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        start, end = int(f[1]), int(f[2])
        score = float(f[8])
        summit = start + int(f[9]) if len(f) >= 10 else (start + end) // 2
        peaks.append(Peak(f[0], start, end, score, summit, sample))
    return peaks


def _point_covered(pos, intervals):
    """True if pos falls in any (start, end) of the sorted intervals list."""
    for s, e in intervals:
        if s <= pos < e:
            return True
        if s > pos:
            break
    return False


def majority_keep(rep_peaks, min_reps):
    """rep_peaks: list (per replicate) of Peak lists. Keep peaks whose summit is
    covered by peaks in >= min_reps replicates of the same condition."""
    rep_index = []
    for peaks in rep_peaks:
        by_chrom = defaultdict(list)
        for p in peaks:
            by_chrom[p.chrom].append((p.start, p.end))
        for c in by_chrom:
            by_chrom[c].sort()
        rep_index.append(by_chrom)

    kept = []
    for peaks in rep_peaks:
        for p in peaks:
            cover = sum(1 for idx in rep_index if _point_covered(p.summit, idx.get(p.chrom, [])))
            if cover >= min_reps:
                kept.append(p)
    return kept


def assign_spm(peaks):
    """Assign per-sample score-per-million to each peak (in place)."""
    total = defaultdict(float)
    for p in peaks:
        total[p.sample] += p.score
    for p in peaks:
        denom = total[p.sample] / 1e6
        p.spm = p.score / denom if denom > 0 else 0.0
    return peaks


def fixed_window(peak, width):
    """Return (chrom, wstart, wend): fixed-`width` window centered on the summit."""
    wstart = peak.summit - width // 2
    if wstart < 0:
        wstart = 0
    return (peak.chrom, wstart, wstart + width)


def _load_bed_intervals(path):
    by_chrom = defaultdict(list)
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        by_chrom[f[0]].append((int(f[1]), int(f[2])))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom


def _overlaps_any(chrom, start, end, by_chrom):
    for s, e in by_chrom.get(chrom, []):
        if s < end and start < e:
            return True
        if s >= end:
            break
    return False


def iterative_overlap_removal(windows, width):
    """Greedily keep highest-SPM window; drop later windows overlapping a kept one.
    All windows are exactly `width` bp, so two overlap iff |start_a - start_b| < width."""
    windows = sorted(windows, key=lambda w: w["spm"], reverse=True)
    kept = []
    starts_by_chrom = defaultdict(list)
    for w in windows:
        starts = starts_by_chrom[w["chrom"]]
        i = bisect.bisect_left(starts, w["start"])
        overlap = (i > 0 and w["start"] - starts[i - 1] < width) or \
                  (i < len(starts) and starts[i] - w["start"] < width)
        if not overlap:
            bisect.insort(starts, w["start"])
            kept.append(w)
    return kept


def _merge_spm_intervals(intervals, touching):
    """Sort {chrom,start,end,spm} and merge into a non-overlapping set carrying the
    max spm. touching=True merges book-ended intervals (s <= prev_end, natural for
    domain unions); touching=False merges only true overlaps (s < prev_end, so
    adjacent fixed-width windows are left untouched)."""
    out = []
    for w in sorted(intervals, key=lambda w: (w["chrom"], w["start"], w["end"])):
        prev = out[-1] if out else None
        joins = prev is not None and prev["chrom"] == w["chrom"] and (
            w["start"] <= prev["end"] if touching else w["start"] < prev["end"])
        if joins:
            prev["end"] = max(prev["end"], w["end"])
            prev["spm"] = max(prev["spm"], w["spm"])
        else:
            out.append(dict(w))
    return out


def _reproducible_peaks(groups, group_method, narrowpeak_paths, idr_paths, min_reps):
    """Surviving peaks per group after the reproducibility filter (majority / IDR /
    single passthrough). Returns the flat list of kept Peak objects."""
    kept = []
    for g, members in groups.items():
        method = group_method[g]
        if method == "majority":
            reps = [load_peaks(narrowpeak_paths[s], s) for s in members]
            kept.extend(majority_keep(reps, min_reps))
        elif method == "single":
            kept.extend(load_peaks(narrowpeak_paths[members[0]], members[0]))
        elif method == "idr":
            kept.extend(load_peaks(idr_paths[g], g))
    return kept


def build_consensus(groups, group_method, narrowpeak_paths, idr_paths,
                    min_reps, width, keep_regex, blacklist_path, group_mode):
    """Mode-aware reproducible consensus -> sorted, named list of
    {chrom,start,end,name,spm}. Reproducibility (majority / IDR / single) selects the
    surviving peaks per group; the GEOMETRY then depends on the group's peak_mode:

      * narrow groups -> fixed-`width` summit windows + SPM iterative overlap removal
        (Corces 2018); the point-source model where MACS2 gives a real summit.
      * broad  groups -> the surviving peaks' ORIGINAL variable-width intervals,
        merged into a non-overlapping union (broad marks are domains; a fixed window
        around the midpoint would discard their extent).

    An all-narrow run is byte-identical to the fixed-width-only behavior. When both
    modes are present the two sets are combined and any residual overlaps merged
    (a narrow window overlapping a broad domain fuses into the domain)."""
    keep_re = re.compile(keep_regex)
    blacklist = _load_bed_intervals(blacklist_path)
    group_mode = group_mode or {}

    kept = _reproducible_peaks(groups, group_method, narrowpeak_paths, idr_paths, min_reps)
    # Partition survivors by their group's mode. A peak is tagged with either its
    # sample id (majority/single) or its group name (IDR peaks), so map both.
    sample_mode = {s: group_mode.get(g, "narrow") for g, members in groups.items() for s in members}
    sample_mode.update({g: group_mode.get(g, "narrow") for g in groups})       # IDR peaks tagged by group
    narrow_kept = [p for p in kept if sample_mode.get(p.sample, "narrow") != "broad"]
    broad_kept  = [p for p in kept if sample_mode.get(p.sample, "narrow") == "broad"]

    assign_spm(narrow_kept + broad_kept)   # per-sample SPM over all survivors

    # Narrow: fixed-width summit windows + SPM iterative overlap removal.
    narrow_windows = []
    for p in narrow_kept:
        if not keep_re.fullmatch(p.chrom):
            continue
        chrom, ws, we = fixed_window(p, width)
        if _overlaps_any(chrom, ws, we, blacklist):
            continue
        narrow_windows.append({"chrom": chrom, "start": ws, "end": we, "spm": p.spm})
    narrow_consensus = iterative_overlap_removal(narrow_windows, width)

    if not broad_kept:                     # all-narrow: preserve the legacy output exactly
        consensus = narrow_consensus
    else:
        # Broad: keep original variable-width intervals, then merge overlapping.
        broad_intervals = []
        for p in broad_kept:
            if not keep_re.fullmatch(p.chrom):
                continue
            if _overlaps_any(p.chrom, p.start, p.end, blacklist):
                continue
            broad_intervals.append({"chrom": p.chrom, "start": p.start, "end": p.end, "spm": p.spm})
        broad_consensus = _merge_spm_intervals(broad_intervals, touching=True)
        # Combine and resolve any narrow/broad overlaps (true overlaps only).
        consensus = _merge_spm_intervals(narrow_consensus + broad_consensus, touching=False)

    consensus.sort(key=lambda w: (w["chrom"], w["start"]))
    for i, w in enumerate(consensus, 1):
        w["name"] = f"consensus_peak_{i}"
    return consensus


def write_bed(consensus, path):
    lines = [f'{w["chrom"]}\t{w["start"]}\t{w["end"]}\t{w["name"]}\t{int(round(w["spm"]))}\t.'
             for w in consensus]
    Path(path).write_text("\n".join(lines) + ("\n" if lines else ""))


def write_saf(consensus, path):
    lines = ["GeneID\tChr\tStart\tEnd\tStrand"]
    for w in consensus:  # SAF is 1-based, inclusive
        lines.append(f'{w["name"]}\t{w["chrom"]}\t{w["start"] + 1}\t{w["end"]}\t.')
    Path(path).write_text("\n".join(lines) + "\n")


if "snakemake" in globals():  # pragma: no cover
    sm = snakemake  # noqa: F821
    _groups = dict(sm.params.groups)
    _method = dict(sm.params.group_method)
    # Per-sample MACS2 peak paths carry the correct narrow/broad extension.
    _npaths = dict(sm.params.narrowpeak_paths)
    _ipaths = dict(sm.params.idr_paths)
    _gmode = dict(sm.params.group_mode)
    _consensus = build_consensus(_groups, _method, _npaths, _ipaths,
                                 int(sm.params.min_reps), int(sm.params.window),
                                 sm.params.keep_regex, str(sm.input.blacklist), _gmode)
    write_bed(_consensus, str(sm.output.bed))
    write_saf(_consensus, str(sm.output.saf))
