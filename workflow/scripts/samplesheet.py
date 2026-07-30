#!/usr/bin/env python3
"""Sample-sheet parsing, derived sample sets, validation and helpers for the
CUT&RUN workflow. Pure module (no snakemake import) so it is unit-testable and
importable from common.smk."""
from itertools import combinations
import pandas as pd

REQUIRED_COLUMNS = ["sample_id", "condition", "replicate", "input_control",
                    "igg_control", "peak_mode", "notes"]
VALID_PEAK_MODES = {"narrow", "broad"}


def load_samples(path):
    df = pd.read_csv(path, dtype=str).fillna("")
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"samples sheet {path} missing columns: {missing}")
    for c in ("sample_id", "condition", "input_control", "igg_control", "peak_mode"):
        df[c] = df[c].astype(str).str.strip()
    df["replicate"] = df["replicate"].astype(str).str.strip()
    return df


def _repro_method(n):
    if n >= 3:
        return "majority"
    if n == 2:
        return "idr"
    return "single"


class SampleSheet:
    def __init__(self, path):
        self.df = load_samples(path)
        d = self.df
        self.all_samples = d["sample_id"].tolist()
        self._peak_mode = dict(zip(d["sample_id"], d["peak_mode"]))
        self._input = dict(zip(d["sample_id"], d["input_control"]))
        self._igg = dict(zip(d["sample_id"], d["igg_control"]))
        self.treatment_samples = [s for s in self.all_samples if self._peak_mode[s]]
        self.control_samples = [s for s in self.all_samples if not self._peak_mode[s]]
        self.narrow_samples = [s for s in self.treatment_samples if self._peak_mode[s] == "narrow"]
        self.broad_samples = [s for s in self.treatment_samples if self._peak_mode[s] == "broad"]
        trt = d[d["sample_id"].isin(self.treatment_samples)]
        self.groups = {g: m["sample_id"].tolist()
                       for g, m in trt.groupby("condition", sort=False)}
        self.group_method = {g: _repro_method(len(m)) for g, m in self.groups.items()}
        self.idr_groups = [g for g in self.groups if self.group_method[g] == "idr"]
        self.idr_samples = [s for g in self.idr_groups for s in self.groups[g]]
        self.idr_pairs = []
        for g, members in self.groups.items():
            for a, b in combinations(members, 2):
                self.idr_pairs.append((g, a, b))

    def peak_mode(self, sample):
        return self._peak_mode.get(sample, "")

    def macs2_ext(self, sample):
        return "broadPeak" if self._peak_mode.get(sample) == "broad" else "narrowPeak"

    def input_control(self, sample):
        return self._input.get(sample, "")

    def igg_control(self, sample):
        return self._igg.get(sample, "")

    def resolved_control(self, sample, control_type):
        """Effective MACS2/SEACR/track control for an IP sample: the column named
        by control_type (input|igg), else the other column, else '' (none)."""
        inp, igg = self._input.get(sample, ""), self._igg.get(sample, "")
        primary, secondary = (inp, igg) if control_type == "input" else (igg, inp)
        return primary or secondary or ""

    def seacr_stringency(self, sample, cfg):
        if self._peak_mode.get(sample) == "broad":
            return cfg["seacr_broad_stringency"]
        return cfg["seacr_narrow_stringency"]

    def validate(self):
        controls = set(self.control_samples)
        for s in self.treatment_samples:
            for col, val in (("input_control", self._input[s]), ("igg_control", self._igg[s])):
                if val and val not in controls:
                    raise ValueError(
                        f"{col} '{val}' for sample '{s}' is not an existing "
                        f"control (empty peak_mode) sample")
        for g, members in self.groups.items():
            modes = {self._peak_mode[s] for s in members}
            if len(modes) > 1:
                raise ValueError(
                    f"condition '{g}' mixes peak_mode values {modes}; all "
                    f"replicates of a condition must share narrow or broad")
        for s in self.all_samples:
            pm = self._peak_mode[s]
            if pm and pm not in VALID_PEAK_MODES:
                raise ValueError(
                    f"sample '{s}' has peak_mode '{pm}'; must be narrow, broad or empty")
