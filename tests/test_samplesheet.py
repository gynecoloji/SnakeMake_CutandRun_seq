import sys
import pytest

sys.path.insert(0, "workflow/scripts")
import samplesheet as ss


def _sheet(tmp_path, rows):
    p = tmp_path / "samples.csv"
    p.write_text("sample_id,condition,replicate,input_control,igg_control,peak_mode,notes\n" + rows)
    return ss.SampleSheet(str(p))


REF = ("t1,cJUN,1,igg,,narrow,a\n"
       "t2,cJUN,2,igg,,narrow,b\n"
       "t3,K27,1,igg,,broad,c\n"
       "igg,IgG,1,,,,ctrl\n")


def test_sets(tmp_path):
    s = _sheet(tmp_path, REF)
    assert s.all_samples == ["t1", "t2", "t3", "igg"]
    assert s.treatment_samples == ["t1", "t2", "t3"]
    assert s.control_samples == ["igg"]
    assert s.narrow_samples == ["t1", "t2"]
    assert s.broad_samples == ["t3"]


def test_groups_and_method(tmp_path):
    s = _sheet(tmp_path, REF)
    assert s.groups == {"cJUN": ["t1", "t2"], "K27": ["t3"]}
    assert s.group_method == {"cJUN": "idr", "K27": "single"}
    assert s.idr_groups == ["cJUN"]
    assert s.idr_samples == ["t1", "t2"]
    assert s.idr_pairs == [("cJUN", "t1", "t2")]


def test_helpers(tmp_path):
    s = _sheet(tmp_path, REF)
    assert s.macs2_ext("t1") == "narrowPeak"
    assert s.macs2_ext("t3") == "broadPeak"
    assert s.input_control("t1") == "igg"
    assert s.input_control("igg") == ""
    cfg = {"seacr_narrow_stringency": "stringent", "seacr_broad_stringency": "relaxed"}
    assert s.seacr_stringency("t1", cfg) == "stringent"
    assert s.seacr_stringency("t3", cfg) == "relaxed"


def test_majority_method(tmp_path):
    s = _sheet(tmp_path, "a,X,1,igg,,narrow,\nb,X,2,igg,,narrow,\nc,X,3,igg,,narrow,\nigg,IgG,1,,,,\n")
    assert s.group_method["X"] == "majority"


def test_validate_bad_control_ref(tmp_path):
    s = _sheet(tmp_path, "t1,cJUN,1,missing,,narrow,\nigg,IgG,1,,,,\n")
    with pytest.raises(ValueError, match="input_control"):
        s.validate()


def test_validate_mixed_peakmode(tmp_path):
    s = _sheet(tmp_path, "t1,cJUN,1,igg,,narrow,\nt2,cJUN,2,igg,,broad,\nigg,IgG,1,,,,\n")
    with pytest.raises(ValueError, match="peak_mode"):
        s.validate()


def test_validate_bad_peakmode_value(tmp_path):
    s = _sheet(tmp_path, "t1,cJUN,1,igg,,wide,\nigg,IgG,1,,,,\n")
    with pytest.raises(ValueError, match="narrow|broad"):
        s.validate()


def test_validate_ok(tmp_path):
    _sheet(tmp_path, REF).validate()  # no raise


DUAL = ("t1,cJUN,1,in1,igg1,narrow,a\n"
        "t2,cJUN,2,,igg1,narrow,b\n"
        "in1,Input,1,,,,c\n"
        "igg1,IgG,1,,,,d\n")

def _dual(tmp_path, rows):
    p = tmp_path / "s.csv"
    p.write_text("sample_id,condition,replicate,input_control,igg_control,peak_mode,notes\n" + rows)
    return ss.SampleSheet(str(p))

def test_resolved_control_igg_primary(tmp_path):
    s = _dual(tmp_path, DUAL)
    assert s.resolved_control("t1", "igg") == "igg1"   # igg primary
    assert s.resolved_control("t1", "input") == "in1"  # input primary

def test_resolved_control_fallback(tmp_path):
    s = _dual(tmp_path, DUAL)
    # t2 has only igg1; input requested -> falls back to igg1
    assert s.resolved_control("t2", "input") == "igg1"
    assert s.resolved_control("t2", "igg") == "igg1"

def test_resolved_control_none(tmp_path):
    s = _dual(tmp_path, "t1,cJUN,1,,,narrow,a\nc1,Ctl,1,,,,x\n")
    assert s.resolved_control("t1", "igg") == ""

def test_validate_bad_igg_ref(tmp_path):
    s = _dual(tmp_path, "t1,cJUN,1,,missing,narrow,a\nigg1,IgG,1,,,,d\n")
    with pytest.raises(ValueError, match="igg_control|control"):
        s.validate()
