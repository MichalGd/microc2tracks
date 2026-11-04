"""Narrow RC4 scientific-contract compatibility; never waive payload validation.

RC5 changes scheduling, resource arguments and output selection, not the RC4
scientific transforms. Only an exact known RC4 source fingerprint may supply
alternative completion signatures. Original artifacts, proofs and records stay
unchanged. Future scientific changes must remove/revise this allowlist.
"""
from pathlib import Path
import json

from workflow_common import VERSION, json_write, truth

RC4_HASHES = {
    "workflow_common.py": "9c19b61d9eb31c34048363ed93e2d6006710781cdb34fec7384e007325a6b497",
    "pairs_processing.py": "6b55a5ee5c5cff2e60b3c51ef7d437bd29173750a7a6ec4ee5b5eed9cd5d6e07",
    "matrix_products.py": "18e2301b88a2feb61e6ae3e2f2ad7fb5c797da6dcf9d0549fa7db20b5f581a0e",
    "pipeline.py": "a7dd094fa0af85cb278e5bea66c3d85b769efefaa4f1b047f6715e22ba27903b",
    "filter_contact_pairs.py": "b6579afbeb980edce24943ec3a1bf49d979ab4b9bd6bce959f4857bf033598e8",
    "prepare_restriction_fragments.py": "f898079ff6185749ec2be0af60709d15cfaa04c8bf34e88edf9fb5620d7b1cee",
    "refilter_inputs.py": "971ca452af12b9859559d35b86f2604eeb855ada0446bdb149df897c5d4026bb",
    "run_downstream.sh": "907037f4790061e84cffa3a52c0ee50ff435d564dbf470fc57433be0056923b2",
    "run_plan.py": "06562ee84829f696e6629adb0722eed3879d8e00510d4618f2cfde91358b2436",
    "artifact_state.py": "5118fe262d6b340e7aa22d0d76baa85cbcccf85df93b95ee497e1792ffe26e88",
}


def enable(plan, current_hashes, preserve=True):
    if VERSION not in {"3.0.0-rc5", "3.0.0-rc6"} or not truth(plan["settings"]["REUSE_RC4_OUTPUTS"]):
        return {}
    metadata = Path(plan.get("metadata_dir", str(Path(plan["roots"]["out"]) / "run_metadata")))
    for filename in ("source_hashes.json", "rc4_source_hashes.json"):
        path = metadata / filename
        if not path.is_file():
            continue
        saved = json.loads(path.read_text())
        if all(saved.get(key) == value for key, value in RC4_HASHES.items()):
            immutable = ("prepare_restriction_fragments.py", "refilter_inputs.py", "run_downstream.sh")
            if any(current_hashes.get(name) != RC4_HASHES[name] for name in immutable):
                raise RuntimeError("RC4 compatibility is not certified for modified scientific modules; disable REUSE_RC4_OUTPUTS and review recomputation")
            preserved = metadata / "rc4_source_hashes.json"
            if preserve and not preserved.exists():
                json_write(preserved, saved)
            return dict(RC4_HASHES)
    previous = metadata / "run_plan.json"
    if previous.is_file() and json.loads(previous.read_text()).get("workflow_version") == "3.0.0-rc4":
        raise RuntimeError("RC4 source fingerprint is missing or unrecognized; reuse cannot be certified. Preserve this tree and use a new OUTDIR or review its provenance.")
    return {}


def translate(value, alternatives):
    if isinstance(value, dict):
        return {key: translate(item, alternatives) for key, item in value.items()}
    if isinstance(value, list):
        return [translate(item, alternatives) for item in value]
    if isinstance(value, tuple):
        return tuple(translate(item, alternatives) for item in value)
    return alternatives.get(value, [value])[0] if isinstance(value, str) else value
