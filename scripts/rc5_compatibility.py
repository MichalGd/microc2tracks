"""Exact RC5 migration allowlist; payload/provenance checks are never waived."""
import json
from pathlib import Path

from workflow_common import json_write, truth

RC5_HASHES = {
    "workflow_common.py": "085d1aef6e9b64f54d7b8cd74bcc9574921c5e37ba4532d6883753fd5ae31a46",
    "pairs_processing.py": "157390df0b8e1473cbf0137e9a2a2c1a714363e183c120dd10d41fa4881fb700",
    "matrix_products.py": "5999ef12a15cb74978885539e33224f05c787776867e7285210e316f6dfa3acc",
    "pipeline.py": "711cf1b248e0e6891fcd86914c0aff02147cdb17d739ec813809f7e7804e118c",
    "filter_contact_pairs.py": "b6579afbeb980edce24943ec3a1bf49d979ab4b9bd6bce959f4857bf033598e8",
    "prepare_restriction_fragments.py": "f898079ff6185749ec2be0af60709d15cfaa04c8bf34e88edf9fb5620d7b1cee",
    "refilter_inputs.py": "971ca452af12b9859559d35b86f2604eeb855ada0446bdb149df897c5d4026bb",
    "run_downstream.sh": "907037f4790061e84cffa3a52c0ee50ff435d564dbf470fc57433be0056923b2",
    "run_plan.py": "7d1bc356e953895713cac2dd04826937d57563c7849a6fb7c9bf381d10d2f8ae",
    "artifact_state.py": "c35d84d10845cf48fff0885c8d4d60168dd052500ec4a22f7c4038de6e5500eb",
}


def enable(plan, preserve=True):
    if not truth(plan["settings"]["REUSE_RC5_OUTPUTS"]):
        return {}
    metadata = Path(plan.get("metadata_dir", str(Path(plan["roots"]["out"]) / "run_metadata")))
    for name in ("source_hashes.json", "rc5_source_hashes.json"):
        path = metadata / name
        if path.is_file():
            saved = json.loads(path.read_text())
            if all(saved.get(k) == v for k, v in RC5_HASHES.items()):
                destination = metadata / "rc5_source_hashes.json"
                if preserve and not destination.exists():
                    json_write(destination, saved)
                return dict(RC5_HASHES)
    return {}
