"""One matrix implementation for all sample, joint-unit and pooled accepted pairs."""
from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

from workflow_common import builds_hic, json_write, text_open, truth


def export_juicer(source: Path, destination: Path) -> int:
    """Project by column name to Juicer short format, with distinct dummy fragments.

    Positions become zero-based: Juicer bins positions by integer division,
    while cooler's pairs reader converts the input's one-based positions.
    No contact selection occurs here (including same-fragment contacts).
    """
    count, columns = 0, None
    with text_open(source) as handle, gzip.open(destination, "wt", encoding="utf-8", newline="\n") as out:
        for line in handle:
            if line.startswith("#columns:"):
                columns = line.split(":", 1)[1].split()
            elif not line.startswith("#") and line.strip():
                if columns is None:
                    raise ValueError("Juicer export requires a #columns schema")
                row = dict(zip(columns, line.rstrip().split("\t")))
                c1, c2, p1, p2 = row["chrom1"], row["chrom2"], int(row["pos1"])-1, int(row["pos2"])-1
                if min(p1, p2) < 0:
                    raise ValueError("Accepted pairs must have positive one-based positions")
                s1, s2 = ("0" if row[f"strand{side}"] == "+" else "16" for side in (1, 2))
                out.write(f"{s1}\t{c1}\t{p1}\t0\t{s2}\t{c2}\t{p2}\t1\n")
                count += 1
    if columns is None:
        raise ValueError("Missing #columns in Juicer export")
    return count


class MatrixProducts:
    def __init__(self, runtime):
        self.rt, self.c = runtime, runtime.c

    def build(self, entity, accepted_key):
        c, rt = self.c, self.rt
        paths, store = entity["paths"], rt.store(entity)
        ref = rt.plan["references"][entity["reference_id"]]
        source, sizes = Path(paths["accepted"]), paths["chrom_sizes"]
        common = {"accepted_signature": accepted_key, "assembly": ref["assembly"], "base": c["BASE_RESOLUTION"]}
        raw_key = rt.key("raw_cool", common, ["cooler"], ["matrix_products.py"])

        def matrix_stage(label, key, action, extra=None):
            dest = Path(paths[label])
            outputs = {"data": dest, "provenance": Path(str(dest) + ".provenance.json"), **(extra or {})}
            def produce(t, work):
                action(t, work)
                if dest.suffix == ".hic":
                    with t["data"].open("rb") as handle:
                        if handle.read(4) != b"HIC\0":
                            raise RuntimeError("Juicer returned an invalid .hic header")
                else:
                    rt.command(["cooler", "ls", str(t["data"])], rt.log(entity, label + "_validate"))
                json_write(t["provenance"], {"schema": 1, "stage": label, "signature": key, **common,
                    "accepted_pairs": str(source), "reference_id": entity["reference_id"], "members": entity["members"],
                    "dedup_units": entity["units"], "tool_versions": rt.versions,
                    "coordinates": "accepted pairs 1-based; cooler converts; Juicer short export explicitly uses pos-1"})
            return store.run(label, key, outputs, produce)

        def ensure_raw():
            def raw(t, work):
                rt.command(["cooler", "cload", "pairix", "--assembly", ref["assembly"], "-p", c["THREADS_MATRIX"],
                            f"{sizes}:{c['BASE_RESOLUTION']}", str(source), str(t["data"])], rt.log(entity, "cooler_cload"))
                rt.command(["cooler", "info", str(t["data"])], rt.log(entity, "cooler_info"), output=t["info"])
            matrix_stage("raw_cool", raw_key, raw, {"info": Path(paths["cool_info"])})

        with rt.slots["matrix"]:
            if truth(c["RUN_MCOOL"]):
                key = rt.key("norm_mcool", {**common, "resolutions": c["MCOOL_RESOLUTIONS"]}, ["cooler"], ["matrix_products.py"])
                dest = Path(paths["norm_mcool"])
                expected = {"data": dest, "provenance": Path(str(dest) + ".provenance.json")}
                if not store.ensure_reusable("norm_mcool", key, expected):
                    ensure_raw()
                matrix_stage("norm_mcool", key, lambda t, work: rt.command(["cooler", "zoomify", "-p", c["THREADS_MATRIX"],
                    "-r", c["MCOOL_RESOLUTIONS"], "--balance", "-o", str(t["data"]), paths["raw_cool"]], rt.log(entity, "zoomify_norm")))
                if truth(c["KEEP_RAW_MCOOL"]):
                    ensure_raw()
                    raw_mcool_key = rt.key("raw_mcool", {**common, "resolutions": c["MCOOL_RESOLUTIONS"]}, ["cooler"], ["matrix_products.py"])
                    matrix_stage("raw_mcool", raw_mcool_key, lambda t, work: rt.command(["cooler", "zoomify", "-p", c["THREADS_MATRIX"],
                        "-r", c["MCOOL_RESOLUTIONS"], "-o", str(t["data"]), paths["raw_cool"]], rt.log(entity, "zoomify_raw")))
            if truth(c["KEEP_SINGLE_RES_COOL"]) or not truth(c["RUN_MCOOL"]):
                ensure_raw()
                key = rt.key("norm_cool", common, ["cooler"], ["matrix_products.py"])
                def balanced(t, work):
                    shutil.copyfile(paths["raw_cool"], t["data"])
                    rt.command(["cooler", "balance", "-p", c["THREADS_MATRIX"], "-f", str(t["data"])], rt.log(entity, "balance"))
                matrix_stage("norm_cool", key, balanced)

        if builds_hic(entity, c):
            hic_args = {**common, "resolutions": c["HIC_RESOLUTIONS"], "jar": rt.jar_identity}
            norm_key = rt.key("norm_hic", {**hic_args, "normalizations": c["HIC_NORMALIZATIONS"]}, ["java"], ["matrix_products.py"])
            final = Path(paths["norm_hic"])
            if not store.ensure_reusable("norm_hic", norm_key, {"data": final, "provenance": Path(str(final) + ".provenance.json")}):
                with rt.slots["hic"]:
                    export_key = rt.key("juicer", accepted_key, [], ["matrix_products.py"])
                    store.run("juicer", export_key, {"pairs": Path(paths["juicer"])},
                              lambda t, work: export_juicer(source, t["pairs"]))
                    pre_key = rt.key("raw_hic", hic_args, ["java"], ["matrix_products.py"])
                    matrix_stage("raw_hic", pre_key, lambda t, work: rt.command(["java", f"-Xmx{c['JAVA_HEAP']}", "-jar", c["JUICER_TOOLS_JAR"],
                        "pre", "-n", "-t", str(work), "-r", c["HIC_RESOLUTIONS"], paths["juicer"], str(t["data"]), sizes], rt.log(entity, "juicer_pre")))
                    def norm(t, work):
                        shutil.copyfile(paths["raw_hic"], t["data"])
                        rt.command(["java", f"-Xmx{c['JAVA_HEAP']}", "-jar", c["JUICER_TOOLS_JAR"], "addNorm", "-j", c["THREADS_HIC_NORM"],
                                    "-k", c["HIC_NORMALIZATIONS"], str(t["data"])], rt.log(entity, "juicer_addNorm"))
                    matrix_stage("norm_hic", norm_key, norm)
            if c["PUBLIC_HIC_DIR"]:
                dest = Path(c["PUBLIC_HIC_DIR"]).resolve() / final.name
                store.run("public_hic", norm_key, {"hic": dest}, lambda t, work: shutil.copyfile(final, t["hic"]))
            if c["PUBLIC_HIC_BASE_URL"]:
                track = Path(paths["matrix_dir"]) / f"{entity['name']}.ucsc.hic.track.txt"
                track.write_text(f"# reference_id={entity['reference_id']} assembly={ref['assembly']}\n"
                    f"track type=hic name=\"{entity['name']}\" bigDataUrl={c['PUBLIC_HIC_BASE_URL'].rstrip('/')}/{final.name}\n", encoding="utf-8")
