#!/usr/bin/env python3
"""Assemble and verify the Mathematics of Computation submission package.

Run from anywhere:  python3 submission_MCOM/tools/assemble_package.py

Steps (all repeatable; nothing is downloaded):
  1. Build main.pdf, cover_letter.pdf and reproducibility_index.pdf with source/build.sh,
     and rebuild them again in an isolated temporary copy of source/ to confirm that the
     page text is identical (build reproducibility).
  2. Preflight every PDF: page count, page size, embedded fonts, black text only, no
     text outside the page box.
  3. Compare source/main.tex with the previous version (author_checks/previous_versions):
     every theorem-like statement, proof, inherited display, Algorithm 1 and the
     stopping-path table must be byte-identical; record counts, labels and citations.
  4. Regenerate the SHA-256 manifest of the reproducibility package (only metadata files
     may differ from the previous manifest), rebuild upload/03_Reproducibility_Materials.zip,
     replay the frozen evidence from the unzipped tree and from a fresh extraction of the
     new ZIP, and record both replay summaries.
  5. Copy the PDFs to upload/, write the author_checks records, regenerate the outer
     MANIFEST_SHA256.json, and run verify_delivery.py.
"""
from __future__ import annotations

import difflib
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]              # submission_MCOM/
SRC = ROOT / "source"
BUILD = SRC / "build"
UPLOAD = ROOT / "upload"
CHECKS = ROOT / "author_checks"
REPRO = ROOT.parent / "reproducibility" / "Supplementary_Materials"
PREVIOUS = CHECKS / "previous_versions" / "main_2026-09-29.tex"
TODAY = "2026-10-01"
NAMES = ("main", "cover_letter", "reproducibility_index")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(cmd: list[str], cwd: Path, log: Path | None = None) -> subprocess.CompletedProcess:
    res = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    if log is not None:
        log.write_text(res.stdout + res.stderr, encoding="utf-8")
    if res.returncode:
        sys.stderr.write(res.stdout + res.stderr)
        raise SystemExit(f"command failed: {' '.join(cmd)}")
    return res


def pdftotext(pdf: Path) -> str:
    return subprocess.run(["pdftotext", "-layout", str(pdf), "-"], text=True, capture_output=True, check=True).stdout


# ----------------------------------------------------------------------------- 1. build
def build_all() -> dict:
    CHECKS.joinpath("compile_logs").mkdir(parents=True, exist_ok=True)
    shutil.rmtree(BUILD, ignore_errors=True)
    run(["bash", "build.sh"], SRC, CHECKS / "isolated_build_console.log")
    for n in NAMES:
        shutil.copy2(BUILD / f"{n}.log", CHECKS / "compile_logs" / f"{n}.log")
    # isolated rebuild
    with tempfile.TemporaryDirectory(prefix="mcom_isolated_") as tmp:
        iso = Path(tmp) / "source"
        shutil.copytree(SRC, iso, ignore=shutil.ignore_patterns("build", "__pycache__", "*.png"))
        run(["bash", "build.sh"], iso)
        identical = {n: pdftotext(iso / "build" / f"{n}.pdf") == pdftotext(BUILD / f"{n}.pdf") for n in NAMES}
    for n in NAMES:
        (CHECKS / f"{n}_extracted.txt").write_text(pdftotext(BUILD / f"{n}.pdf"), encoding="utf-8")
    return identical


# ------------------------------------------------------------------------- 2. preflight
def preflight(pdf: Path, isolated_identical: bool) -> dict:
    import pymupdf as fitz

    doc = fitz.open(pdf)
    rec = {
        "pages": doc.page_count,
        "page_size_points": [round(v, 1) for v in doc[0].rect],
        "all_text_black": True,
        "text_outside_page": [],
        "all_font_programs_embedded": True,
        "font_count": 0,
        "compile_errors_and_layout_warnings": [],
        "isolated_rebuild_page_text_identical": isolated_identical,
    }
    fonts = set()
    for pno, page in enumerate(doc, 1):
        for f in page.get_fonts(full=True):
            fonts.add(f[3])
        d = page.get_text("dict")
        for block in d["blocks"]:
            for line in block.get("lines", []):
                for span in line["spans"]:
                    if span["color"] != 0:
                        rec["all_text_black"] = False
                    bbox = fitz.Rect(span["bbox"])
                    if not page.rect.contains(bbox):
                        rec["text_outside_page"].append({"page": pno, "text": span["text"][:40]})
    rec["font_count"] = len(fonts)
    pf = subprocess.run(["pdffonts", str(pdf)], text=True, capture_output=True, check=True).stdout.splitlines()[2:]
    rec["all_font_programs_embedded"] = all(line.split()[-5] == "yes" for line in pf if line.strip())
    log = (BUILD / f"{pdf.stem}.log").read_text(encoding="utf-8", errors="replace")
    for pat in (r"^!.*", r"LaTeX Warning: (Reference|Citation).*undefined", r"Overfull \\hbox.*", r"Missing character:.*"):
        rec["compile_errors_and_layout_warnings"] += re.findall(pat, log, re.M)
    rec["other_log_warnings"] = sorted(set(re.findall(r"(Underfull \\vbox.*|pdfTeX warning \(font expansion\).*)", log)))
    return rec


# ------------------------------------------------------------------- 3. source integrity
def envs(text: str, name: str) -> list[str]:
    return re.findall(r"\\begin\{%s\}.*?\\end\{%s\}" % (name, name), text, re.S)


def displays(text: str) -> list[str]:
    return re.findall(r"\\begin\{(equation|align|gather|multline)\}.*?\\end\{\1\}|\\\[.*?\\\]", text, re.S)


def source_integrity() -> dict:
    old = PREVIOUS.read_text(encoding="utf-8")
    new = (SRC / "main.tex").read_text(encoding="utf-8")
    rec: dict = {"baseline": PREVIOUS.name, "baseline_sha256": sha256(PREVIOUS)}
    for name in ("theorem", "lemma", "proposition", "corollary", "remark", "proof", "algorithm"):
        a, b = envs(old, name), envs(new, name)
        rec[name] = {"baseline": len(a), "new": len(b), "byte_identical": a == b}
    da, db = displays(old), displays(new)
    pos, ordered = 0, True
    for d in da:
        i = new.find(d, pos)
        if i < 0:
            ordered = False
            break
        pos = i + len(d)
    rec["display_math"] = {"baseline_blocks": len(da), "new_blocks": len(db), "all_inherited_blocks_byte_identical_in_order": ordered}
    ta, tb = re.findall(r"\\begin\{table\}.*?\\end\{table\}", old, re.S), re.findall(r"\\begin\{table\}.*?\\end\{table\}", new, re.S)
    rec["table"] = {"baseline": len(ta), "new": len(tb), "baseline_table_present_verbatim": all(t in new for t in ta)}
    rec["figure"] = {"baseline": len(envs(old, "figure")), "new": len(envs(new, "figure"))}
    lo, ln = set(re.findall(r"\\label\{([^}]*)\}", old)), set(re.findall(r"\\label\{([^}]*)\}", new))
    refs = set(re.findall(r"\\(?:eqref|ref|cref)\{([^}]*)\}", new))
    rec["cross_reference_labels"] = {"baseline": len(lo), "new": len(ln), "added": sorted(ln - lo), "removed": sorted(lo - ln), "unresolved": sorted(refs - ln)}
    bo, bn = re.findall(r"\\bibitem\{([^}]*)\}", old), re.findall(r"\\bibitem\{([^}]*)\}", new)
    cites = {c.strip() for m in re.findall(r"\\cite(?:\[[^\]]*\])?\{([^}]*)\}", new) for c in m.split(",")}
    body_old = {k: re.search(r"\\bibitem\{%s\}(.*)" % re.escape(k), old).group(1) for k in bo}
    body_new = {k: re.search(r"\\bibitem\{%s\}(.*)" % re.escape(k), new).group(1) for k in bn}
    rec["bibliography"] = {
        "baseline_count": len(bo), "new_count": len(bn), "added": sorted(set(bn) - set(bo)), "removed": sorted(set(bo) - set(bn)),
        "edited_existing_entries": sorted(k for k in bo if k in body_new and body_new[k] != body_old[k]),
        "uncited": sorted(set(bn) - cites), "missing": sorted(cites - set(bn)),
    }
    abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", new, re.S).group(1)
    rec["abstract_words_whitespace"] = len(abstract.split())
    rec["source_sha256"] = sha256(SRC / "main.tex")
    diff = difflib.unified_diff(old.splitlines(keepends=True), new.splitlines(keepends=True), "MCOM_2026-09-29/source/main.tex", "MCOM_2026-10-01/source/main.tex")
    (CHECKS / "MCOM_2026-09-29_to_2026-10-01.diff").write_text("".join(diff), encoding="utf-8")
    return rec


# ----------------------------------------------------------- 4. reproducibility package
def repro_manifest_and_zip() -> dict:
    manifest_path = REPRO / "MANIFEST_SHA256.json"
    previous = json.loads(manifest_path.read_text(encoding="utf-8"))["files"]
    files = sorted(p for p in REPRO.rglob("*") if p.is_file() and p.name != "MANIFEST_SHA256.json"
                   and "__pycache__" not in p.parts and "replay_output" not in p.parts)
    current = {p.relative_to(REPRO).as_posix(): sha256(p) for p in files}
    changed = sorted(k for k in current if k in previous and previous[k] != current[k])
    added = sorted(set(current) - set(previous))
    removed = sorted(set(previous) - set(current))
    allowed = {"SECTION_MAP.txt", "TRANSFER_NOTE.txt"}
    if added or removed or not set(changed) <= allowed:
        raise SystemExit(f"unexpected supplement changes: changed={changed} added={added} removed={removed}")
    manifest_path.write_text(json.dumps({"algorithm": "SHA-256", "files": current}, indent=2) + "\n", encoding="utf-8")
    # deterministic zip
    zpath = UPLOAD / "03_Reproducibility_Materials.zip"
    UPLOAD.mkdir(exist_ok=True)
    members = sorted(p for p in REPRO.rglob("*") if p.is_file() and "__pycache__" not in p.parts and "replay_output" not in p.parts)
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for p in members:
            info = zipfile.ZipInfo("Supplementary_Materials/" + p.relative_to(REPRO).as_posix(), date_time=(2026, 10, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            z.writestr(info, p.read_bytes())
    with zipfile.ZipFile(zpath) as z:
        assert z.testzip() is None
        n_members = len(z.namelist())
    # replay from the unzipped tree and from a fresh extraction
    for out in ("evidence_replay", "zip_evidence_replay"):
        shutil.rmtree(CHECKS / out, ignore_errors=True)
    run([sys.executable, "run_all.py", "--output", str(CHECKS / "evidence_replay")], REPRO, CHECKS / "evidence_replay_console.log")
    with tempfile.TemporaryDirectory(prefix="mcom_zip_replay_") as tmp:
        with zipfile.ZipFile(zpath) as z:
            z.extractall(tmp)
        run([sys.executable, "run_all.py", "--output", str(CHECKS / "zip_evidence_replay")], Path(tmp) / "Supplementary_Materials", CHECKS / "zip_replay_console.log")
    summary = json.loads((CHECKS / "zip_evidence_replay" / "SUMMARY.json").read_text())
    unchanged_payload = sum(1 for k in current if k.startswith("code/"))
    return {
        "new_archive_sha256": sha256(zpath), "zip_members": n_members, "manifest_entries": len(current),
        "metadata_files_changed_since_2026-09-29": changed, "code_and_data_files": unchanged_payload,
        "code_and_data_files_changed": [k for k in changed if k.startswith("code/")], "zip_crc": "PASS",
        "fresh_extraction_replay": {k: summary[k] for k in ("status", "independent_certificate_tasks", "checks_including_rejections", "rejection_checks", "manifest_files_checked")},
    }


# --------------------------------------------------------------------- 5. deliverables
def outer_manifest() -> int:
    skip_dirs = {"build", "__pycache__", ".git"}
    files = []
    for p in sorted(ROOT.rglob("*")):
        if not p.is_file() or set(p.relative_to(ROOT).parts) & skip_dirs or p.name == "MANIFEST_SHA256.json" or p.suffix == ".png":
            continue
        files.append(p)
    entries = {p.relative_to(ROOT).as_posix(): {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in files}
    (ROOT / "MANIFEST_SHA256.json").write_text(json.dumps({"algorithm": "SHA-256", "files": entries}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return len(entries)


def main() -> None:
    t0 = time.perf_counter()
    identical = build_all()
    print("built; isolated rebuild identical:", identical)
    pre = {n: preflight(BUILD / f"{n}.pdf", identical[n]) for n in NAMES}
    (CHECKS / "pdf_preflight.json").write_text(json.dumps(pre, indent=2) + "\n", encoding="utf-8")
    for n, rec in pre.items():
        if rec["compile_errors_and_layout_warnings"] or not rec["all_font_programs_embedded"] or not rec["all_text_black"] or rec["text_outside_page"]:
            raise SystemExit(f"preflight failed for {n}: {rec}")
    print("preflight:", {n: (r["pages"], r["font_count"]) for n, r in pre.items()})
    integ = source_integrity()
    (CHECKS / "source_integrity.json").write_text(json.dumps(integ, indent=2) + "\n", encoding="utf-8")
    for name in ("theorem", "lemma", "proposition", "corollary", "remark", "proof", "algorithm"):
        if not integ[name]["byte_identical"]:
            raise SystemExit(f"{name} environments changed")
    if not integ["display_math"]["all_inherited_blocks_byte_identical_in_order"] or not integ["table"]["baseline_table_present_verbatim"] or integ["cross_reference_labels"]["unresolved"] or integ["bibliography"]["missing"]:
        raise SystemExit(f"integrity failed: {integ}")
    print("source integrity: PASS; abstract words", integ["abstract_words_whitespace"])
    sup = repro_manifest_and_zip()
    (CHECKS / "supplement_transfer.json").write_text(json.dumps(sup, indent=2) + "\n", encoding="utf-8")
    print("supplement:", sup["fresh_extraction_replay"])
    for n, out in zip(NAMES, ("01_Main_Manuscript.pdf", "02_Cover_Letter.pdf", "04_Reproducibility_Index.pdf")):
        shutil.copy2(BUILD / f"{n}.pdf", UPLOAD / out)
    verification = {
        "status": "PASS", "date": TODAY, "target_journal": "Mathematics of Computation",
        "title": "Frequency Saturation and Certified Computation of Relaxation Functionals",
        "basis": {"previous_version": "MCOM transfer package of 2026-09-29", "previous_main_tex_sha256": integ["baseline_sha256"], "author": "Haotian Zhong"},
        "changes": [
            "abstract and introduction rewritten: application setting, problem statement, main results, classical and applied related work, organization and notation",
            "lead-in paragraphs for Sections 2-6; definition of the baseline prediction distance stated before Theorem 4.2",
            "Section 7: Tables 1, 2, 4 and Figure 1 added from archived certified quantities; stopping-path table is now Table 3",
            "conclusions extended with limitations; AI-use declaration reworded with the same scope",
            "five references added (Karlin-Studden, Krein-Nudel'man, Orazem-Tribollet, Ciucci-Chen, Wan et al.); Wasilkowski-Wozniakowski entry format harmonized",
            "all theorem, lemma, proposition, corollary, remark and proof environments, all inherited displays, Algorithm 1 and the stopping table are byte-identical to the previous version",
        ],
        "source_integrity": integ, "pdf_preflight": pre, "supplement": sup,
        "upload_files": {p.name: {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in sorted(UPLOAD.iterdir())},
        "seconds": round(time.perf_counter() - t0, 1),
    }
    (CHECKS / "TRANSFER_VERIFICATION.json").write_text(json.dumps(verification, indent=2) + "\n", encoding="utf-8")
    n = outer_manifest()
    print("outer manifest entries:", n)
    run([sys.executable, "verify_delivery.py"], ROOT)
    print("verify_delivery.py: PASS")


if __name__ == "__main__":
    main()
