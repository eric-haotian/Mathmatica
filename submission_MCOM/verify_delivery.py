#!/usr/bin/env python3
"""Verify the delivered transfer package without executing its scientific code."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import sys
import zipfile

def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "MANIFEST_SHA256.json").read_text(encoding="utf-8"))
    failed = []
    for name, expected in manifest["files"].items():
        rel = PurePosixPath(name)
        if rel.is_absolute() or ".." in rel.parts:
            failed.append(name + ": invalid path")
            continue
        path = root.joinpath(*rel.parts)
        if not path.is_file() or path.is_symlink():
            failed.append(name + ": missing or symbolic file")
            continue
        got = hashlib.sha256(path.read_bytes()).hexdigest()
        if got != expected["sha256"] or path.stat().st_size != expected["bytes"]:
            failed.append(name + ": byte/hash mismatch")
    archive = root / "upload/03_Reproducibility_Materials.zip"
    if archive.is_file():
        with zipfile.ZipFile(archive) as z:
            bad = z.testzip()
            if bad is not None:
                failed.append("supplement CRC: " + bad)
    if failed:
        print("FAIL\n" + "\n".join(failed), file=sys.stderr)
        return 1
    print("PASS: " + str(len(manifest["files"])) + " delivered file hashes; supplement CRC.")
    print("To rebuild PDFs: cd source && bash build.sh")
    print("To replay scientific certificates: extract the supplement and run python3 run_all.py")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
