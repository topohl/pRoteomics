"""Preserve the analysis version held in the pRoteomics working tree (normalisation review 2026-10-21).

Copies exports/, data/processed/ and results/ to a sibling archive folder with robocopy, hashes every copied file
against its source (SHA256_MANIFEST.csv), and stores the uncommitted code state as a patch against a base commit
plus copies of the untracked code files, so the version survives without a commit. Nothing in the working tree is
changed. Run before regenerating the working tree under another sample centring (config/normalisation.yml).

USAGE  python tools/normalisation_review/snapshot_working_version.py --dest <folder> [--base b478a13] [--exclude <rel dir>]...
       add --code-only to refresh only the code state of an existing snapshot
       --trees <rel dir>... copies only those trees of the working tree (default: exports, data/processed, results);
       --external <path>::<rel dest> adds a folder from outside the repository (repeatable). Crash dumps (*.stackdump)
       and Excel lock files (~$*) are not copied.
"""
import argparse, csv, hashlib, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

SRC = r"S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress\Analysis\proteomics"
TREES = ["exports", os.path.join("data", "processed"), "results"]


def lp(p):
    p = os.path.abspath(p)
    return p if p.startswith("\\\\?\\") else ("\\\\?\\UNC\\" + p[2:] if p.startswith("\\\\") else "\\\\?\\" + p)


def sha(path):
    h = hashlib.sha256()
    with open(lp(path), "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", required=True)
    ap.add_argument("--base", default="b478a13")
    ap.add_argument("--exclude", action="append", default=[])
    ap.add_argument("--code-only", action="store_true")
    ap.add_argument("--trees", nargs="+", default=TREES)
    ap.add_argument("--external", action="append", default=[])
    a = ap.parse_args()
    trees = [t.replace("/", os.sep) for t in a.trees]
    ext = [tuple(e.split("::", 1)) for e in a.external]
    jobs = [(os.path.join(SRC, t), t) for t in trees] + [(src, rel.replace("/", os.sep)) for src, rel in ext]
    os.makedirs(a.dest, exist_ok=True)
    logs = os.path.join(a.dest, "_robocopy_logs"); os.makedirs(logs, exist_ok=True)
    excl = [os.path.join(SRC, e.replace("/", os.sep)) for e in a.exclude]
    for src, t in ([] if a.code_only else jobs):
        cmd = ["robocopy", src, os.path.join(a.dest, t), "/E", "/COPY:DAT", "/DCOPY:T", "/R:2", "/W:2",
               "/MT:16", "/NFL", "/NDL", "/NP", "/XF", "*.stackdump", "~$*", "/LOG:" + os.path.join(logs, t.replace(os.sep, "_") + ".log")]
        if excl:
            cmd += ["/XD"] + excl
        rc = subprocess.call(cmd)
        print(t, "robocopy exit", rc, flush=True)
        if rc >= 8:
            sys.exit(f"robocopy failed for {t} (exit {rc})")
    # code state: patch against the base commit (tracked changes, binary-safe) and the untracked code files
    code = os.path.join(a.dest, "code"); os.makedirs(code, exist_ok=True)
    git = lambda *args: subprocess.run(["git", "-C", SRC, *args], capture_output=True)
    open(os.path.join(code, "HEAD.txt"), "wb").write(git("rev-parse", "HEAD").stdout)
    open(os.path.join(code, "status_porcelain.txt"), "wb").write(git("status", "--porcelain", "--untracked-files=all").stdout)
    open(os.path.join(code, f"working_tree_vs_{a.base}.patch"), "wb").write(git("diff", "--binary", a.base).stdout)
    untracked = [u for u in git("ls-files", "--others", "--exclude-standard").stdout.decode().splitlines()
                 if u.split("/")[0] in ("R", "analysis", "tools", "config", "archive", "audits", "docs", "tests", "inst")
                 or u in ("pipeline.yml", "run_dataset_pipeline.R", "RUN_ORDER.md")]
    for u in untracked:
        d = os.path.join(code, "untracked", u.replace("/", os.sep)); os.makedirs(os.path.dirname(d), exist_ok=True)
        open(d, "wb").write(open(os.path.join(SRC, u.replace("/", os.sep)), "rb").read())
    open(os.path.join(code, "untracked_files.txt"), "w", encoding="utf-8").write("\n".join(untracked) + "\n")
    if a.code_only:
        print("code state refreshed:", len(untracked), "untracked code files; manifest left unchanged", flush=True)
        return
    # hash-verify every copied file against its source
    rels = []; source_of = {}
    for src, t in jobs:
        root = os.path.join(a.dest, t)
        for dp, _, fs in os.walk(lp(root)):
            for f in fs:
                full = os.path.join(dp, f)
                rel = os.path.relpath(full, lp(a.dest)); rels.append(rel)
                source_of[rel] = os.path.join(src, os.path.relpath(full, lp(root)))
    print("files to verify:", len(rels), flush=True)

    def job(rel):
        try:
            ha = sha(os.path.join(a.dest, rel)); size = os.path.getsize(lp(os.path.join(a.dest, rel)))
        except Exception as e:
            return rel, -1, "", "", f"archive_error:{e}"
        try:
            hs = sha(source_of[rel])
        except Exception as e:
            return rel, size, ha, "", f"source_error:{e}"
        return rel, size, ha, hs, "match" if ha == hs else "mismatch"

    counts = {}
    with open(os.path.join(a.dest, "SHA256_MANIFEST.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["relative_path", "bytes", "sha256", "source_sha256", "status"])
        with ThreadPoolExecutor(16) as ex:
            for i, row in enumerate(ex.map(job, rels), 1):
                w.writerow(row); counts[row[4].split(":")[0]] = counts.get(row[4].split(":")[0], 0) + 1
                if i % 20000 == 0:
                    print(f"{i}/{len(rels)} {counts}", flush=True)
    summary = (f"Snapshot of {SRC} trees {trees}, external {ext} (excluded: {a.exclude}, *.stackdump, ~$*).\nFiles: {len(rels)}; status counts: {counts}\n"
               f"Manifest: SHA256_MANIFEST.csv\nCode: code/HEAD.txt, code/working_tree_vs_{a.base}.patch, code/untracked/ ({len(untracked)} files)\n")
    open(os.path.join(a.dest, "ARCHIVE_SUMMARY.txt"), "w", encoding="utf-8").write(summary)
    print("DONE", counts, flush=True)


if __name__ == "__main__":
    main()
