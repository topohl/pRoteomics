"""Re-hash a random sample of a preserved archive against its SHA256 manifest (immutability spot check).

Writes one row per sampled file (path, manifest hash, hash now, status) to the given CSV.
USAGE  python tools/normalisation_review/spot_check_archive.py --archive <folder> --out <csv> [--n 3000] [--seed 20261021]
"""
import argparse, csv, hashlib, os, random


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
    ap.add_argument("--archive", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=20261021)
    a = ap.parse_args()
    with open(os.path.join(a.archive, "SHA256_MANIFEST.csv"), newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    cols = list(rows[0].keys())
    pcol = next(c for c in cols if "path" in c.lower())
    hcol = next(c for c in cols if c.lower() in ("sha256", "archive_sha256"))
    random.seed(a.seed)
    sample = random.sample(rows, min(a.n, len(rows)))
    bad = 0
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["relative_path", "manifest_sha256", "sha256_now", "status"])
        for r in sample:
            try:
                h = sha(os.path.join(a.archive, r[pcol]))
                st = "match" if h == r[hcol] else "mismatch"
            except Exception as e:
                h, st = "", f"error: {e}"
            bad += st != "match"
            w.writerow([r[pcol], r[hcol], h, st])
    print(f"archive files in manifest: {len(rows)}; sampled: {len(sample)}; mismatching or unreadable: {bad}")


if __name__ == "__main__":
    main()
