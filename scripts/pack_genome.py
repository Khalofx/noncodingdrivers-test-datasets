#!/usr/bin/env python3
"""Pack a bgdata genome directory into a reproducible .tar.gz for the pipeline test profile.
Entries are sorted, owned by 0:0 with fixed permissions and mtime, and the gzip header carries no
timestamp or file name, so the same directory always gives a byte-identical archive.
Usage: pack_genome.py <genome_dir> <out.tar.gz>"""
import gzip, io, pathlib, sys, tarfile

MTIME = 1792454400  # 2026-10-20T00:00:00Z
src, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
buf = io.BytesIO()
with tarfile.open(fileobj=buf, mode="w", format=tarfile.GNU_FORMAT) as tar:
    for path in [src] + sorted(src.rglob("*")):
        info = tar.gettarinfo(str(path), arcname=str(pathlib.Path(src.name) / path.relative_to(src)))
        info.uid = info.gid = 0
        info.uname = info.gname = ""
        info.mtime = MTIME
        info.mode = 0o755 if info.isdir() else 0o644
        if info.isfile():
            with open(path, "rb") as fh:
                tar.addfile(info, fh)
        else:
            tar.addfile(info)
with open(out, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as gz:
    gz.write(buf.getvalue())
print(f"{out}: {out.stat().st_size} bytes")
