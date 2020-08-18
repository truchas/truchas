from datetime import datetime, timezone
from json import load, dump
import os
import sys

dest_type = sys.argv[1]
version = sys.argv[2]
commit_sha = sys.argv[3]

url_base = "https://github.com/truchas/truchas_releases/releases/download"
tarball="truchas-%s-Linux.tar.bz2" % version
tarball_url = "{url_base}/{version}/{tarball}".format(
        url_base=url_base,
        version=version,
        tarball=tarball,
        )

assert dest_type in ["dev", "release"]
now_utc = datetime.now(timezone.utc)

filename = "data.json"
if os.path.exists(filename):
    d = load(open(filename))
else:
    d = {"dev": [], "release": []}
entry = {
    "url": tarball_url,
    "filename": tarball,
    "version": version,
    "commit_sha": commit_sha,
    "created": str(now_utc)
}
d[dest_type].append(entry)
print("Saving to %s." % filename)
with open(filename, "w") as f:
    dump(d, f, indent=4, ensure_ascii=False, sort_keys=True)
