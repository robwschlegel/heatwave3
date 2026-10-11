"""Regenerate the figures and printed outputs used by getting-started.md.

Runs the Getting Started examples on the test data bundled with the
repository and writes the results to python/docs/_generated/, which the page
includes. Run from anywhere after installing heatwave3[plot]:

    python python/docs/generate_examples.py
"""

import os
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

import heatwave3 as hw3  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "_generated"
OUT.mkdir(exist_ok=True)
sst_file = str(REPO / "inst" / "extdata" / "sst_test.nc")
pd.set_option("display.width", 100, "display.max_columns", 8)


def save_text(name, text):
    (OUT / f"{name}.txt").write_text(str(text).rstrip() + "\n")


def save_fig(ax, name):
    ax.figure.savefig(OUT / f"{name}.png", dpi=110, bbox_inches="tight")
    plt.close(ax.figure)


stem = os.path.join(tempfile.mkdtemp(), "demo")
out = hw3.detect3(sst_file, stem, ("1982-01-01", "2011-12-31"), category=True,
                  n_threads=2, quiet=True)
clim_file, event_file = out["clim"], out["events"]

events = hw3.hw3_export(event_file)
save_text("events_head", events[["lon", "lat", "event_no", "date_start", "date_end",
                                 "duration", "intensity_max", "category"]].head())

save_fig(hw3.event_line3(sst_file, clim_file, lon=26.525, lat=-34.125,
                         start_date="2010-01-01", end_date="2012-12-31"), "event_line")
save_fig(hw3.plot_metric3(event_file, metric="intensity_max", summary="mean",
                          coastline=False), "metric_map")

cats = hw3.category3(event_file, clim_file)
save_text("category_head", cats.head())
save_text("category_counts", cats["category"].value_counts())
save_text("block_head", hw3.block_average3(event_file)[
    ["lon", "lat", "year", "count", "duration_mean", "intensity_max_max", "total_days"]].head())

blobs = hw3.detect_blob3(sst_file, clim_file, minVoxels=20, topN=3,
                         return_=("event", "daily"))
save_text("blob_events", blobs["event"][["blob", "duration", "date_start", "date_end",
                                         "peakArea_km2", "cumI_km2_day", "rank"]])
fig, ax = plt.subplots(figsize=(8, 3))
for _, ev in blobs["event"].iterrows():
    d = blobs["daily"][blobs["daily"]["blob"] == ev["blob"]]
    ax.plot(d["date"], d["area_km2"], marker="o", ms=2, label=f"#{ev['rank']}")
ax.set(ylabel="Area (km²)", title="Daily area of the top 3 blobs")
ax.legend(title="rank")
save_fig(ax, "blob_area")

print(f"Wrote {len(list(OUT.iterdir()))} files to {OUT}")
