"""Simplified basemap for the landing-page ancestry map.

Input: work/map/countries-50m.json (Natural Earth via world-atlas@2.0.2, public).
Output: work/map/basemap.json  {viewbox, highlight, paths: {<country>..., other}, project: {...}}
The map frame and highlighted countries come from family.json "ancestry" -> "map":
  {"bbox": [lon_west, lon_east, lat_south, lat_north], "countries": [Natural Earth names], "px_per_degree": 22}
Without it the whole world (minus Antarctica) is drawn with nothing highlighted.

Usage: python3 -I scripts/basemap.py
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import FAMILY, WORK  # noqa: E402
SRC = os.path.join(WORK, "map", "countries-50m.json")
OUT = os.path.join(WORK, "map", "basemap.json")

MAP = FAMILY.get("ancestry", {}).get("map", {})
LON0, LON1, LAT0, LAT1 = (float(v) for v in MAP.get("bbox", (-180, 180, -58, 84)))
K = float(MAP.get("px_per_degree", 22 if "bbox" in MAP else 2.5))  # px per degree latitude
COS = math.cos(math.radians((LAT0 + LAT1) / 2))
# Microstates drawn as part of the country that surrounds them.
ENCLAVES = {"San Marino": "Italy", "Vatican": "Italy", "Monaco": "France", "Lesotho": "South Africa"}
HIGHLIGHT = {c: c for c in MAP.get("countries", [])}
HIGHLIGHT.update({e: c for e, c in ENCLAVES.items() if c in HIGHLIGHT and e not in HIGHLIGHT})


def project(lon, lat):
    return (lon - LON0) * COS * K, (LAT1 - lat) * K


def decode_arcs(topo):
    sx, sy = topo["transform"]["scale"]
    tx, ty = topo["transform"]["translate"]
    arcs = []
    for arc in topo["arcs"]:
        x = y = 0
        pts = []
        for dx, dy in arc:
            x += dx
            y += dy
            pts.append((x * sx + tx, y * sy + ty))
        arcs.append(pts)
    return arcs


def ring(arcs, idxs):
    pts = []
    for i in idxs:
        a = arcs[i] if i >= 0 else arcs[~i][::-1]
        pts.extend(a if not pts else a[1:])
    return pts


def path_for(arcs, geom):
    polys = geom["arcs"] if geom["type"] == "MultiPolygon" else [geom["arcs"]]
    out = []
    for poly in polys:
        for r in poly:
            pts = ring(arcs, r)
            if not any(LON0 - 3 <= lo <= LON1 + 3 and LAT0 - 3 <= la <= LAT1 + 3 for lo, la in pts):
                continue
            xy, last = [], None
            for lo, la in pts:
                p = tuple(round(v, 1) for v in project(lo, la))
                if p != last:
                    xy.append(p)
                    last = p
            if len(xy) > 2:
                out.append("M" + "L".join(f"{x:g} {y:g}" for x, y in xy) + "Z")
    return "".join(out)


def main():
    topo = json.load(open(SRC))
    arcs = decode_arcs(topo)
    paths = {c: "" for c in MAP.get("countries", [])}
    paths["other"] = ""
    for g in topo["objects"]["countries"]["geometries"]:
        if g.get("type") not in ("Polygon", "MultiPolygon"):
            continue
        paths[HIGHLIGHT.get(g["properties"]["name"], "other")] += path_for(arcs, g)
    w, h = project(LON1, LAT0)
    json.dump({"viewbox": f"0 0 {w:.0f} {h:.0f}", "highlight": list(MAP.get("countries", [])), "paths": paths,
               "project": {"lon0": LON0, "lat1": LAT1, "cos": COS, "k": K},
               "source": "Natural Earth 1:50m via world-atlas@2.0.2"}, open(OUT, "w"))
    print(OUT, {k: len(v) for k, v in paths.items()})


if __name__ == "__main__":
    main()
