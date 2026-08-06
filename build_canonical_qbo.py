#!/usr/bin/env python3
"""Build the 'canonical' observed QBO cycle + occupancy statistics from the FU-Berlin
radiosonde series (Canton/Gan/Singapore, 70-10 hPa, monthly).

Method: deseasonalized anomalies -> 2-EOF phase space -> phase angle theta(t).
The composite cycle u(z, theta) is mapped to a TIME axis via the observed dwell time
per phase bin (occupancy x mean period), so the canonical cycle carries the real QBO's
temporal asymmetry: fast westerly descent, stalling easterlies, unequal quadrant
occupancy. Output (both written next to input/):
  input/qbo_canonical_cycle.txt -- u(level, month) anomaly table over one mean period,
                                   for the tabulated fake_qbo nudging target
  printed occupancy / duration / descent statistics for the paper

Usage: build_canonical_qbo.py qbo_fub.dat  (from https://www.geo.fu-berlin.de/met/ag/
strat/produkte/qbo/qbo.dat)
"""
import sys

import numpy as np

LEVS = np.array([70.0, 50.0, 40.0, 30.0, 20.0, 15.0, 10.0])  # hPa
NBIN = 56


def parse_fub(path):
    rec = {}
    for line in open(path, errors="ignore"):
        t = line.split()
        if len(t) < 9 or not t[0].isdigit() or len(t[1]) != 4 or not t[1].isdigit():
            continue
        yy, mm = int(t[1][:2]), int(t[1][2:])
        if not 1 <= mm <= 12:
            continue
        year = 1900 + yy if yy >= 53 else 2000 + yy
        vals = t[2::2] if len(t) >= 16 else t[2:9]
        if len(vals) != 7:
            continue
        try:
            u = np.array([float(v) for v in vals]) / 10.0  # 0.1 m/s -> m/s
        except ValueError:
            continue
        rec[(year, mm)] = u  # later stations overwrite earlier at overlaps
    keys = sorted(rec)
    t0, t1 = keys[0], keys[-1]
    n = (t1[0] - t0[0]) * 12 + (t1[1] - t0[1]) + 1
    U = np.full((n, 7), np.nan)
    months = np.empty(n, dtype=int)
    for (y, m), u in rec.items():
        i = (y - t0[0]) * 12 + (m - t0[1])
        U[i] = u
        months[i] = m
    for i in range(n):  # month-of-year for every slot (incl. any gaps)
        months[i] = (t0[1] - 1 + i) % 12 + 1
    return U, months, t0, t1


def main():
    U, moy, t0, t1 = parse_fub(sys.argv[1] if len(sys.argv) > 1 else "qbo_fub.dat")
    ok = ~np.isnan(U).any(axis=1)
    print(f"record {t0[0]}-{t0[1]:02d} .. {t1[0]}-{t1[1]:02d}: {ok.sum()}/{len(U)} complete months")
    # deseasonalize per level
    A = U.copy()
    for m in range(1, 13):
        sel = ok & (moy == m)
        A[moy == m] -= np.nanmean(U[sel], axis=0)
    A[~ok] = np.nan

    # EOF phase space on standardized anomalies (complete months only)
    S = (A[ok] - A[ok].mean(0)) / A[ok].std(0)
    _, _, Vt = np.linalg.svd(S, full_matrices=False)
    pcs = S @ Vt[:2].T
    pcs /= pcs.std(0)
    theta = np.unwrap(np.arctan2(pcs[:, 1], pcs[:, 0]))
    if np.median(np.diff(theta)) < 0:  # orient so phase advances forward in time
        theta = -theta
    ncyc = (theta[-1] - theta[0]) / (2 * np.pi)
    period = ok.sum() / ncyc
    print(f"{ncyc:.1f} cycles -> mean period {period:.1f} months")

    # occupancy + composite per phase bin
    th = np.mod(theta, 2 * np.pi)
    bins = np.floor(th / (2 * np.pi) * NBIN).astype(int)
    occ = np.bincount(bins, minlength=NBIN) / len(bins)
    comp = np.array([A[ok][bins == b].mean(0) for b in range(NBIN)])  # (NBIN, 7)

    # roll so the cycle starts at the westerly-onset at 10 hPa (u crosses +)
    u10 = comp[:, LEVS == 10.0].ravel()
    start = int(np.argmax((u10 > 0) & (np.roll(u10, 1) <= 0)))
    comp, occ = np.roll(comp, -start, axis=0), np.roll(occ, -start)

    # map phase bins -> canonical (asymmetric) time axis, resample to monthly grid
    tedge = np.concatenate([[0.0], np.cumsum(occ)]) * period
    tmid = 0.5 * (tedge[:-1] + tedge[1:])
    nmon = int(round(period))
    tgrid = (np.arange(nmon) + 0.5) * period / nmon
    tper = np.concatenate([tmid - period, tmid, tmid + period])
    cyc = np.vstack([np.interp(tgrid, tper, np.tile(comp[:, j], 3)) for j in range(7)])
    cyc -= cyc.mean(axis=1, keepdims=True)  # exactly zero time-mean forcing per level

    out = "input/qbo_canonical_cycle.txt"
    with open(out, "w") as f:
        f.write("# canonical observed QBO cycle (FUB radiosonde composite, anomaly u m/s)\n")
        f.write(f"# built by build_canonical_qbo.py; mean period {period:.1f} mo, "
                f"{ncyc:.1f} cycles {t0[0]}-{t1[0]}; phase origin: westerly onset at 10 hPa\n")
        f.write(f"{len(LEVS)} {nmon}\n")
        f.write(" ".join(f"{p:.0f}" for p in LEVS) + "\n")
        for j in range(len(LEVS)):
            f.write(" ".join(f"{v:7.2f}" for v in cyc[j]) + "\n")
    print(f"wrote {out}  ({len(LEVS)} levels x {nmon} months)")

    # --- occupancy / asymmetry statistics ---
    quad = np.floor(np.mod(np.arctan2(pcs[:, 1], pcs[:, 0]), 2 * np.pi) / (np.pi / 2)).astype(int)
    qf = np.bincount(quad, minlength=4) / len(quad)
    print("\nphase-quadrant occupancy (symmetric would be 25% each):")
    print("  " + "  ".join(f"Q{i+1} {100*f:.1f}%" for i, f in enumerate(qf)))
    for p in (50.0, 30.0, 20.0):
        a = A[ok][:, LEVS == p].ravel()
        w = a > 0
        runs = np.diff(np.flatnonzero(np.concatenate([[True], np.diff(w), [True]])))
        labs = w[np.concatenate([[0], np.cumsum(runs)[:-1]])]
        wd = runs[labs].mean(); ed = runs[~labs].mean()
        print(f"  {p:.0f} hPa: westerly {100*w.mean():.1f}% of months "
              f"(mean spell {wd:.1f} mo) | easterly {100*(1-w.mean()):.1f}% (mean spell {ed:.1f} mo)")
    # descent-rate asymmetry from the canonical cycle: zero-crossing lag 20 -> 50 hPa
    for lab, sgn in (("westerly", 1.0), ("easterly", -1.0)):
        lag = {}
        for p in (20.0, 50.0):
            uu = cyc[np.flatnonzero(LEVS == p)[0]]
            cross = np.flatnonzero((sgn * uu > 0) & (np.roll(sgn * uu, 1) <= 0))
            lag[p] = tgrid[cross[0]] if len(cross) else np.nan
        dt = (lag[50.0] - lag[20.0]) % period
        print(f"  canonical {lab} onset descent 20->50 hPa: {dt:.1f} months")


if __name__ == "__main__":
    main()
