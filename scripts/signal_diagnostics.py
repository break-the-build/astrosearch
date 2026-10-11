#!/usr/bin/env python3
"""Report event-pattern groups, depth consistency and variability conflicts
for one target's candidates (see vetting/signal_diagnostics.py).

Inputs are read only. The report goes to new files; candidate rows are never
modified, filtered or re-ranked, and no verdict is a classification.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from vetting import signal_diagnostics as sd
from vetting.io import atomic_write_text, read_table

csv.field_size_limit(2**31 - 1)


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def markdown(report: dict, title: str) -> str:
    lines = [f"# Signal diagnostics: {title}", "",
             "Report-only. Candidate rows are unchanged; verdicts name competing explanations, "
             "not classifications.", ""]
    for group in report["event_groups"]:
        members = group["members"]
        if len(members) == 1 and group["depth_check"]["status"] == "consistent":
            continue
        lines.append(f"## Event pattern {group['group']} ({len(members)} detections)")
        lines.append("")
        lines.append("| Detection | Period (d) | Events | Reported depth (ppm) | Estimator | Box depth | Local depth |")
        lines.append("|---|---|---|---|---|---|---|")
        for m in members:
            fmt = lambda v: "" if v is None else f"{v:.1f}"  # noqa: E731
            lines.append(f"| {m['detection']} | {m['period_days']:.8f} | {m['observed_events']} | "
                         f"{m['reported_depth_ppm']:.1f} | {m['depth_estimator']} | {fmt(m['box_depth_ppm'])} | "
                         f"{fmt(m['local_depth_ppm'])} |")
        lines.append("")
        for r in group["relations"]:
            lines.append(f"- {r['a']} and {r['b']}: {r['shared_events']} shared events "
                         f"({100 * r['overlap_fraction']:.0f}% of the smaller set), {r['period_relation']}")
        lines.append(f"- Depth check: **{group['depth_check']['status']}**")
        lines.extend(f"  - {issue}" for issue in group["depth_check"]["issues"])
        lines.append("")
    if report.get("interleaved_pairs"):
        lines += ["## Interleaved detections", ""]
        for pair in report["interleaved_pairs"]:
            a, b = pair["detections"]
            lines.append(f"- {a} and {b}: second events at phase {pair['phase_of_second_events']:.3f} of "
                         f"{pair['reference_period_days']:.8f} d, box depths {pair['box_depths_ppm'][0]:.1f} and "
                         f"{pair['box_depths_ppm'][1]:.1f} ppm. {pair['interpretation']}.")
        lines.append("")
    tested = [c for c in report["variability_conflicts"] if c.get("modulation_test", {}).get("status") == "OK"]
    if tested:
        lines += ["## Variability conflicts", "",
                  "| Detection | Period (d) | Modulation period (d) | Depth before (ppm) | Depth after (ppm) | Retained | Verdict |",
                  "|---|---|---|---|---|---|---|"]
        for c in tested:
            t = c["modulation_test"]
            retained = "" if t["retained_fraction"] is None else f"{100 * t['retained_fraction']:.0f}%"
            lines.append(f"| {c['detection']} | {c['period_days']:.8f} | {t['variability_period_days']:.8f} | "
                         f"{t['depth_before_ppm']:.1f} | {t['depth_after_ppm']:.1f} | {retained} | {t['verdict']} |")
        lines.append("")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--candidates", type=Path, required=True, help="astrosearch candidates.csv")
    p.add_argument("--lightcurve", type=Path, help="searched light curve (time, flux, ...); enables the "
                   "modulation test and TLS event matching")
    p.add_argument("--tls", type=Path, help="tls_search.py output CSV for the same light curve")
    p.add_argument("--output", type=Path, required=True, help="report JSON; a .md summary is written beside it")
    p.add_argument("--title", default="")
    a = p.parse_args(argv)
    if a.output.suffix.lower() != ".json":
        p.error("--output must end in .json")
    if a.tls and not a.lightcurve:
        p.error("--tls needs --lightcurve to locate observed TLS events")
    if a.output.resolve() in {x.resolve() for x in (a.candidates, a.lightcurve, a.tls) if x}:
        p.error("--output must not overwrite an input")
    time = flux = None
    if a.lightcurve:
        data = read_table(a.lightcurve)
        time, flux = data[:, 0], data[:, 1]
    report = sd.diagnose(read_rows(a.candidates), read_rows(a.tls) if a.tls else None, time, flux)
    report["inputs"] = {"candidates": str(a.candidates), "lightcurve": str(a.lightcurve or ""),
                        "tls": str(a.tls or "")}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_text(a.output, json.dumps(report, indent=1) + "\n")
    atomic_write_text(a.output.with_suffix(".md"), markdown(report, a.title or a.candidates.parent.name))
    inconsistent = sum(g["depth_check"]["status"] == "INCONSISTENT" for g in report["event_groups"])
    explained = sum("competing explanation" in c.get("modulation_test", {}).get("verdict", "")
                    for c in report["variability_conflicts"])
    print(f"{report['detections']} detections, {len(report['event_groups'])} event patterns, "
          f"{inconsistent} with inconsistent depths, {explained} dips explained by related modulation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
