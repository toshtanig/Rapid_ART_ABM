#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Progress monitor for the revision (R1) runs, with a remaining-time estimate.

Reads the output folders under revision_R1_runs/ (synchronised across the workstations used) and
reports, for every analysis, setting and terminal, how many allocated seeds are
done, which seed is running, and which seeds have stalled (folder present but
no recent write and no CEA/break-even output). Minutes per seed are estimated
for each terminal from the spacing of its recent completions, and the expected
finish time of each analysis is projected assuming each terminal works through
RUN_ORDER without pauses.

The seed lists, seed targets and setting names are read directly from
run_revision_scenarios.py (no import, so the model and numba are not loaded):
this script needs only the Python standard library and runs from VS Code on
Windows or WSL. Run as-is; all options are in CONFIG.
"""

import os
import sys
import ast
import time
import statistics

# =============================================================================
# CONFIG
# =============================================================================
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import REVISION_RUNS

DRIVER_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_revision_scenarios.py")

RUN_ORDER = ["A11_timing_convention", "A12_age_offset", "A1_factorial", "A8_status_quo_empirical", "A9_observed_ineligible_wait", "A10_conservative_combined",
             "A5_post_diagnosis_behaviour", "A2_foi_trajectory", "A7_delayed_start", "A4_time_to_suppression"]

STALL_MINUTES = 90            # unfinished seed with no write for longer than this -> "stalled"
IDLE_MINUTES = 60             # terminal with remaining work and no write for longer -> "idle?"
RECENT_COMPLETIONS = 20       # number of recent completions used for the per-terminal speed
MAX_GAP_FACTOR = 3.0          # gaps longer than this x median are idle time, not run time
DEFAULT_MINUTES_PER_SEED = {  # used until a terminal has at least 3 completions
    "node1": 14, "node2": 23, "node3": 10, "node4": 21, "node5": 20,
}
WRITE_REPORT = True           # also write revision_R1_runs/progress_status.md
# =============================================================================


def _read_driver_config(path):
    """Return SEED_SETS, SEED_TARGET_BY_ANALYSIS and {analysis: [settings]} from the driver source."""
    tree = ast.parse(open(path, encoding="utf-8").read())
    seed_sets = targets = None
    targets_by_setting = {}
    settings = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        name = getattr(node.targets[0], "id", None)
        if name == "SEED_SETS":
            seed_sets = ast.literal_eval(node.value)
        elif name == "SEED_TARGET_BY_ANALYSIS":
            targets = ast.literal_eval(node.value)
        elif name == "SEED_TARGET_BY_SETTING":
            targets_by_setting = ast.literal_eval(node.value)
        elif name == "ANALYSES" and isinstance(node.value, ast.Dict):
            for k, v in zip(node.value.keys, node.value.values):
                an = ast.literal_eval(k)
                st_names = []
                if isinstance(v, ast.Dict):
                    for kk, vv in zip(v.keys, v.values):
                        if ast.literal_eval(kk) == "settings" and isinstance(vv, ast.Dict):
                            st_names = [ast.literal_eval(x) for x in vv.keys]
                settings[an] = st_names
    if seed_sets is None or targets is None or not settings:
        raise RuntimeError(f"Could not read SEED_SETS / SEED_TARGET_BY_ANALYSIS / ANALYSES from {path}")
    return seed_sets, targets, settings, targets_by_setting


def _target(targets, targets_by_setting, an, st):
    t = targets_by_setting.get((an, st))
    if t is None:
        t = targets.get(an)
    if t is None:
        t = max([v for (a, _), v in targets_by_setting.items() if a == an] or [100])
    return int(t)


def _allocate(seed_sets, target):
    """Same allocation as run_revision_scenarios.allocate_seeds (largest remainder)."""
    names = list(seed_sets)
    sizes = {n: len(seed_sets[n]) for n in names}
    total = sum(sizes.values())
    raw = {n: target * sizes[n] / total for n in names}
    k = {n: int(raw[n]) for n in names}
    for n in sorted(names, key=lambda n: raw[n] - k[n], reverse=True)[: target - sum(k.values())]:
        k[n] += 1
    return {n: list(seed_sets[n][: min(k[n], sizes[n])]) for n in names}


def _last_write(d):
    latest = 0.0
    for root, _, files in os.walk(d):
        for f in files:
            try:
                latest = max(latest, os.path.getmtime(os.path.join(root, f)))
            except OSError:
                pass
    if latest == 0.0:
        try:
            latest = os.path.getmtime(d)
        except OSError:
            pass
    return latest


def _seed_status(d, now):
    """('done', t_done) / ('running', t_last) / ('stalled', t_last) / ('todo', None)."""
    be = os.path.join(d, "breakeven", "breakeven_report_ready.json")
    cea = os.path.join(d, "cea", "cea_report_ready.json")
    het = os.path.join(d, "heterogeneity_summary.csv")
    if os.path.exists(be) and os.path.exists(cea) and os.path.exists(het):
        return "done", os.path.getmtime(be)
    if os.path.isdir(d):
        t = _last_write(d)
        return ("running" if (now - t) / 60.0 <= STALL_MINUTES else "stalled"), t
    return "todo", None


def _minutes_per_seed(done_times, terminal):
    ts = sorted(done_times)[-(RECENT_COMPLETIONS + 1):]
    gaps = [(b - a) / 60.0 for a, b in zip(ts, ts[1:]) if b > a]
    if len(gaps) >= 2:
        med = statistics.median(gaps)
        run_gaps = [g for g in gaps if g <= MAX_GAP_FACTOR * med]
        if run_gaps:
            return statistics.median(run_gaps), f"measured ({len(run_gaps)} intervals)"
    return float(DEFAULT_MINUTES_PER_SEED.get(terminal, 20)), "default"


def _fmt_t(t):
    return time.strftime("%m-%d %H:%M", time.localtime(t)) if t else "-"


def _fmt_dur(minutes):
    h = minutes / 60.0
    return f"{h:.1f} h" if h < 48 else f"{h/24:.1f} d"


def main():
    now = time.time()
    seed_sets, targets, settings, targets_by_setting = _read_driver_config(DRIVER_FILE)
    order = [a for a in RUN_ORDER if a in settings] + [a for a in settings if a not in RUN_ORDER and a != "A0_regression_check"]
    terminals = list(seed_sets)

    out = []
    def p(line=""):
        print(line)
        out.append(line)

    p(f"# Revision R1 progress  ({_fmt_t(now)})")
    p(f"Root: {REVISION_RUNS}")

    # ---- per-seed status ----------------------------------------------------
    status = {}            # (an, st, terminal) -> list of (seed, state, t)
    done_times = {t: [] for t in terminals}
    last_write = {t: 0.0 for t in terminals}
    for an in order:
        for st in settings[an]:
            alloc = _allocate(seed_sets, _target(targets, targets_by_setting, an, st))
            for term in terminals:
                rows = []
                for s in alloc[term]:
                    state, t = _seed_status(os.path.join(REVISION_RUNS, an, st, f"seed_{s}"), now)
                    rows.append((s, state, t))
                    if state == "done":
                        done_times[term].append(t)
                    if t:
                        last_write[term] = max(last_write[term], t)
                status[(an, st, term)] = rows

    # ---- speed and remaining work per terminal ------------------------------
    speed, speed_src = {}, {}
    for term in terminals:
        speed[term], speed_src[term] = _minutes_per_seed(done_times[term], term)

    p("")
    p("## Progress by analysis and setting (done / allocated)")
    for an in order:
        p(f"\n### {an}")
        for st in settings[an]:
            tot_done = sum(1 for term in terminals for _, s_, _ in status[(an, st, term)] if s_ == "done")
            tot_all = sum(len(status[(an, st, term)]) for term in terminals)
            parts = []
            for term in terminals:
                rows = status[(an, st, term)]
                dn = sum(1 for _, s_, _ in rows if s_ == "done")
                extra = [f"run s{s}" for s, s_, _ in rows if s_ == "running"] + \
                        [f"STALLED s{s} (last write {_fmt_t(t)})" for s, s_, t in rows if s_ == "stalled"]
                parts.append(f"{term} {dn}/{len(rows)}" + (f" [{', '.join(extra)}]" if extra else ""))
            mark = "  complete" if tot_done == tot_all else ""
            p(f"- {st}: {tot_done}/{tot_all} (target {_target(targets, targets_by_setting, an, st)}){mark}")
            p("    " + " | ".join(parts))

    # ---- projection ----------------------------------------------------------
    p("")
    p("## Terminals and projected finish (each terminal works through RUN_ORDER without pauses)")
    finish = {}            # terminal -> projected finish time
    an_finish = {an: 0.0 for an in order}
    remaining_total = 0
    for term in terminals:
        cum = 0.0
        rem_term = 0
        for an in order:
            rem = 0.0
            for st in settings[an]:
                for s, s_, t in status[(an, st, term)]:
                    if s_ == "done":
                        continue
                    if s_ == "running":
                        rem += 0.5   # in progress: count as half a seed
                    else:
                        rem += 1.0   # todo or stalled (a stalled seed must be re-run)
            rem_term += rem
            cum += rem * speed[term]
            if rem > 0:
                an_finish[an] = max(an_finish[an], now + cum * 60.0)
        remaining_total += rem_term
        finish[term] = now + cum * 60.0
        running_now = any(s_ == "running" for k, rows in status.items() if k[2] == term for _, s_, _ in rows)
        idle = (rem_term > 0 and not running_now and (now - last_write[term]) / 60.0 > IDLE_MINUTES)
        flag = "  <-- IDLE? remaining work but no recent writes" if idle else ""
        p(f"- {term}: {speed[term]:.0f} min/seed ({speed_src[term]}); remaining {rem_term:.1f} seed-runs = "
          f"{_fmt_dur(cum)}; finish ~{_fmt_t(finish[term])}; last write {_fmt_t(last_write[term])}{flag}")

    stalled = [(k, s, t) for k, rows in status.items() for s, s_, t in rows if s_ == "stalled"]
    if stalled:
        p("")
        p("## Stalled seeds (re-run with ANALYSIS and TERMINAL below; completed seeds are skipped)")
        for (an, st, term), s, t in stalled:
            p(f"- {an} / {st} / seed {s}: ANALYSIS = \"{an}\", TERMINAL = \"{term}\"  (last write {_fmt_t(t)})")

    p("")
    p("## Projected completion by analysis (latest terminal)")
    for an in order:
        p(f"- {an}: {'complete' if an_finish[an] == 0.0 else '~' + _fmt_t(an_finish[an])}")

    if remaining_total > 0:
        rate = sum(1.0 / speed[t] for t in terminals)       # seed-runs per minute, all terminals
        bal = now + remaining_total / rate * 60.0
        slow = max(terminals, key=lambda t: finish[t])
        fast = min(terminals, key=lambda t: finish[t])
        p("")
        p(f"Overall: {remaining_total:.0f} seed-runs left; last terminal ({slow}) ~{_fmt_t(finish[slow])}; "
          f"if perfectly balanced ~{_fmt_t(bal)}.")
        if finish[slow] - bal > 3 * 3600:
            p(f"Rebalancing: when {fast} finishes, run it with TERMINAL = \"{slow}\" and SEED_ORDER = \"reverse\" "
              f"(it works from the end of {slow}'s list while {slow} works from the front; watch the one seed "
              f"where they meet).")

    if WRITE_REPORT:
        path = os.path.join(REVISION_RUNS, "progress_status.md")
        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write("\n".join(out) + "\n")
            print(f"\n(report written to {path})")
        except OSError as e:
            print(f"\n(could not write report: {e})")


if __name__ == "__main__":
    main()
