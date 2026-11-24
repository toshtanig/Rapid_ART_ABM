#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Progress monitor for the queue of run_revision_scenarios.py, with a remaining-time estimate.

Reads the output folders under outputs/revision_R2_runs/ (synchronised between the terminals) and reports, for
every queue item, how many seeds are done, which seeds are running on which terminal, and which seeds were
started and abandoned. The speed of each terminal is taken from the run times recorded with its completed
seeds (run_info.json), and the expected finish time of every queue item is projected from the combined
speed of the terminals that are currently active.

The queue, the seed list and the terminal table are read from run_revision_scenarios.py without importing
it (so the model and numba are not loaded): this script needs only the Python standard library and runs
from VS Code on any terminal, at any time. Run as-is; all options are in CONFIG.
"""

import os
import sys
import ast
import glob
import json
import time
import statistics

# =============================================================================
# CONFIG
# =============================================================================
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import REVISION_RUNS

DRIVER_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_revision_scenarios.py")

IDLE_MINUTES = 90             # terminal without a completion or a new claim for longer than this -> not active
RECENT_COMPLETIONS = 20       # number of recent completions used for the speed of a terminal
THREE_ARM_FACTOR = 1.15       # run time of a three-scenario seed relative to a two-scenario seed, until measured
WRITE_REPORT = True           # also write outputs/revision_R2_runs/progress_status.md
# =============================================================================


def _read_driver_config(path):
    """MAIN_QUEUE, PUBLICATION_SEEDS, MINUTES_PER_SEED_RUN, CLAIM_TIMEOUT_MINUTES and the number of
    scenario slots of every setting, from the driver source."""
    tree = ast.parse(open(path, encoding="utf-8").read())
    cfg = {}
    names = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        name = getattr(node.targets[0], "id", None)
        if name in ("MAIN_QUEUE", "PUBLICATION_SEEDS", "MINUTES_PER_SEED_RUN", "CLAIM_TIMEOUT_MINUTES",
                    "_THREE_ARMS", "_TWO_ARMS"):
            cfg[name] = ast.literal_eval(node.value)
        elif name == "ANALYSES" and isinstance(node.value, ast.Dict):
            for k, v in zip(node.value.keys, node.value.values):
                an = ast.literal_eval(k)
                for kk, vv in zip(v.keys, v.values):
                    if ast.literal_eval(kk) != "settings":
                        continue
                    for sk, sv in zip(vv.keys, vv.values):
                        for fk, fv in zip(sv.keys, sv.values):
                            if ast.literal_eval(fk) == "scenarios":
                                names[(an, ast.literal_eval(sk))] = fv
    arms = {}
    for key, node in names.items():
        if isinstance(node, ast.Name):
            arms[key] = len(cfg.get(node.id, [0, 0]))
        else:
            arms[key] = len(ast.literal_eval(node))
    missing = [k for k in ("MAIN_QUEUE", "PUBLICATION_SEEDS", "MINUTES_PER_SEED_RUN", "CLAIM_TIMEOUT_MINUTES") if k not in cfg]
    if missing:
        raise RuntimeError(f"Could not read {missing} from {path}")
    return cfg, arms


def _owner_map(seeds, minutes):
    """Same rule as run_revision_scenarios.build_owner_map."""
    minutes = {t: float(m) for t, m in minutes.items()}
    count = {t: 0 for t in minutes}
    owner = {}
    for s in seeds:
        t = min(minutes, key=lambda k: ((count[k] + 1) * minutes[k], k))
        count[t] += 1
        owner[s] = t
    return owner


def _seed_status(d, now, timeout_min):
    """dict(state = done / running / abandoned / todo, terminal, t, minutes, md5)."""
    done = (os.path.exists(os.path.join(d, "breakeven", "breakeven_report_ready.json"))
            and os.path.exists(os.path.join(d, "cea", "cea_report_ready.json"))
            and os.path.exists(os.path.join(d, "heterogeneity_summary.csv")))
    if done:
        info = {}
        try:
            with open(os.path.join(d, "run_info.json"), encoding="utf-8") as f:
                info = json.load(f)
        except Exception:
            pass
        t = info.get("finished_epoch")
        if t is None:
            try:
                t = os.path.getmtime(os.path.join(d, "breakeven", "breakeven_report_ready.json"))
            except OSError:
                t = None
        return {"state": "done", "terminal": info.get("terminal"), "t": t, "minutes": info.get("minutes"),
                "md5": info.get("model_file_md5")}
    claims = []
    for path in glob.glob(os.path.join(d, "_claim_*.json")):
        term = os.path.basename(path)[len("_claim_"):-len(".json")]
        try:
            with open(path, encoding="utf-8") as f:
                epoch = float(json.load(f)["epoch"])
        except Exception:
            try:
                epoch = os.path.getmtime(path)
            except OSError:
                continue
        claims.append((epoch, term))
    if not claims:
        return {"state": "todo", "terminal": None, "t": None, "minutes": None, "md5": None}
    epoch, term = max(claims)
    fresh = (now - epoch) / 60.0 < timeout_min
    return {"state": "running" if fresh else "abandoned", "terminal": term, "t": epoch, "minutes": None, "md5": None}


def _fmt_t(t):
    return time.strftime("%m-%d %H:%M", time.localtime(t)) if t else "-"


def _fmt_dur(minutes):
    h = minutes / 60.0
    return f"{h:.1f} h" if h < 48 else f"{h / 24:.1f} d"


def main():
    now = time.time()
    cfg, arms = _read_driver_config(DRIVER_FILE)
    queue = [tuple(q) for q in cfg["MAIN_QUEUE"]]
    seeds_all = cfg["PUBLICATION_SEEDS"]
    table = cfg["MINUTES_PER_SEED_RUN"]
    timeout_min = float(cfg["CLAIM_TIMEOUT_MINUTES"])
    terminals = list(table)
    owner = _owner_map(seeds_all, table)

    out = []
    def p(line=""):
        print(line)
        out.append(line)

    p(f"# Progress of the R2 runs  ({_fmt_t(now)})")
    p(f"Root: {REVISION_RUNS}")

    # ---- status of every seed of every (analysis, setting) -----------------------
    target = {}
    for an, st, n in queue:
        target[(an, st)] = max(int(n), target.get((an, st), 0))
    status = {}
    for (an, st), n in target.items():
        for s in seeds_all[:n]:
            status[(an, st, s)] = _seed_status(os.path.join(REVISION_RUNS, an, st, f"seed_{s}"), now, timeout_min)

    # ---- speed and activity of the terminals ----------------------------------------
    runs = {t: [] for t in terminals}          # (finished, minutes, arms)
    last_activity = {t: 0.0 for t in terminals}
    running_now = {t: [] for t in terminals}
    md5s = {}
    for (an, st, s), r in status.items():
        t = r["terminal"]
        if r["state"] == "done":
            md5s[r["md5"]] = md5s.get(r["md5"], 0) + 1
            if t in runs and r["minutes"] is not None and r["t"] is not None:
                runs[t].append((r["t"], float(r["minutes"]), arms.get((an, st), 2)))
        if t in last_activity and r["t"]:
            last_activity[t] = max(last_activity[t], r["t"])
        if r["state"] == "running" and t in running_now:
            running_now[t].append(f"{an}/{st} seed {s} (since {_fmt_t(r['t'])})")
    ratio3 = []
    speed, speed_src = {}, {}
    for t in terminals:
        recent = sorted(runs[t])[-RECENT_COMPLETIONS:]
        two = [m for _, m, a in recent if a == 2]
        three = [m for _, m, a in recent if a == 3]
        if two and three and statistics.median(two) > 0:
            ratio3.append(statistics.median(three) / statistics.median(two))
    factor3 = statistics.median(ratio3) if ratio3 else THREE_ARM_FACTOR
    for t in terminals:
        recent = sorted(runs[t])[-RECENT_COMPLETIONS:]
        eq = [m / (factor3 if a == 3 else 1.0) for _, m, a in recent]      # two-scenario equivalents
        if len(eq) >= 2 and statistics.median(eq) > 0:
            speed[t], speed_src[t] = statistics.median(eq), f"measured, {len(eq)} runs"
        else:
            speed[t], speed_src[t] = float(table[t]), "table value"
    started = any(last_activity[t] > 0 for t in terminals)
    active = {t: (not started) or (now - last_activity[t]) / 60.0 <= IDLE_MINUTES for t in terminals}

    # ---- progress by queue item -------------------------------------------------------
    p("")
    p("## Queue (done / seeds; R = running, A = abandoned)")
    counted = set()
    remaining_by_item = []
    for i, (an, st, n) in enumerate(queue, 1):
        keys = [(an, st, s) for s in seeds_all[: int(n)]]
        states = [status[k]["state"] for k in keys]
        n_done = states.count("done")
        new_keys = [k for k in keys if k not in counted]
        counted.update(new_keys)
        rem = sum((0.5 if status[k]["state"] == "running" else 1.0) for k in new_keys if status[k]["state"] != "done")
        weight = factor3 if arms.get((an, st), 2) == 3 else 1.0
        remaining_by_item.append(rem * weight)
        extra = []
        if states.count("running"):
            extra.append(f"R {states.count('running')}")
        if states.count("abandoned"):
            extra.append(f"A {states.count('abandoned')}")
        mark = "  complete" if n_done == len(keys) else ""
        p(f"{i:2d}. {an} / {st}: {n_done}/{len(keys)}" + (f"  [{', '.join(extra)}]" if extra else "") + mark)

    # ---- terminals ----------------------------------------------------------------------
    p("")
    p("## Terminals")
    for t in terminals:
        n_done = len(runs[t])
        n_taken = sum(1 for (an, st, s), r in status.items() if r["state"] == "done" and r["terminal"] == t and owner[s] != t)
        if last_activity[t] == 0:
            flag = "  (not started)"
        else:
            flag = "" if active[t] else f"  <-- NOT ACTIVE (no completion or new seed for more than {IDLE_MINUTES} min)"
        p(f"- {t}: {n_done} seed-runs done ({n_taken} taken over from other terminals); "
          f"{speed[t]:.0f} min per two-scenario seed ({speed_src[t]}); last activity {_fmt_t(last_activity[t])}{flag}")
        for r in running_now[t]:
            p(f"    running: {r}")

    abandoned = [(k, r) for k, r in status.items() if r["state"] == "abandoned"]
    if abandoned:
        p("")
        p("## Abandoned seeds (started, not finished, and the claim has timed out; any running terminal takes them over)")
        for (an, st, s), r in abandoned:
            p(f"- {an} / {st} / seed {s}: started on {r['terminal']} at {_fmt_t(r['t'])}")

    # ---- projection -----------------------------------------------------------------------
    p("")
    total_rem = sum(remaining_by_item)
    rate = sum(1.0 / speed[t] for t in terminals if active[t])      # two-scenario seeds per minute
    if total_rem == 0:
        p("## All queue items are complete.")
    elif rate <= 0:
        p(f"## {total_rem:.0f} seed-runs left (two-scenario equivalents), but no terminal is active.")
    else:
        p(f"## Projected completion (active terminals: {', '.join(t for t in terminals if active[t])}; "
          f"three-scenario seeds count {factor3:.2f})")
        cum = 0.0
        for i, ((an, st, n), rem) in enumerate(zip(queue, remaining_by_item), 1):
            cum += rem
            if rem > 0:
                p(f"{i:2d}. {an} / {st} ({n}): ~{_fmt_t(now + cum / rate * 60.0)}")
        p(f"Overall: {total_rem:.0f} seed-runs left (two-scenario equivalents), about {_fmt_dur(total_rem / rate)}; "
          f"finish ~{_fmt_t(now + total_rem / rate * 60.0)}.")

    # ---- consistency checks -----------------------------------------------------------------
    known = {k: v for k, v in md5s.items() if k}
    if len(known) > 1:
        p("")
        p("## WARNING: completed seeds were run with different versions of the model file")
        for k, v in sorted(known.items(), key=lambda kv: -kv[1]):
            p(f"- md5 {k}: {v} seeds")
    errors = []
    for path in sorted(glob.glob(os.path.join(REVISION_RUNS, "_logs", "run_*.log"))):
        try:
            with open(path, encoding="utf-8", errors="replace") as f:
                errors += [line.strip() for line in f if ": ERROR " in line or "] STOPPED: " in line]
        except OSError:
            pass
    if errors:
        p("")
        p(f"## Errors and stops recorded in _logs ({len(errors)}; the last five are shown)")
        for line in errors[-5:]:
            p(f"- {line[:300]}")

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
