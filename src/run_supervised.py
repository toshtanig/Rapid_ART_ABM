#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Start run_revision_scenarios.py and start it again if it ends abnormally.

Run this file instead of run_revision_scenarios.py (TERMINAL is still set in the CONFIG block of
run_revision_scenarios.py; nothing else changes). The driver runs as a separate process with the Python
interpreter that runs this file, in the same terminal, so its lines appear as usual.

- If the driver is ended by a crash of the interpreter or is killed by the system (segmentation fault,
  abort, out of memory), it is started again after RESTART_WAIT_SECONDS. The seed that was running is run
  again from its start; completed seeds are skipped as always.
- If the driver ends normally (queue complete, or it stopped itself after errors), or is stopped by the user
  (Ctrl+C, or pkill -f run_revision_scenarios), it is not started again.
- If the driver ends abnormally MAX_ABNORMAL_ENDS times within WINDOW_MINUTES, this script stops: a crash
  that repeats that quickly needs a look before any further start.

Every start and end is written to _logs/supervisor_<TERMINAL>.log; an abnormal end is also written to
_logs/run_<TERMINAL>.log as an ERROR line, so that check_progress.py lists it.

Run as-is in VS Code; all settings are in the CONFIG block, no command-line arguments.
"""
import os
import re
import sys
import time
import signal
import socket
import datetime
import subprocess

# =============================================================================
# CONFIG
# =============================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _paths import REVISION_RUNS

DRIVER = os.path.join(HERE, "run_revision_scenarios.py")
LOG_DIR = os.path.join(REVISION_RUNS, "_logs")
RESTART_WAIT_SECONDS = 60        # pause before the driver is started again
MAX_ABNORMAL_ENDS = 3            # stop if the driver ends abnormally this many times ...
WINDOW_MINUTES = 120             # ... within this time
# The driver is started again only if it was ended by one of these signals
RESTART_AFTER_SIGNALS = ["SIGSEGV", "SIGBUS", "SIGILL", "SIGFPE", "SIGABRT", "SIGKILL"]
# =============================================================================


def terminal_name():
    """TERMINAL as set in the CONFIG block of the driver (read from the file; the driver is not imported)."""
    with open(DRIVER, encoding="utf-8") as f:
        m = re.search(r'^TERMINAL\s*=\s*["\'](.*?)["\']', f.read(), flags=re.M)
    return m.group(1) if m else ""


TERMINAL = terminal_name()


def log(msg, run_log=False):
    line = f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} [{TERMINAL}] {msg}"
    print(line, flush=True)
    names = [f"supervisor_{TERMINAL}.log"] + ([f"run_{TERMINAL}.log"] if run_log else [])
    for name in names:
        try:
            os.makedirs(LOG_DIR, exist_ok=True)
            with open(os.path.join(LOG_DIR, name), "a", encoding="utf-8") as f:
                f.write(line + "\n")
        except OSError:
            pass


def signal_name(returncode):
    """Name of the signal that ended the process, or None if it ended by itself."""
    if returncode is None or returncode >= 0:
        return None
    try:
        return signal.Signals(-returncode).name
    except ValueError:
        return f"signal {-returncode}"


def main():
    if not TERMINAL:
        raise SystemExit('Set TERMINAL in the CONFIG block of run_revision_scenarios.py to the name of this machine.')
    log(f"supervisor started on {socket.gethostname()} (process {os.getpid()}); driver: {DRIVER}")
    abnormal = []                                   # times of the abnormal ends
    starts = 0
    while True:
        starts += 1
        t0 = time.time()
        proc = subprocess.Popen([sys.executable, DRIVER], cwd=HERE)
        log(f"driver start {starts} (process {proc.pid})")
        try:
            rc = proc.wait()
        except KeyboardInterrupt:
            # Ctrl+C reaches the driver as well; it withdraws its claim and ends by itself
            try:
                rc = proc.wait(timeout=120)
            except (subprocess.TimeoutExpired, KeyboardInterrupt):
                proc.kill()
                rc = proc.wait()
            log(f"stopped by the user (driver ended with code {rc}); not started again")
            return
        minutes = (time.time() - t0) / 60.0
        sig = signal_name(rc)
        if sig is None:
            log(f"driver ended by itself after {minutes:.1f} min (exit code {rc}); not started again")
            return
        if sig not in RESTART_AFTER_SIGNALS:
            log(f"driver stopped by {sig} after {minutes:.1f} min; not started again")
            return
        abnormal.append(time.time())
        recent = [t for t in abnormal if time.time() - t <= WINDOW_MINUTES * 60.0]
        if len(recent) >= MAX_ABNORMAL_ENDS:
            log(f"driver: ERROR ended abnormally ({sig}) after {minutes:.1f} min; {len(recent)} abnormal ends within "
                f"{WINDOW_MINUTES} min (MAX_ABNORMAL_ENDS = {MAX_ABNORMAL_ENDS}): not started again. "
                f"See fault_{TERMINAL}.log and memory_{TERMINAL}.log in {LOG_DIR}", run_log=True)
            return
        log(f"driver: ERROR ended abnormally ({sig}) after {minutes:.1f} min; started again in "
            f"{RESTART_WAIT_SECONDS} s (see fault_{TERMINAL}.log and memory_{TERMINAL}.log)", run_log=True)
        try:
            time.sleep(RESTART_WAIT_SECONDS)
        except KeyboardInterrupt:
            log("stopped by the user while waiting; not started again")
            return


if __name__ == "__main__":
    main()
