"""Ported from cna_maste `python/cnamaste/config/clock.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.

The run's clock: every logger's zero.
"""

import time

start_time = time.time()
"""When `cnamaste.config` was first imported: the zero of every logged elapsed time."""
