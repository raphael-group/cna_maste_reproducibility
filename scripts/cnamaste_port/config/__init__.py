"""Ported from cna_maste `python/cnamaste/config/__init__.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim, with their imports pointed at `scripts.cnamaste_port`.

The run's configuration: `YAMLConfig` (`cnamaste.config.yaml_config`), its defaults and choices
(`cnamaste.config.defaults`), and the one configuration a run holds.
"""

from scripts.cnamaste_port.config.clock import start_time

__all__ = ["start_time"]
