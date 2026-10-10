"""Process-wide network timeouts for MAST, Gaia and archive clients.

requests (used by astroquery and lightkurve) waits forever when no timeout
is passed, and http.client (used by the Gaia TAP client) falls back to the
socket default, which is also unlimited. One stalled MAST connection hung
the known-systems benchmark for 45 minutes in poll(). Installing these
defaults turns such stalls into exceptions that callers already handle.
"""
from __future__ import annotations

import functools
import socket

DEFAULT_TIMEOUT_SECONDS = 300.0
_installed: dict = {}


def install_default_timeouts(seconds: float = DEFAULT_TIMEOUT_SECONDS) -> None:
    """Give every requests call and new socket a finite timeout unless the
    caller passed one explicitly. Idempotent; the latest value wins."""
    import requests

    if seconds <= 0:
        raise ValueError("timeout must be positive")
    socket.setdefaulttimeout(seconds)
    _installed["seconds"] = seconds
    if "original" in _installed:
        return
    original = requests.Session.request
    _installed["original"] = original

    @functools.wraps(original)
    def request(self, method, url, **kwargs):
        if kwargs.get("timeout") is None:
            kwargs["timeout"] = _installed["seconds"]
        return original(self, method, url, **kwargs)

    requests.Session.request = request


def uninstall_default_timeouts() -> None:
    """Restore the libraries' defaults (for tests)."""
    import requests

    if "original" in _installed:
        requests.Session.request = _installed.pop("original")
    _installed.pop("seconds", None)
    socket.setdefaulttimeout(None)
