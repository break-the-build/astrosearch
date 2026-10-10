"""A stalled server must raise, not hang: requests (astroquery, lightkurve)
and http.client (Gaia TAP) both get the installed default timeout."""
import http.client
import socket
import sys
import threading
import time
import unittest
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from vetting import network  # noqa: E402


class SilentServer:
    """Accepts TCP connections and never sends a byte."""

    def __enter__(self):
        self.sock = socket.socket()
        self.sock.bind(("127.0.0.1", 0))
        self.sock.listen(8)
        self.clients = []
        self.running = True
        threading.Thread(target=self._accept, daemon=True).start()
        return self.sock.getsockname()[1]

    def _accept(self):
        self.sock.settimeout(0.2)
        while self.running:
            try:
                self.clients.append(self.sock.accept()[0])
            except OSError:
                continue

    def __exit__(self, *exc):
        self.running = False
        for client in self.clients:
            client.close()
        self.sock.close()


class NetworkTimeoutTests(unittest.TestCase):
    def tearDown(self):
        network.uninstall_default_timeouts()

    def test_requests_without_timeout_fail_fast(self):
        network.install_default_timeouts(0.5)
        with SilentServer() as port:
            start = time.monotonic()
            with self.assertRaises(requests.exceptions.Timeout):
                requests.Session().get(f"http://127.0.0.1:{port}/")
            self.assertLess(time.monotonic() - start, 5)

    def test_explicit_timeout_is_preserved(self):
        network.install_default_timeouts(30)
        with SilentServer() as port:
            start = time.monotonic()
            with self.assertRaises(requests.exceptions.Timeout):
                requests.get(f"http://127.0.0.1:{port}/", timeout=0.3)
            self.assertLess(time.monotonic() - start, 5)

    def test_http_client_uses_socket_default(self):
        network.install_default_timeouts(0.5)
        with SilentServer() as port:
            connection = http.client.HTTPConnection("127.0.0.1", port)
            start = time.monotonic()
            with self.assertRaises((socket.timeout, TimeoutError)):
                connection.request("GET", "/")
                connection.getresponse()
            self.assertLess(time.monotonic() - start, 5)

    def test_install_is_idempotent_and_validates(self):
        network.install_default_timeouts(10)
        wrapped = requests.Session.request
        network.install_default_timeouts(20)
        self.assertIs(requests.Session.request, wrapped)
        with self.assertRaises(ValueError):
            network.install_default_timeouts(0)


if __name__ == "__main__":
    unittest.main()
