"""Run a function in a separate process with a wall-clock deadline.

Threads cannot be interrupted, so a stalled computation (for example a
TRICERATOPS run that never converges or waits on a network call) is run in
a spawned child that is terminated at the deadline. The function and its
arguments must be picklable (module-level functions, plain data, numpy).
"""
from __future__ import annotations

import multiprocessing
import queue as queue_module


def _child(fn, args, kwargs, results) -> None:
    try:
        results.put(("ok", fn(*args, **kwargs)))
    except BaseException as error:  # reported to the parent, never raised there blindly
        results.put(("error", f"{type(error).__name__}: {error}"))


def run_with_deadline(fn, args: tuple = (), kwargs: dict | None = None, timeout: float = 3600.0):
    """Return fn(*args, **kwargs) computed in a child process.

    Raises TimeoutError if it does not finish within `timeout` seconds (the
    child is terminated) and RuntimeError carrying the child's error text if
    it raised."""
    if timeout <= 0:
        raise ValueError("timeout must be positive")
    context = multiprocessing.get_context("spawn")
    results = context.Queue()
    process = context.Process(target=_child, args=(fn, args, kwargs or {}, results), daemon=True)
    process.start()
    try:
        status, value = results.get(timeout=timeout)
    except queue_module.Empty:
        raise TimeoutError(f"{getattr(fn, '__name__', fn)} exceeded {timeout:g} s") from None
    finally:
        if process.is_alive():
            process.terminate()
            process.join(5)
            if process.is_alive():
                process.kill()
        process.join(1)
        results.close()
    if status == "error":
        raise RuntimeError(value)
    return value
