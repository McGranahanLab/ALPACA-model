try:
    from importlib.metadata import version  # type: ignore
except Exception:
    try:
        from importlib_metadata import version  # type: ignore
    except Exception:
        def version(pkg_name="alpaca"):
            return "unknown"

import os as _os
import tempfile as _tempfile


def _is_writable_dir(path):
    try:
        _os.makedirs(path, exist_ok=True)
        return _os.access(path, _os.W_OK)
    except OSError:
        return False


def _configure_writable_caches():
    """Point matplotlib/fontconfig caches at a writable dir (read-only $HOME in containers).

    Uses ALPACA_CACHE_DIR if set; otherwise only overrides when the default location is unwritable.
    """
    explicit = _os.environ.get("ALPACA_CACHE_DIR")
    if explicit:
        base = _os.path.abspath(_os.path.expanduser(explicit))
    else:
        uid = _os.getuid() if hasattr(_os, "getuid") else "user"
        base = _os.path.join(_tempfile.gettempdir(), f"alpaca-cache-{uid}")

    home = _os.path.expanduser("~")
    xdg_cache = _os.environ.get("XDG_CACHE_HOME") or _os.path.join(home, ".cache")
    mpl_dir = _os.path.join(
        _os.environ.get("XDG_CONFIG_HOME") or _os.path.join(home, ".config"), "matplotlib"
    )

    if not _os.environ.get("MPLCONFIGDIR") and (explicit or not _is_writable_dir(mpl_dir)):
        target = _os.path.join(base, "matplotlib")
        if _is_writable_dir(target):
            _os.environ["MPLCONFIGDIR"] = target
    if explicit or not _is_writable_dir(_os.path.join(xdg_cache, "fontconfig")):
        target = _os.path.join(base, "xdg")
        if _is_writable_dir(target):
            _os.environ["XDG_CACHE_HOME"] = target


_configure_writable_caches()

__all__ = ["version"]
