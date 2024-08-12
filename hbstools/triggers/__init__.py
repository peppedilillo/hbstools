from typing import Any


class TriggerAlgorithm:
    """Meant to be subclassed. Implements utilities"""

    def _asdict(self) -> dict[str, Any]:
        """From Ramalho's Fluent Python 2ed."""
        return {
            name: getattr(self, name) for name, attr in self.__class__.__dict__.items()
        }

    def __repr__(self) -> str:
        """From Ramalho's Fluent Python 2ed."""
        kwargs = ", ".join(f"{key}={value!r}" for key, value in self._asdict().items())
        return f"{self.__class__.__name__}({kwargs})"


def _library_path(libname: str) -> str:
    """Gets the path to the C shared library, with appropriate extension.
    Assumes the library and this file to be installed in the same directory"""
    from pathlib import Path
    import sys

    if sys.platform.startswith("win32"):
        suffix = ".dll"
    elif sys.platform.startswith("linux"):
        suffix = ".so"
    elif sys.platform.startswith("darwin"):
        suffix = ".dylib"
    else:
        raise OSError("System not supported")

    return str(Path(__file__).parent / f"{libname}{suffix}")


_LIBCFOCUS = _library_path(".sharedlibs/lib-pfocus")
