def _library_path(libname: str) -> str:
    """Gets the path to the C shared library, with appropriate extension.
    Assumes the library and this file to be installed in the same directory"""
    from pathlib import Path
    import sys

    if sys.platform.startswith('win32'):
        suffix = ".dll"
    elif sys.platform.startswith('linux'):
        suffix = ".so"
    elif sys.platform.startswith('darwin'):
        suffix = ".dylib"
    else:
        raise OSError("System not supported")

    return str(Path(__file__).parent / f"{libname}{suffix}")


_LIBCFOCUS = _library_path(".sharedlibs/lib-pfocus")
