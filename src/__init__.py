"""KES v7 — paket za analitiko, GA-optimizacijo in FEMM verifikacijo
6-polnega sinhronskega motorja z vzbujalnim navitjem.

Ob uvozu paketa preklopimo standardna izhoda na UTF-8. Brez tega na
Windows konzoli (privzeto cp1250) `print` z znaki kot je `η`, `δ`, `³`
sproži `UnicodeEncodeError` in skripta pade še pred zapisom izhodnih
datotek. `reconfigure` je na voljo od Pythona 3.7 naprej; če iz kakršnega
koli razloga ni (npr. preusmerjen tok brez `reconfigure`), napako tiho
požremo, da uvoz paketa nikoli ne pade.
"""

from __future__ import annotations

import sys


def _force_utf8_stdio() -> None:
    """Zagotovi UTF-8 izhod, da unicode znaki ne sprožijo izjeme."""
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            try:
                reconfigure(encoding="utf-8")
            except (ValueError, OSError):
                pass


_force_utf8_stdio()
