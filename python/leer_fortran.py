# -*- coding: utf-8 -*-
"""Lector de .dat de Fortran: exponentes de 3 digitos sin 'E' (0.12-102) y cabecera '# clave = valor'."""
import re
import numpy as np


def leer_dat(path):
    """Devuelve (filas como array, dict con las claves numericas de la cabecera)."""
    meta, filas = {}, []
    for linea in open(path, encoding="utf-8"):
        linea = linea.strip()
        if not linea:
            continue
        if linea.startswith("#"):
            if "=" in linea:
                k, v = linea.lstrip("#").split("=", 1)
                try:
                    meta[k.strip()] = float(v)
                except ValueError:
                    pass
            continue
        linea = re.sub(r"(\d)([-+])(\d{2,3})\b", r"\1E\2\3", linea)
        filas.append([float(x) for x in linea.split()])
    return np.array(filas), meta
