"""
Utilidades Generales
====================

Funciones de normalización, formateo y clasificación.
"""

from __future__ import annotations
import re
import unicodedata
import pandas as pd
import numpy as np

from .config import TOP_OPEN


# =============================================================================
# NORMALIZACIÓN DE TEXTO
# =============================================================================

def normalize(text: str) -> str:
    """
    Normaliza texto: remueve acentos, convierte a minúsculas, limpia espacios.
    
    Args:
        text: Texto a normalizar
        
    Returns:
        Texto normalizado
    """
    if text is None:
        return ""
    text = unicodedata.normalize("NFKD", str(text))
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    text = re.sub(r"\s+", " ", text.lower().strip())
    return text


def find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """
    Busca una columna en el DataFrame por nombre parcial (normalizado).
    
    Args:
        df: DataFrame donde buscar
        candidates: Lista de nombres candidatos (parciales)
        
    Returns:
        Nombre original de la columna encontrada, o None
    """
    normalized_cols = {normalize(c): c for c in df.columns}
    for want in candidates:
        for ncol, orig in normalized_cols.items():
            if want in ncol:
                return orig
    return None


# =============================================================================
# CLASIFICACIÓN DE PREVISIÓN
# =============================================================================

def classify_prevision(value: str) -> str:
    """
    Clasifica el sistema de previsión en FONASA, FFAA o ISAPRE.
    
    Args:
        value: Valor de previsión del registro
        
    Returns:
        'FONASA', 'FFAA' o 'ISAPRE'
    """
    v = normalize(value)
    if "fonasa" in v:
        return "FONASA"
    if "ffaa" in v or "fuerzas armadas" in v:
        return "FFAA"
    return "ISAPRE"


# =============================================================================
# NORMALIZACIÓN DE SEXO
# =============================================================================

def normalize_sex(value: str) -> str:
    """
    Normaliza valores de sexo a 'F' o 'M'.
    
    Args:
        value: Valor de sexo del registro
        
    Returns:
        'F', 'M' o valor original si no reconocido
    """
    v = normalize(value)
    if v in ("f", "fem", "femenino", "mujer"):
        return "F"
    if v in ("m", "masc", "masculino", "hombre"):
        return "M"
    return str(value).strip() if pd.notna(value) else ""


# =============================================================================
# GRUPOS DE EDAD
# =============================================================================

def to_age_group_5y(age: float | int, top_open: int = TOP_OPEN) -> str:
    """
    Convierte edad a grupo quinquenal (0-4, 5-9, ..., 75-79, 80+).
    
    Args:
        age: Edad en años
        top_open: Límite superior abierto (default 80)
        
    Returns:
        Grupo de edad como string (ej: '30-34', '80+')
    """
    try:
        a = int(age)
    except Exception:
        return "Unknown"
    if a < 0:
        return "Unknown"
    if a >= top_open:
        return f"{top_open}+"
    low = (a // 5) * 5
    high = low + 4
    return f"{low}-{high}"


def to_age_group_10y(age: float | int, top_open: int = TOP_OPEN) -> str:
    """
    Convierte edad a grupo decenal (0-9, 10-19, ..., 70-79, 80+).
    
    Args:
        age: Edad en años
        top_open: Límite superior abierto (default 80)
        
    Returns:
        Grupo de edad como string (ej: '30-39', '80+')
    """
    try:
        a = int(age)
    except Exception:
        return "Unknown"
    if a < 0:
        return "Unknown"
    if a >= top_open:
        return f"{top_open}+"
    low = (a // 10) * 10
    high = low + 9
    return f"{low}-{high}"


def age_group_sortkey(g: str) -> tuple:
    """
    Clave de ordenamiento para grupos de edad.
    Orden lógico: 0-4, 5-9, …, 75-79, 80+.
    
    Args:
        g: Grupo de edad como string
        
    Returns:
        Tupla para ordenamiento
    """
    if not isinstance(g, str):
        return (9999, 0)
    g = g.strip()
    m = re.match(r"^(\d+)\s*-\s*(\d+)$", g)
    if m:
        return (int(m.group(1)), 0)
    m2 = re.match(r"^(\d+)\+$", g)
    if m2:
        return (int(m2.group(1)), 1)
    return (9999, 1)


def sort_age_index(idx):
    """Ordena un índice de grupos de edad."""
    return sorted(list(idx), key=age_group_sortkey)


# =============================================================================
# FORMATEADORES PARA EJES DE GRÁFICOS
# =============================================================================

def fmt_lon(x, pos=None):
    """Formatea longitud para ejes de mapas."""
    return f"{abs(x):.0f}°{'W' if x < 0 else 'E'}"


def fmt_lat(y, pos=None):
    """Formatea latitud para ejes de mapas."""
    return f"{abs(y):.0f}°{'S' if y < 0 else 'N'}"

