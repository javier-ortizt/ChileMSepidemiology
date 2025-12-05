"""
Pesos Estándar WHO
==================

Pesos de la Población Estándar Mundial WHO (2000-2025)
para el cálculo de tasas estandarizadas por edad (ASR).
"""

from __future__ import annotations
from .config import TOP_OPEN


# =============================================================================
# PESOS WHO - GRUPOS QUINQUENALES (5 AÑOS)
# =============================================================================

def who_world_standard_weights_collapsed(top_open: int = TOP_OPEN) -> dict:
    """
    Pesos estándar WHO por millón, colapsando grupos ≥80 años en un solo grupo abierto.
    
    La población estándar WHO (2000-2025) define pesos por millón para cada
    grupo quinquenal. Esta función colapsa los grupos 80-84, 85-89, 90-94,
    95-99 y 100+ en un solo grupo '80+'.
    
    Args:
        top_open: Límite inferior del grupo abierto (default 80)
        
    Returns:
        Dict con grupos de edad como claves y pesos normalizados (suman 1) como valores
    """
    # Pesos originales WHO por millón
    per_million = {
        "0-4": 88569,
        "5-9": 86870,
        "10-14": 85970,
        "15-19": 84670,
        "20-24": 82171,
        "25-29": 79272,
        "30-34": 76073,
        "35-39": 71475,
        "40-44": 65877,
        "45-49": 60379,
        "50-54": 53681,
        "55-59": 45484,
        "60-64": 37187,
        "65-69": 29590,
        "70-74": 22092,
        "75-79": 15195,
        "80-84": 9097,
        "85-89": 4398,
        "90-94": 1500,
        "95-99": 400,
        "100+": 50
    }
    
    # Suma de grupos ≥80 años
    top_sum = (per_million["80-84"] + per_million["85-89"] + 
               per_million["90-94"] + per_million["95-99"] + per_million["100+"])
    
    # Remover grupos originales y agregar grupo colapsado
    keep = {k: v for k, v in per_million.items() 
            if k not in ["80-84", "85-89", "90-94", "95-99", "100+"]}
    keep[f"{top_open}+"] = top_sum
    
    # Normalizar a proporciones (suman 1)
    total = sum(keep.values())
    return {k: v / total for k, v in keep.items()}


# =============================================================================
# PESOS WHO - GRUPOS DECENALES (10 AÑOS)
# =============================================================================

def who_world_standard_weights_10y(top_open: int = TOP_OPEN) -> dict:
    """
    Pesos estándar WHO agregados a grupos decenales (0-9, 10-19, ..., 70-79, 80+).
    
    Suma los pesos quinquenales para obtener grupos de 10 años,
    útil cuando los datos solo están disponibles en ese formato.
    
    Args:
        top_open: Límite inferior del grupo abierto (default 80)
        
    Returns:
        Dict con grupos decenales y pesos normalizados
    """
    # Partir de pesos quinquenales
    w5 = who_world_standard_weights_collapsed(top_open=top_open)
    out = {}

    for ag, w in w5.items():
        if ag.endswith("+"):
            # Grupo abierto se mantiene
            key = f"{top_open}+"
        else:
            lo = int(ag.split("-")[0])
            if lo >= top_open:
                key = f"{top_open}+"
            else:
                # Agrupar en décadas: 0-9, 10-19, etc.
                g10 = (lo // 10) * 10
                hi10 = g10 + 9
                key = f"{g10}-{hi10}"
        
        out[key] = out.get(key, 0.0) + w

    # Re-normalizar
    total = sum(out.values())
    if total > 0:
        out = {k: v / total for k, v in out.items()}
    
    return out

