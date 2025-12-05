"""
Paquete de Análisis Epidemiológico de Esclerosis Múltiple (EM) - Chile
======================================================================

Este paquete contiene herramientas para el análisis de:
- Incidencia y prevalencia de EM
- Tasas estandarizadas WHO (ASR)
- Análisis por región, sexo, edad y sistema previsional
- Correlación con latitud
- Visualizaciones y mapas coropléticos

Estructura del paquete:
-----------------------
- config.py          : Configuraciones y constantes
- utils.py           : Utilidades generales (normalización, formateo)
- io_data.py         : Carga de datos (Excel, shapefiles)
- who_standards.py   : Pesos estándar WHO para estandarización
- incidence.py       : Cálculos de incidencia
- prevalence.py      : Cálculos de prevalencia
- statistics.py      : Análisis estadísticos (Poisson, Pearson, CI)
- geography.py       : Helpers geográficos y de región
- plotting.py        : Visualizaciones básicas
- maps.py            : Mapas coropléticos
- insurance.py       : Análisis por sistema previsional (ISAPRE/FONASA)
- main.py            : Orquestación del análisis completo
"""

from .config import *
from .main import main

