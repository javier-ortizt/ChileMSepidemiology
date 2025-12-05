"""
Configuración y Constantes
==========================

Parámetros globales del análisis de EM.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from cycler import cycler
import warnings

# Suprimir advertencias no críticas
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=r".*GroupBy\.apply operated on the grouping columns.*")

# =============================================================================
# ARCHIVOS DE ENTRADA
# =============================================================================
EXCEL_FILE = "BBDD EM.xlsx"                                    # Base de datos principal de casos
INE_POP_FILE = "estimaciones-y-proyecciones-2002-2035-comunas.xlsx"  # Proyecciones INE
PREV_REGION_FILE = "Previsionporregion.xlsx"                   # Población por previsión y región
BEN_SYS_FILE = "beneficiariosTramoEtario.xlsx"                 # Beneficiarios por tramo etario
SHAPEFILE_PATH = "Regional.shp"                                # Shapefile de regiones de Chile

# =============================================================================
# PARÁMETROS TEMPORALES
# =============================================================================
MIN_YEAR = 2010          # Año inicio para PREVALENCIA
INC_START_YEAR = 2010    # Año base para INCIDENCIA (2010 es baseline; reporte desde 2011+)
MAX_YEAR = 2023          # Año final del análisis

# =============================================================================
# PARÁMETROS DE GRUPOS DE EDAD
# =============================================================================
TOP_OPEN = 80            # Grupo abierto superior (80+)

# =============================================================================
# NOMBRES DE ISAPRES PARA GRÁFICOS DETALLADOS
# =============================================================================
ISAPRE_NAMES = [
    "Banmédica",
    "Colmena Golden Cross",
    "Consalud",
    "Cruz Blanca",
    "Isalud",
    "Isapre Fundación",
    "Nueva Masvida",
    "Vida Tres",
]

# =============================================================================
# EXTENSIÓN GEOGRÁFICA DE CHILE CONTINENTAL
# =============================================================================
CHILE_CONTINENTAL_EXTENT = (-74.2, -66.4, -56.0, -17.5)  # (lon_min, lon_max, lat_min, lat_max)

# =============================================================================
# CONFIGURACIÓN DE MATPLOTLIB
# =============================================================================
plt.rcParams.update({
    "figure.dpi": 200,
    "savefig.dpi": 200,
    "axes.prop_cycle": cycler(color=[
        "#3b82f6",  # Azul
        "#ef4444",  # Rojo
        "#10b981",  # Verde
        "#f59e0b",  # Ámbar
        "#8b5cf6",  # Violeta
        "#06b6d4",  # Cian
        "#f97316",  # Naranja
        "#22c55e",  # Verde claro
        "#e11d48",  # Rosa
    ]),
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
})

# =============================================================================
# DIRECTORIOS DE SALIDA
# =============================================================================
def get_output_dirs():
    """Crea y retorna los directorios de salida."""
    outdir = Path("salidas_bdem")
    outdir.mkdir(exist_ok=True)
    gdir = outdir / "graficos"
    gdir.mkdir(exist_ok=True)
    return outdir, gdir

