"""Test Suite para GFF Statistics Parser

Este módulo contiene pruebas unitarias (usando pytest) para validar que el módulo
gff_stats.py funciona correctamente. Las pruebas cubren:

1. Cálculo correcto de estadísticas desde archivos GFF
2. Filtrado por tipo de feature
3. Generación correcta de archivos JSON
4. Validación de valores calculados (promedios, distribuciones, etc.)

Estructura:
- load_module(): función auxiliar para importar dinámicamente el módulo bajo prueba
- test_*: funciones de prueba (pytest detecta automáticamente funciones que comienzan con test_)

Datos de prueba:
- data/sample.gff: archivo GFF con 6 features para pruebas
  Contenido: 2 genes, 3 CDS, 1 mRNA (combinación de + y -)
"""
import json
import importlib.util
from pathlib import Path


# CONFIGURACIÓN DE RUTAS
# Obtener la ruta raíz del proyecto (directorio padre de tests/)
ROOT = Path(__file__).resolve().parents[1]
# Ruta al módulo bajo prueba
SRC_PATH = ROOT / "src" / "gff_stats.py"


def load_module():
    """Carga dinámicamente el módulo gff_stats.py usando importlib.
    
    Esta función es necesaria porque:
    - pytest ejecuta tests desde un directorio diferente
    - Necesitamos cargar el módulo desde una ruta específica
    - Permite recargar el módulo en cada test (evita efectos secundarios entre tests)
    
    Internals (importlib):
        1. spec_from_file_location(): crea una especificación del módulo desde archivo
        2. module_from_spec(): crea un módulo vacío basado en la especificación
        3. exec_module(): ejecuta el código del módulo en el namespace del módulo
    
    Returns:
        module: Módulo gff_stats cargado y listo para usar
        
    Ejemplo:
        >>> mod = load_module()
        >>> stats = mod.compute_stats_from_gff("datos.gff")
    """
    # Crear especificación del módulo desde la ruta del archivo
    spec = importlib.util.spec_from_file_location("gff_stats", str(SRC_PATH))
    # Crear módulo vacío basado en la especificación
    mod = importlib.util.module_from_spec(spec)
    # Ejecutar el código del módulo en el namespace del módulo
    spec.loader.exec_module(mod)
    return mod


def test_compute_stats_from_sample():
    """Test: Verificar cálculo correcto de estadísticas sin filtro.
    
    Objetivo:
        Validar que compute_stats_from_gff() procesa correctamente un archivo GFF
        y calcula todas las estadísticas esperadas.
    
    Datos de entrada:
        - data/sample.gff: archivo con 6 features
          * 2 genes (longitudes 900 cada uno, 1 +, 1 -)
          * 3 CDS (longitudes 101, 51, 101, mezcla de + y -)
          * 1 mRNA (longitud 150, +)
    
    Validaciones:
        1. total_features: debe ser 6 (todas las líneas no comentadas)
        2. by_type: debe contar correctamente por tipo
        3. avg_length: debe calcular promedios correctos (redondeado a 1 decimal)
        4. strand_distribution: debe calcular porcentajes de + y -
    
    Cálculos manuales (verificación):
        - gene: (900 + 900) / 2 = 900.0
        - CDS: (101 + 51 + 101) / 3 = 253 / 3 = 84.333... ≈ 84.3
        - mRNA: 150 / 1 = 150.0
        - Strands: + = 4 features, - = 2 features
          * % +: 4 / 6 * 100 = 66.67% ≈ 67 (redondeado)
          * % -: 2 / 6 * 100 = 33.33% ≈ 33 (redondeado)
    """
    # Cargar módulo dinámicamente
    mod = load_module()
    # Obtener ruta del archivo de prueba
    sample = ROOT / "data" / "sample.gff"
    # Llamar función bajo prueba (sin filtro)
    stats = mod.compute_stats_from_gff(str(sample))

    # ASERCIONES: verificar resultados
    # Verificar total de features
    assert stats["total_features"] == 6
    assert stats["by_type"] == {"gene": 2, "CDS": 3, "mRNA": 1}
    # averages: gene lengths 900 and 900 -> 900.0; CDS lengths 101,51,101 -> 84.3
    assert stats["avg_length"]["gene"] == 900.0
    assert stats["avg_length"]["CDS"] == 84.3
    assert stats["avg_length"]["mRNA"] == 150.0
    # strand distribution: + = 4, - = 2 -> 67 and 33 (rounded)
    assert stats["strand_distribution"]["+"] in (66, 67, 67)
    assert stats["strand_distribution"]["-"] in (33, 34, 33)


def test_compute_stats_with_filter():
    """Test: Verificar filtrado por tipo de feature.
    
    Objetivo:
        Validar que compute_stats_from_gff() filtra correctamente cuando se especifica
        un tipo de feature, y solo calcula estadísticas para ese tipo.
    
    Casos de prueba:
        1. Filtro por CDS: debe retornar solo 3 features (los 3 CDS)
        2. Filtro por gene: debe retornar solo 2 features (los 2 genes)
    
    Validaciones:
        - total_features: debe ser solo el conteo del tipo filtrado
        - by_type: debe contener solo el tipo filtrado
        - avg_length: debe contener solo promedio del tipo filtrado
        - filter_type: debe estar presente en el resultado
        - strand_distribution: se recalcula solo para features filtrados
    """
    # Cargar módulo dinámicamente
    mod = load_module()
    # Obtener ruta del archivo de prueba
    sample = ROOT / "data" / "sample.gff"
    
    # TEST 1: Filtrar por tipo CDS
    stats = mod.compute_stats_from_gff(str(sample), filter_type="CDS")
    # Verificaciones para CDS
    assert stats["total_features"] == 3
    assert stats["by_type"] == {"CDS": 3}
    assert stats["avg_length"]["CDS"] == 84.3
    assert stats["filter_type"] == "CDS"
    
    # TEST 2: Filtrar por tipo gene
    stats = mod.compute_stats_from_gff(str(sample), filter_type="gene")
    # Verificaciones para gene
    assert stats["total_features"] == 2
    assert stats["by_type"] == {"gene": 2}
    assert stats["avg_length"]["gene"] == 900.0
    assert stats["filter_type"] == "gene"


def test_cli_writes_json(tmp_path):
    """Test: Verificar que CLI genera correctamente archivo JSON.
    
    Objetivo:
        Validar que la función cli() procesa argumentos correctamente,
        llama a compute_stats_from_gff(), y escribe JSON válido al archivo.
    
    Parámetros:
        tmp_path: fixture de pytest que proporciona directorio temporal único
                  Útil para evitar conflictos de archivos en tests paralelos
    
    Flujo:
        1. Cargar módulo
        2. Preparar rutas: archivo de entrada (sample.gff) y salida (JSON temporal)
        3. Llamar cli() con argumentos personalizados
        4. Verificar que:
           - cli() retorna 0 (éxito)
           - Archivo de salida existe
           - JSON es válido y contiene campos esperados
    """
    # Cargar módulo dinámicamente
    mod = load_module()
    # Ruta del archivo de entrada
    sample = ROOT / "data" / "sample.gff"
    # Ruta del archivo de salida (en directorio temporal proporcionado por pytest)
    out = tmp_path / "out.json"
    # call cli with argv
    res = mod.cli([f"--gff={str(sample)}", f"--out={str(out)}"])
    assert res == 0
    assert out.exists()
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["total_features"] == 6
    assert "by_type" in data and "avg_length" in data and "strand_distribution" in data


def test_cli_with_filter(tmp_path):
    """Test: Verificar que CLI respeta argumento --filter-type.
    
    Objetivo:
        Validar que el argumento --filter-type se propaga correctamente desde
        cli() a compute_stats_from_gff() y afecta el resultado JSON.
    
    Parámetros:
        tmp_path: fixture de pytest para directorio temporal
    
    Flujo:
        1. Llamar cli() con argumentos incluyendo --filter-type=CDS
        2. Verificar que:
           - Archivo JSON se crea correctamente
           - total_features = 3 (solo CDS)
           - by_type contiene solo CDS
           - filter_type está presente en JSON
    """
    # Cargar módulo dinámicamente
    mod = load_module()
    # Ruta del archivo de entrada
    sample = ROOT / "data" / "sample.gff"
    # Ruta del archivo de salida (en directorio temporal)
    out = tmp_path / "out_filtered.json"
    # call cli with filter
    res = mod.cli([f"--gff={str(sample)}", f"--out={str(out)}", "--filter-type=CDS"])
    assert res == 0
    assert out.exists()
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["total_features"] == 3
    assert data["by_type"] == {"CDS": 3}
    assert data["filter_type"] == "CDS"
