#!/usr/bin/env python3
"""GFF Statistics Parser

Este módulo proporciona herramientas para analizar archivos GFF (General Feature Format)
y generar estadísticas en formato JSON. Implementa:

Features:
- argparse: interfaz CLI con argumentos opcionales
- Diccionarios y comprensión de listas: estructuras de datos para cálculos
- Manejo básico de archivos: lectura de GFF y escritura de JSON
- Filtrado opcional: permite analizar solo tipos específicos de features

Uso:
    python3 gff_stats.py --gff archivo.gff --out reporte.json [--filter-type TYPE]

Ejemplo:
    # Analizar todos los features
    python3 gff_stats.py --gff genes.gff --out stats.json
    
    # Analizar solo CDS
    python3 gff_stats.py --gff genes.gff --out stats.json --filter-type CDS
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from typing import Dict, Tuple


def compute_stats_from_gff(path: str, filter_type: str | None = None) -> Dict:
    """Analiza un archivo GFF y calcula estadísticas de features.
    
    Esta función realiza un análisis completo del archivo GFF, extrayendo información
    sobre features (genes, CDS, mRNA, etc.) y generando estadísticas agregadas.
    
    Estructura GFF (Tab-separated, 9 columnas):
        1. seqname: nombre del cromosoma/secuencia
        2. source: fuente de la predicción
        3. feature: tipo de feature (gene, CDS, mRNA, etc.) ← USADO PARA FILTRADO
        4. start: posición inicial (1-based) ← USADO PARA CALCULAR LONGITUD
        5. end: posición final ← USADO PARA CALCULAR LONGITUD
        6. score: puntuación
        7. strand: cadena (+, -, .) ← USADO PARA DISTRIBUCIÓN
        8. frame: marco de lectura
        9. attributes: atributos adicionales
    
    Procesamiento:
        - Ignora líneas vacías y comentarios (que comienzan con #)
        - Valida que haya mínimo 9 columnas
        - Calcula longitud = end - start + 1
        - Agrupa datos por tipo de feature
        - Opcionalmente filtra por tipo específico
    
    Args:
        path (str): Ruta al archivo GFF de entrada
        filter_type (str | None): Si se proporciona, solo analiza features de este tipo.
                                  Ejemplo: "gene", "CDS", "mRNA"
                                  Default: None (analiza todos los tipos)
    
    Returns:
        Dict: Diccionario con las siguientes claves:
            - total_features (int): Número total de features no comentados
            - by_type (dict): Conteo de features por tipo
              Ejemplo: {"gene": 210, "CDS": 290, "mRNA": 12}
            - avg_length (dict): Longitud promedio por tipo (redondeada a 1 decimal)
              Ejemplo: {"gene": 890.3, "CDS": 320.7}
            - strand_distribution (dict): Distribución de strands en porcentajes
              Ejemplo: {"+": 61, "-": 39}
            - filter_type (str, opcional): Se añade solo si se aplicó un filtro
    
    Ejemplo:
        >>> stats = compute_stats_from_gff("genes.gff")
        >>> print(stats["total_features"])
        512
        
        >>> stats_cds = compute_stats_from_gff("genes.gff", filter_type="CDS")
        >>> print(stats_cds["by_type"])
        {"CDS": 290}
    """
    # Inicializar contenedores de datos usando defaultdict para evitar KeyError
    total = 0  # Contador total de features válidos
    counts = defaultdict(int)  # Conteo de features por tipo: {"gene": 2, "CDS": 3, ...}
    length_sums = defaultdict(int)  # Suma acumulada de longitudes por tipo
    strand_counts = defaultdict(int)  # Conteo de strands: {"+": N, "-": M, ".": K}

    # Abrir y procesar el archivo GFF línea por línea
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            # Eliminar espacios en blanco al inicio/final (incluye \n)
            line = line.strip()
            
            # Saltar líneas vacías y líneas de comentario (comienzan con #)
            if not line or line.startswith("#"):
                continue
            
            # Dividir por tabulaciones para obtener las 9 columnas GFF
            cols = line.split("\t")
            
            # Validar que haya suficientes columnas
            if len(cols) < 9:
                # Saltar líneas malformadas silenciosamente
                continue
            
            # Extraer campos relevantes (GFF es 1-indexed en la documentación, 0-indexed aquí)
            # Columna 2 (índice 2): feature type
            feature_type = cols[2]
            
            # Aplicar filtro si se especificó: saltar si el tipo no coincide
            if filter_type and feature_type != filter_type:
                continue
            
            # Columnas 3 y 4 (índices 3 y 4): start y end positions
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                # Saltar si start o end no son números válidos
                continue
            
            # Calcular longitud: end - start + 1 (ambos inclusive)
            length = end - start + 1
            
            # Columna 6 (índice 6): strand (+ o -)
            strand = cols[6]

            # Actualizar contadores
            total += 1  # Incrementar contador total
            counts[feature_type] += 1  # Incrementar conteo para este tipo
            length_sums[feature_type] += length  # Acumular longitud para promedios
            
            # Incrementar conteo de strand (distingue +, -, .)
            if strand in ("+", "-"):
                strand_counts[strand] += 1
            else:
                # Para strands inválidos (.), inicializar si no existe
                strand_counts.setdefault(".", 0)

    # CÁLCULOS AGREGADOS USANDO COMPRENSIÓN DE LISTAS

    # 1. Calcular longitud promedio por tipo
    # Comprensión de diccionario: {tipo: promedio_redondeado}
    # Solo incluye tipos con al menos 1 feature (evita división por cero)
    avg_length = {
        k: round(length_sums[k] / counts[k], 1) for k in counts
    } if counts else {}
    # Resultado ejemplo: {"gene": 900.0, "CDS": 84.3, "mRNA": 150.0}

    # 2. Calcular distribución de strands como porcentajes
    # Suma total de strands válidos (+ y -, excluyendo . si existe)
    total_strands = sum(v for k, v in strand_counts.items() if k in ("+", "-"))
    
    if total_strands:
        # Calcular porcentaje para cada strand y redondear a entero
        # Fórmula: (conteo / total) * 100, redondeado al entero más cercano
        strand_distribution = {
            "+": int(round(100 * strand_counts.get("+", 0) / total_strands)),
            "-": int(round(100 * strand_counts.get("-", 0) / total_strands)),
        }
    else:
        # Si no hay strands válidos, inicializar con 0
        strand_distribution = {"+": 0, "-": 0}
    # Resultado ejemplo: {"+": 67, "-": 33}

    # 3. Construir diccionario de resultado
    result = {
        "total_features": total,  # Número total de features procesados
        "by_type": dict(counts),  # Convertir defaultdict a dict normal
        "avg_length": avg_length,  # Diccionario de promedios por tipo
        "strand_distribution": strand_distribution,  # Diccionario de porcentajes
    }
    
    # 4. Añadir información de filtro si se aplicó
    if filter_type:
        result["filter_type"] = filter_type

    return result


def cli(argv: list[str] | None = None) -> int:
    """Interfaz de línea de comandos (CLI) para analizar archivos GFF.
    
    Esta función configura el parser de argumentos, procesa las opciones del usuario,
    llama a compute_stats_from_gff() y escribe el resultado en JSON.
    
    Argumentos CLI (todos opcionales):
        --gff GFF_FILE
            Ruta del archivo GFF de entrada
            Default: "input.gff"
            Ejemplo: --gff genes.gff
        
        --out OUTPUT_FILE
            Ruta del archivo JSON de salida
            Default: "output.json"
            Ejemplo: --out reporte.json
        
        --filter-type FEATURE_TYPE
            Filtro opcional para analizar solo un tipo de feature
            Default: None (analiza todos los tipos)
            Ejemplo: --filter-type CDS
    
    Flujo de ejecución:
        1. Parsear argumentos de línea de comandos
        2. Llamar a compute_stats_from_gff() con los argumentos parseados
        3. Serializar el resultado a formato JSON
        4. Escribir el JSON en el archivo de salida
        5. Retornar código de estado (0 = éxito, no-cero = error)
    
    Args:
        argv (list[str] | None): Lista de argumentos del CLI
                                 Si es None, usa sys.argv (argumentos reales del sistema)
                                 Útil para testing: pasar una lista personalizada
    
    Returns:
        int: Código de estado
             0 = éxito
             No se retornan códigos de error específicos; las excepciones se propagan
    
    Ejemplo de uso (terminal):
        # Analizar todos los features
        $ python3 gff_stats.py --gff datos.gff --out stats.json
        
        # Analizar solo CDS
        $ python3 gff_stats.py --gff datos.gff --out stats.json --filter-type CDS
        
        # Con rutas absolutas
        $ python3 gff_stats.py --gff /ruta/genes.gff --out /ruta/reporte.json
    
    Ejemplo de uso (programático con testing):
        >>> cli(["--gff", "test.gff", "--out", "test_out.json"])
        0
        >>> # El archivo test_out.json ahora contiene las estadísticas en JSON
    """
    # 1. CONFIGURAR ARGUMENTOS CON argparse
    parser = argparse.ArgumentParser(
        description="Compute stats from GFF file and write JSON output"
    )
    
    # Argumento --gff: ruta del archivo GFF de entrada
    parser.add_argument(
        "--gff",
        default="input.gff",  # Valor por defecto si no se especifica
        help="Input GFF file (default: input.gff)"
    )
    
    # Argumento --out: ruta del archivo JSON de salida
    parser.add_argument(
        "--out",
        default="output.json",  # Valor por defecto si no se especifica
        help="Output JSON file (default: output.json)"
    )
    
    # Argumento --filter-type: filtro opcional por tipo de feature
    parser.add_argument(
        "--filter-type",
        default=None,  # No hay filtro por defecto
        help="Filter statistics by feature type (e.g., gene, CDS, mRNA)"
    )
    
    # 2. PARSEAR LOS ARGUMENTOS
    # Si argv es None, argparse usará sys.argv automáticamente
    args = parser.parse_args(argv)

    # 3. LLAMAR A LA FUNCIÓN DE CÁLCULO DE ESTADÍSTICAS
    # Pasar los argumentos parseados a compute_stats_from_gff()
    stats = compute_stats_from_gff(args.gff, filter_type=args.filter_type)
    
    # 4. ESCRIBIR RESULTADO A ARCHIVO JSON
    # Abrir archivo en modo escritura con codificación UTF-8
    with open(args.out, "w", encoding="utf-8") as outfh:
        # json.dump(): serializar el diccionario de estadísticas a JSON
        # indent=2: indentar con 2 espacios para legibilidad
        # ensure_ascii=False: permitir caracteres no-ASCII (ej: ñ, á)
        json.dump(stats, outfh, indent=2, ensure_ascii=False)

    # 5. RETORNAR CÓDIGO DE ESTADO
    return 0  # 0 = éxito


if __name__ == "__main__":
    """Punto de entrada del programa cuando se ejecuta como script.
    
    Esta sección permite ejecutar el módulo directamente desde la línea de comandos:
        python3 gff_stats.py --gff entrada.gff --out salida.json
    
    Explicación:
        - __name__ == "__main__": se ejecuta solo cuando se llama como script directo
          (no cuando se importa como módulo en otro archivo)
        - raise SystemExit(cli()): captura el código de retorno de cli() y lo usa como
          código de salida del proceso (shell exit code)
          - exit code 0: éxito
          - exit code no-cero: error (si cli() retorna otro valor)
    """
    raise SystemExit(cli())
