#!/usr/bin/env python3
"""
h5ad_optimize.py

Analyze and optimize h5ad file size by removing unnecessary data
and applying better compression.

Usage:
    python h5ad_optimize.py input.h5ad output.h5ad [--analyze-only]
"""

import argparse
import os
from typing import Any, Dict, List, Optional

import anndata as ad
import numpy as np
import scipy.sparse as sp


def format_size(num_bytes: int) -> str:
    """Return a human-readable size string."""
    size: float = float(num_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024.0:
            return f"{size:3.1f}{unit}"
        size /= 1024.0
    return f"{size:.1f}PB"


def analyze_adata(adata: ad.AnnData, filepath: Optional[str] = None) -> Dict[str, Any]:
    """Analyze an AnnData object to understand size contributors."""

    analysis = {
        "shape": adata.shape,
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
    }

    if filepath and os.path.exists(filepath):
        analysis["file_size"] = os.path.getsize(filepath)

    # Check X matrix
    analysis["X_type"] = type(adata.X).__name__
    analysis["X_dtype"] = str(adata.X.dtype) if hasattr(adata.X, "dtype") else "unknown"

    if sp.issparse(adata.X):
        analysis["X_sparse"] = True
        analysis["X_sparsity"] = 1.0 - (adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1]))
        analysis["X_size_bytes"] = adata.X.data.nbytes + adata.X.indices.nbytes + adata.X.indptr.nbytes
    else:
        analysis["X_sparse"] = False
        if hasattr(adata.X, "shape") and adata.X.size > 0:
            zero_count = np.sum(adata.X == 0)
            analysis["X_sparsity"] = zero_count / adata.X.size
        if hasattr(adata.X, "nbytes"):
            analysis["X_size_bytes"] = adata.X.nbytes

    # Check for raw data
    analysis["has_raw"] = adata.raw is not None
    if adata.raw is not None:
        analysis["raw_shape"] = (adata.raw.n_obs, adata.raw.n_vars)

    # Check layers
    analysis["layers"] = {}
    if adata.layers:
        for key, layer in adata.layers.items():
            layer_info = {
                "shape": layer.shape,
                "dtype": str(layer.dtype),
                "sparse": sp.issparse(layer),
            }
            if sp.issparse(layer):
                layer_info["sparsity"] = 1.0 - (layer.nnz / (layer.shape[0] * layer.shape[1]))
            elif hasattr(layer, "size") and layer.size > 0:
                layer_info["sparsity"] = np.sum(layer == 0) / layer.size
            analysis["layers"][key] = layer_info

    # Check obsm/varm with size estimates
    analysis["obsm_keys"] = {}
    analysis["obsm_total_size"] = 0
    if adata.obsm:
        for key, val in adata.obsm.items():
            info = {
                "shape": val.shape if hasattr(val, "shape") else "unknown",
                "dtype": str(val.dtype) if hasattr(val, "dtype") else "unknown",
            }
            # Estimate memory size
            if hasattr(val, "nbytes"):
                info["size_bytes"] = val.nbytes
                analysis["obsm_total_size"] += val.nbytes
            elif sp.issparse(val):
                info["size_bytes"] = val.data.nbytes + val.indices.nbytes + val.indptr.nbytes
                analysis["obsm_total_size"] += info["size_bytes"]
            analysis["obsm_keys"][key] = info

    analysis["varm_keys"] = {}
    analysis["varm_total_size"] = 0
    if adata.varm:
        for key, val in adata.varm.items():
            info = {
                "shape": val.shape if hasattr(val, "shape") else "unknown",
                "dtype": str(val.dtype) if hasattr(val, "dtype") else "unknown",
            }
            # Estimate memory size
            if hasattr(val, "nbytes"):
                info["size_bytes"] = val.nbytes
                analysis["varm_total_size"] += val.nbytes
            elif sp.issparse(val):
                info["size_bytes"] = val.data.nbytes + val.indices.nbytes + val.indptr.nbytes
                analysis["varm_total_size"] += info["size_bytes"]
            analysis["varm_keys"][key] = info

    # Check metadata
    analysis["n_obs_columns"] = len(adata.obs.columns)
    analysis["n_var_columns"] = len(adata.var.columns)

    return analysis


def print_analysis(analysis: Dict[str, Any]) -> None:
    """Pretty print the analysis results."""
    print("=" * 80)
    print("H5AD FILE ANALYSIS")
    print("=" * 80)

    if "file_size" in analysis:
        print(f"File size: {format_size(analysis['file_size'])}")

    print(f"\nData shape: {analysis['n_obs']:,} cells × {analysis['n_vars']:,} genes")

    print("\n[X Matrix]")
    print(f"  Type: {analysis['X_type']} ({'sparse' if analysis['X_sparse'] else 'dense'})")
    print(f"  Dtype: {analysis['X_dtype']}")
    print(f"  Sparsity: {analysis['X_sparsity']*100:.1f}% zeros")
    if "X_size_bytes" in analysis:
        print(f"  Estimated size: {format_size(analysis['X_size_bytes'])}")

    if analysis["has_raw"]:
        print("\n[Raw Data] ⚠️")
        print(f"  Shape: {analysis['raw_shape'][0]:,} × {analysis['raw_shape'][1]:,}")
        print("  Consider removing if not needed!")

    if analysis["layers"]:
        print(f"\n[Layers] ({len(analysis['layers'])} total)")
        for key, info in analysis["layers"].items():
            sparse_str = "sparse" if info["sparse"] else "dense"
            print(f"  '{key}': {info['shape']} {info['dtype']} ({sparse_str})")
            if "sparsity" in info:
                print(f"    Sparsity: {info['sparsity']*100:.1f}% zeros")

    if analysis["obsm_keys"]:
        print(f"\n[Embeddings/obsm] ({len(analysis['obsm_keys'])} total)")
        print(f"  Total estimated size: {format_size(analysis['obsm_total_size'])}")
        # Sort by size descending
        sorted_obsm = sorted(analysis["obsm_keys"].items(), key=lambda x: x[1].get("size_bytes", 0), reverse=True)
        for key, info in sorted_obsm:
            size_str = f" - {format_size(info['size_bytes'])}" if "size_bytes" in info else ""
            print(f"  '{key}': {info['shape']}{size_str}")

    if analysis["varm_keys"]:
        print(f"\n[varm] ({len(analysis['varm_keys'])} total)")
        print(f"  Total estimated size: {format_size(analysis['varm_total_size'])}")
        # Sort by size descending
        sorted_varm = sorted(analysis["varm_keys"].items(), key=lambda x: x[1].get("size_bytes", 0), reverse=True)
        for key, info in sorted_varm:
            size_str = f" - {format_size(info['size_bytes'])}" if "size_bytes" in info else ""
            print(f"  '{key}': {info['shape']}{size_str}")

    print("\n[Metadata]")
    print(f"  obs columns: {analysis['n_obs_columns']}")
    print(f"  var columns: {analysis['n_var_columns']}")


def suggest_optimizations(analysis: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Suggest optimizations based on analysis."""
    suggestions = []

    # Raw data
    if analysis["has_raw"]:
        suggestions.append(
            {"action": "remove_raw", "description": "Remove raw data (can save significant space)", "impact": "high"}
        )

    # Large embeddings (obsm)
    total_obsm_size = analysis.get("obsm_total_size", 0)
    x_size = analysis.get("X_size_bytes", 1)
    if total_obsm_size > x_size * 0.1:  # If embeddings are >10% of X matrix size
        suggestions.append(
            {
                "action": "review_embeddings",
                "description": f"Review embeddings in obsm ({format_size(total_obsm_size)}) - consider removing unused ones",
                "impact": "high" if total_obsm_size > x_size * 0.5 else "medium",
            }
        )

    # Identify largest embeddings
    for key, info in analysis["obsm_keys"].items():
        size = info.get("size_bytes", 0)
        if size > 1_000_000:  # > 1 MB
            suggestions.append(
                {
                    "action": f"remove_obsm_{key}",
                    "description": f"Consider removing embedding '{key}' ({format_size(size)})",
                    "impact": "medium",
                }
            )

    # Dense to sparse conversion
    if not analysis["X_sparse"] and analysis["X_sparsity"] > 0.5:
        suggestions.append(
            {
                "action": "convert_X_sparse",
                "description": f'Convert X to sparse (currently {analysis["X_sparsity"]*100:.1f}% zeros)',
                "impact": "high",
            }
        )

    # Layer conversions
    for layer_name, info in analysis["layers"].items():
        if not info["sparse"] and info.get("sparsity", 0) > 0.5:
            suggestions.append(
                {
                    "action": f"convert_layer_sparse_{layer_name}",
                    "description": f"Convert layer '{layer_name}' to sparse ({info.get('sparsity', 0)*100:.1f}% zeros)",
                    "impact": "medium",
                }
            )

    # Always suggest recompression
    suggestions.append(
        {"action": "recompress", "description": "Rewrite with optimal compression (gzip level 7-9)", "impact": "medium"}
    )

    return suggestions


def optimize_adata(
    adata: ad.AnnData,
    remove_raw: bool = False,
    convert_sparse: bool = True,
    remove_layers: Optional[List[str]] = None,
    remove_obsm: Optional[List[str]] = None,
    remove_varm: Optional[List[str]] = None,
    sparsity_threshold: float = 0.5,
) -> ad.AnnData:
    """Apply optimizations to AnnData object."""

    print("\nApplying optimizations...")

    # Remove raw data
    if remove_raw and adata.raw is not None:
        print("  ✓ Removing raw data")
        adata.raw = None

    # Remove specified embeddings from obsm
    if remove_obsm and adata.obsm:
        for key in remove_obsm:
            if key in adata.obsm:
                size = adata.obsm[key].nbytes if hasattr(adata.obsm[key], "nbytes") else 0
                print(f"  ✓ Removing obsm['{key}'] ({format_size(size)})")
                del adata.obsm[key]

    # Remove specified embeddings from varm
    if remove_varm and adata.varm:
        for key in remove_varm:
            if key in adata.varm:
                size = adata.varm[key].nbytes if hasattr(adata.varm[key], "nbytes") else 0
                print(f"  ✓ Removing varm['{key}'] ({format_size(size)})")
                del adata.varm[key]

    # Convert X to sparse if needed
    if convert_sparse and not sp.issparse(adata.X) and hasattr(adata.X, "size") and adata.X.size > 0:
        sparsity = np.sum(adata.X == 0) / adata.X.size
        if sparsity > sparsity_threshold:
            print(f"  ✓ Converting X to sparse CSR matrix ({sparsity*100:.1f}% zeros)")
            adata.X = sp.csr_matrix(adata.X)

    # Convert layers to sparse
    if convert_sparse and adata.layers:
        for layer_name, layer in adata.layers.items():
            if not sp.issparse(layer) and hasattr(layer, "size") and layer.size > 0:
                sparsity = np.sum(layer == 0) / layer.size
                if sparsity > sparsity_threshold:
                    print(f"  ✓ Converting layer '{layer_name}' to sparse ({sparsity*100:.1f}% zeros)")
                    adata.layers[layer_name] = sp.csr_matrix(layer)

    # Remove specified layers
    if remove_layers:
        for layer_name in remove_layers:
            if layer_name in adata.layers:
                print(f"  ✓ Removing layer '{layer_name}'")
                del adata.layers[layer_name]

    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze and optimize h5ad file size.")
    parser.add_argument("input", help="Input h5ad file")
    parser.add_argument("output", nargs="?", help="Output h5ad file (required unless --analyze-only)")
    parser.add_argument("--analyze-only", action="store_true", help="Only analyze the file, don't optimize")
    parser.add_argument("--remove-raw", action="store_true", help="Remove raw data")
    parser.add_argument("--no-sparse-conversion", action="store_true", help="Don't convert dense matrices to sparse")
    parser.add_argument(
        "--compression", type=int, default=7, choices=range(0, 10), help="Gzip compression level (0-9, default: 7)"
    )
    parser.add_argument("--remove-layers", nargs="+", help="List of layer names to remove")
    parser.add_argument("--remove-obsm", nargs="+", help="List of obsm (embedding) keys to remove")
    parser.add_argument("--remove-varm", nargs="+", help="List of varm keys to remove")

    args = parser.parse_args()

    if not args.analyze_only and not args.output:
        parser.error("output file is required unless --analyze-only is specified")

    if not os.path.isfile(args.input):
        raise SystemExit(f"File not found: {args.input}")

    print(f"Reading {args.input}...")
    adata = ad.read_h5ad(args.input)

    # Analyze
    analysis = analyze_adata(adata, args.input)
    print_analysis(analysis)

    # Suggest optimizations
    suggestions = suggest_optimizations(analysis)
    if suggestions:
        print("\n" + "=" * 80)
        print("OPTIMIZATION SUGGESTIONS")
        print("=" * 80)
        for i, sug in enumerate(suggestions, 1):
            impact_symbol = "🔴" if sug["impact"] == "high" else "🟡"
            print(f"{i}. {impact_symbol} {sug['description']}")

    if args.analyze_only:
        print("\nAnalysis complete. Use without --analyze-only to create optimized file.")
        return

    # Optimize
    print("\n" + "=" * 80)
    adata = optimize_adata(
        adata,
        remove_raw=args.remove_raw,
        convert_sparse=not args.no_sparse_conversion,
        remove_layers=args.remove_layers,
        remove_obsm=args.remove_obsm,
        remove_varm=args.remove_varm,
    )

    # Write optimized file
    print(f"\nWriting optimized file to {args.output}...")
    print(f"  Compression: gzip level {args.compression}")
    adata.write_h5ad(args.output, compression="gzip", compression_opts=args.compression)

    # Compare sizes
    input_size = os.path.getsize(args.input)
    output_size = os.path.getsize(args.output)
    reduction = (input_size - output_size) / input_size * 100

    print("\n" + "=" * 80)
    print("RESULTS")
    print("=" * 80)
    print(f"Original size: {format_size(input_size)}")
    print(f"Optimized size: {format_size(output_size)}")
    print(f"Reduction: {format_size(input_size - output_size)} ({reduction:.1f}%)")
    print("\n✓ Optimization complete!")


if __name__ == "__main__":
    main()
