#!/usr/bin/env python3
"""
Generate Table S1: Complete Experimental Configuration Matrix

This script generates the comprehensive Table S1 that lists all experimental
configurations tested in the research. It creates a CSV file with one row per
configuration, including all factor combinations.

Usage:
    python generate_s1.py --config config.json --output artigo_cientifico/fase5_suplementar/tabela_s1_configuracoes.csv

Author: MegaPrompt v2.0 Framework
Date: 2025-12-26
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Any
from itertools import product


def load_config(config_path: str) -> Dict[str, Any]:
    """Load configuration from JSON file."""
    with open(config_path, 'r', encoding='utf-8') as f:
        return json.load(f)


def get_experimental_factors() -> Dict[str, List[str]]:
    """
    Define all experimental factors and their levels.
    
    This should match exactly what was implemented in the code.
    Update these lists based on the actual experiments run.
    
    Returns:
        Dictionary mapping factor names to lists of levels
    """
    factors = {
        'dataset': ['moons', 'circles', 'blobs', 'iris'],
        'ansatz': [
            'BasicEntangler',
            'StronglyEntangling',
            'HardwareEfficient',
            'TreeTensorNetwork',
            'QAOAInspired',
            'AlternatingLayers',
            'RandomEntangling'
        ],
        'noise_type': [
            'None',
            'Depolarizing',
            'AmplitudeDamping',
            'PhaseDamping',
            'BitFlip',
            'PhaseFlip',
            'GeneralizedAmplitudeDamping'
        ],
        'noise_level': [0.0, 0.001, 0.005, 0.01, 0.05, 0.1],
        'schedule': ['Constant', 'Linear', 'Cosine'],
        'initialization': [
            'random',
            'zeros',
            'xavier_uniform',
            'he_normal',
            'small_random',
            'pi_4',
            'pi_6',
            'pi_8'
        ],
        'seed': [42, 123, 456, 789, 1024]
    }
    return factors


def calculate_total_configurations(factors: Dict[str, List[str]]) -> int:
    """Calculate total number of configurations."""
    total = 1
    for levels in factors.values():
        total *= len(levels)
    return total


def generate_configuration_matrix(factors: Dict[str, List[str]]) -> List[Dict[str, Any]]:
    """
    Generate all possible combinations of factors.
    
    Args:
        factors: Dictionary of factor names to level lists
        
    Returns:
        List of configuration dictionaries
    """
    factor_names = list(factors.keys())
    factor_levels = [factors[name] for name in factor_names]
    
    configurations = []
    config_id = 1
    
    for combination in product(*factor_levels):
        config = {'config_id': config_id}
        for name, value in zip(factor_names, combination):
            config[name] = value
        
        # Add computed fields
        config['noise_schedule'] = f"{config['schedule']}_{config['noise_level']}"
        config['framework'] = 'PennyLane'  # Update if multi-framework
        config['n_qubits'] = 4  # Update based on actual config
        config['n_layers'] = 2  # Update based on actual config
        config['optimizer'] = 'Adam'  # Update based on actual config
        config['learning_rate'] = 0.01  # Update based on actual config
        config['n_epochs'] = 100  # Update based on actual config
        config['batch_size'] = 32  # Update based on actual config
        
        configurations.append(config)
        config_id += 1
    
    return configurations


def write_csv(configurations: List[Dict[str, Any]], output_path: str):
    """Write configurations to CSV file."""
    if not configurations:
        print("Error: No configurations to write")
        return
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Get all field names
    fieldnames = list(configurations[0].keys())
    
    # Write CSV
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(configurations)
    
    print(f"✓ Table S1 generated: {output_path}")
    print(f"  Total configurations: {len(configurations)}")
    print(f"  Columns: {len(fieldnames)}")


def generate_summary(factors: Dict[str, List[str]], configurations: List[Dict[str, Any]]):
    """Print summary of generated configurations."""
    print("\n" + "="*60)
    print("TABLE S1 GENERATION SUMMARY")
    print("="*60)
    
    print("\nFactor Breakdown:")
    for factor, levels in factors.items():
        print(f"  {factor}: {len(levels)} levels")
        print(f"    → {levels}")
    
    total_calculated = calculate_total_configurations(factors)
    print(f"\nTotal Configurations: {len(configurations)} (calculated: {total_calculated})")
    
    if len(configurations) != total_calculated:
        print("⚠️  WARNING: Mismatch between calculated and generated configurations!")
    else:
        print("✓ Configuration count verified")
    
    print("\nColumns in CSV:")
    if configurations:
        for col in configurations[0].keys():
            print(f"  - {col}")
    
    print("\n" + "="*60)


def main():
    parser = argparse.ArgumentParser(
        description='Generate Table S1: Complete experimental configuration matrix'
    )
    parser.add_argument(
        '--config',
        type=str,
        default='config.json',
        help='Path to config.json file (default: config.json)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='artigo_cientifico/fase5_suplementar/tabela_s1_configuracoes.csv',
        help='Output CSV file path'
    )
    parser.add_argument(
        '--verify-only',
        action='store_true',
        help='Only verify configuration count without generating file'
    )
    
    args = parser.parse_args()
    
    # Load configuration
    try:
        config = load_config(args.config)
        print(f"✓ Loaded configuration from {args.config}")
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)
    
    # Get experimental factors
    factors = get_experimental_factors()
    
    # Calculate total
    total = calculate_total_configurations(factors)
    print(f"\n📊 Experimental Design:")
    print(f"   Total Configurations: {total:,}")
    
    if args.verify_only:
        generate_summary(factors, [])
        return
    
    # Generate configuration matrix
    print("\n⚙️  Generating configuration matrix...")
    configurations = generate_configuration_matrix(factors)
    
    # Write to CSV
    write_csv(configurations, args.output)
    
    # Print summary
    generate_summary(factors, configurations)
    
    print("\n✅ Table S1 generation complete!")
    print(f"   Output: {args.output}")
    print(f"\n💡 Next steps:")
    print("   1. Review the generated CSV file")
    print("   2. Match against actual experimental results")
    print("   3. Add columns for results (accuracy, loss, etc.)")
    print("   4. Include in supplementary material")


if __name__ == '__main__':
    main()
