#!/usr/bin/env python3
"""
Enhanced Code Analyzer for Scientific Article Generation
Extracts detailed information from Beneficial Quantum Noise repository
"""

import ast
import re
from pathlib import Path
from typing import Dict, List, Tuple, Set
import json

class CodeAnalyzer:
    """Analyzes Python code to extract components for scientific article."""
    
    def __init__(self, repo_path: str):
        self.repo_path = Path(repo_path)
        self.python_files = list(self.repo_path.glob("**/*.py"))
        
    def extract_ansatze(self) -> List[Dict]:
        """Extract quantum ansatz/circuit definitions."""
        ansatze = []
        patterns = [
            r'BasicEntangler',
            r'StronglyEntangling',
            r'RandomEntangling',
            r'SimplifiedTwoDesign',
            r'IQPEmbedding',
            r'AngleEmbedding',
            r'AmplitudeEmbedding'
        ]
        
        for py_file in self.python_files:
            content = py_file.read_text(encoding='utf-8', errors='ignore')
            for pattern in patterns:
                if re.search(pattern, content):
                    # Extract context
                    matches = re.finditer(rf'{pattern}.*', content)
                    for match in matches:
                        ansatze.append({
                            'name': pattern,
                            'file': str(py_file.relative_to(self.repo_path)),
                            'line': content[:match.start()].count('\n') + 1,
                            'context': match.group(0)[:100]
                        })
                        break  # One per file is enough
        
        # Remove duplicates by name
        seen = set()
        unique_ansatze = []
        for a in ansatze:
            if a['name'] not in seen:
                seen.add(a['name'])
                unique_ansatze.append(a)
        
        return unique_ansatze
    
    def extract_noise_models(self) -> List[Dict]:
        """Extract noise model implementations."""
        noise_models = []
        
        # Search for class definitions with exact class names
        class_patterns = [
            ('Depolarizing', r'class\s+RuidoDepolarizante'),
            ('AmplitudeDamping', r'class\s+RuidoAmplitudeDamping'),
            ('PhaseDamping', r'class\s+RuidoPhaseDamping'),
            ('BitFlip', r'class\s+RuidoBitFlip'),
            ('PhaseFlip', r'class\s+RuidoPhaseFlip'),
        ]
        
        for py_file in self.python_files:
            content = py_file.read_text(encoding='utf-8', errors='ignore')
            for name, pattern in class_patterns:
                matches = list(re.finditer(pattern, content))
                if matches:
                    match = matches[0]
                    noise_models.append({
                        'name': name,
                        'file': str(py_file.relative_to(self.repo_path)),
                        'line': content[:match.start()].count('\n') + 1,
                        'context': match.group(0)
                    })
        
        # Remove duplicates
        seen = set()
        unique = []
        for nm in noise_models:
            if nm['name'] not in seen:
                seen.add(nm['name'])
                unique.append(nm)
        
        return unique
    
    def extract_datasets(self) -> List[Dict]:
        """Extract dataset information."""
        datasets = []
        patterns = [
            (r'make_moons', 'Moons', 'sklearn.datasets'),
            (r'make_circles', 'Circles', 'sklearn.datasets'),
            (r'make_blobs', 'Blobs', 'sklearn.datasets'),
            (r'load_iris', 'Iris', 'sklearn.datasets'),
            (r'load_wine', 'Wine', 'sklearn.datasets'),
            (r'load_breast_cancer', 'Breast Cancer', 'sklearn.datasets')
        ]
        
        for py_file in self.python_files:
            content = py_file.read_text(encoding='utf-8', errors='ignore')
            for pattern, name, source in patterns:
                if re.search(pattern, content):
                    # Try to extract parameters
                    param_match = re.search(rf'{pattern}\s*\([^)]*n_samples\s*=\s*(\d+)', content)
                    n_samples = param_match.group(1) if param_match else "[não especificado]"
                    
                    datasets.append({
                        'name': name,
                        'source': source,
                        'n_samples': n_samples,
                        'file': str(py_file.relative_to(self.repo_path)),
                        'line': content[:content.find(pattern)].count('\n') + 1
                    })
                    break
        
        # Remove duplicates
        seen = set()
        unique = []
        for ds in datasets:
            if ds['name'] not in seen:
                seen.add(ds['name'])
                unique.append(ds)
        
        return unique
    
    def extract_metrics(self) -> List[Dict]:
        """Extract evaluation metrics."""
        metrics = []
        patterns = [
            ('Accuracy', r'accuracy_score|acc_test|acc_train'),
            ('F1-Score', r'f1_score'),
            ('Precision', r'precision_score'),
            ('Recall', r'recall_score'),
            ('Loss/Cross-Entropy', r'cross_entropy|loss_fn|calculate_cost')
        ]
        
        for py_file in self.python_files:
            content = py_file.read_text(encoding='utf-8', errors='ignore')
            for name, pattern in patterns:
                if re.search(pattern, content, re.IGNORECASE):
                    metrics.append({
                        'name': name,
                        'file': str(py_file.relative_to(self.repo_path))
                    })
                    break
        
        # Remove duplicates
        seen = set()
        unique = []
        for m in metrics:
            if m['name'] not in seen:
                seen.add(m['name'])
                unique.append(m)
        
        return unique
    
    def extract_schedules(self) -> List[Dict]:
        """Extract noise schedule implementations."""
        schedules = []
        # Search for specific schedule names referenced in code
        patterns = [
            ('Static', r"schedule.*None|sem.*schedule"),
            ('Linear', r"'linear'"),
            ('Exponential', r"'exponencial'"),
            ('Cosine', r"'cosine'"),
        ]
        
        for py_file in self.python_files:
            content = py_file.read_text(encoding='utf-8', errors='ignore')
            for name, pattern in patterns:
                matches = list(re.finditer(pattern, content))
                if matches:
                    match = matches[0]
                    schedules.append({
                        'name': name,
                        'file': str(py_file.relative_to(self.repo_path)),
                        'line': content[:match.start()].count('\n') + 1,
                        'context': match.group(0)[:50]
                    })
        
        # Remove duplicates
        seen = set()
        unique = []
        for s in schedules:
            if s['name'] not in seen:
                seen.add(s['name'])
                unique.append(s)
        
        return unique
    
    def extract_seeds(self) -> List[int]:
        """Extract random seeds used in code."""
        seeds = []
        
        for py_file in self.python_files:
            content = py_file.read_text(encoding='utf-8', errors='ignore')
            # Look for seed assignments
            seed_matches = re.finditer(r'seed\s*=\s*(\d+)|random_state\s*=\s*(\d+)', content, re.IGNORECASE)
            for match in seed_matches:
                seed = int(match.group(1) or match.group(2))
                if seed not in seeds:
                    seeds.append(seed)
        
        return sorted(seeds)
    
    def calculate_configurations(self) -> Dict:
        """Calculate total experimental configurations."""
        ansatze = self.extract_ansatze()
        noise_models = self.extract_noise_models()
        schedules = self.extract_schedules()
        datasets = self.extract_datasets()
        seeds = self.extract_seeds()
        
        # Try to extract noise parameter values
        noise_params = []
        for py_file in self.python_files:
            content = py_file.read_text(encoding='utf-8', errors='ignore')
            # Look for noise parameter ranges
            param_matches = re.finditer(r'noise.*?=.*?\[([\d\.,\s]+)\]', content, re.IGNORECASE)
            for match in param_matches:
                params_str = match.group(1)
                params = [float(p.strip()) for p in params_str.split(',')]
                noise_params.extend(params)
        
        # Default if not found
        if not noise_params:
            noise_params = [0.0, 0.01, 0.05, 0.1]  # Common defaults
        
        # Assume initializations (common values)
        n_inits = 8  # Default assumption
        
        total = (len(datasets) * 
                len(ansatze) * 
                len(noise_models) * 
                len(noise_params) *
                len(schedules) *
                n_inits *
                max(len(seeds), 1))
        
        return {
            'total': total,
            'factors': {
                'datasets': len(datasets),
                'ansatze': len(ansatze),
                'noise_models': len(noise_models),
                'noise_params': len(noise_params),
                'schedules': len(schedules),
                'initializations': n_inits,
                'seeds': max(len(seeds), 1)
            },
            'formula': f"{len(datasets)} datasets × {len(ansatze)} ansätze × {len(noise_models)} noise types × {len(noise_params)} noise params × {len(schedules)} schedules × {n_inits} inits × {max(len(seeds), 1)} seeds"
        }
    
    def generate_full_report(self) -> Dict:
        """Generate complete analysis report."""
        return {
            'ansatze': self.extract_ansatze(),
            'noise_models': self.extract_noise_models(),
            'datasets': self.extract_datasets(),
            'metrics': self.extract_metrics(),
            'schedules': self.extract_schedules(),
            'seeds': self.extract_seeds(),
            'configurations': self.calculate_configurations()
        }


def main():
    """Main entry point for code analysis."""
    import sys
    
    repo_path = sys.argv[1] if len(sys.argv) > 1 else "."
    
    analyzer = CodeAnalyzer(repo_path)
    report = analyzer.generate_full_report()
    
    print("=" * 80)
    print("CODE ANALYSIS REPORT")
    print("=" * 80)
    print()
    
    print(f"📊 Ansätze Found: {len(report['ansatze'])}")
    for ansatz in report['ansatze']:
        print(f"  - {ansatz['name']} ({ansatz['file']}:L{ansatz['line']})")
    print()
    
    print(f"🔊 Noise Models Found: {len(report['noise_models'])}")
    for noise in report['noise_models']:
        print(f"  - {noise['name']} ({noise['file']}:L{noise['line']})")
    print()
    
    print(f"📁 Datasets Found: {len(report['datasets'])}")
    for ds in report['datasets']:
        print(f"  - {ds['name']} (n={ds['n_samples']}) [{ds['file']}]")
    print()
    
    print(f"📏 Metrics Found: {len(report['metrics'])}")
    for metric in report['metrics']:
        print(f"  - {metric['name']}")
    print()
    
    print(f"📅 Schedules Found: {len(report['schedules'])}")
    for sched in report['schedules']:
        print(f"  - {sched['name']}")
    print()
    
    print(f"🎲 Seeds Found: {report['seeds']}")
    print()
    
    print("🔢 Configuration Calculation:")
    config = report['configurations']
    print(f"  Formula: {config['formula']}")
    print(f"  Total: {config['total']} configurations")
    print()
    
    # Save to JSON
    output_file = Path(repo_path) / "code_analysis_report.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    
    print(f"✅ Report saved to: {output_file}")


if __name__ == '__main__':
    main()
