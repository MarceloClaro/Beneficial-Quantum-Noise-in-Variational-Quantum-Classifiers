#!/usr/bin/env python3
"""
Audit Checklist - 0-100 Points Scoring System

This script implements the comprehensive audit checklist for evaluating
manuscript quality according to Qualis A1 standards.

Usage:
    python audit_checklist.py --config config.json --manuscript artigo_cientifico/

Author: MegaPrompt v2.0 Framework
Date: 2025-12-26
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from datetime import datetime


class AuditChecklist:
    """Comprehensive audit checklist for Qualis A1 compliance."""
    
    def __init__(self, config_path: str, manuscript_dir: str):
        self.config_path = Path(config_path)
        self.manuscript_dir = Path(manuscript_dir)
        self.config = self.load_config()
        self.scores = {}
        self.max_scores = {}
        self.checks = {}
        
    def load_config(self) -> Dict:
        """Load configuration from JSON file."""
        with open(self.config_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    def check_reproducibility(self) -> Tuple[int, int]:
        """
        Category 1: Reproducibility (30 points)
        
        Returns:
            Tuple of (score, max_score)
        """
        print("\n📋 Category 1: Reproducibility")
        score = 0
        max_score = 30
        checks = []
        
        # Check 1.1: Environment documented (10 pts)
        print("  Checking environment documentation...")
        env_files = [
            'requirements.txt',
            'environment.yml',
            'Pipfile',
            'pyproject.toml'
        ]
        env_found = any((Path('.') / f).exists() for f in env_files)
        
        if env_found:
            score += 10
            checks.append("✓ Environment file found (10/10)")
        else:
            checks.append("✗ No environment file found (0/10)")
        
        # Check 1.2: Seeds fixed and reported (10 pts)
        print("  Checking seed configuration...")
        seed_files = [
            'qualis_a1_modules/reproducibility.py',
            self.manuscript_dir / 'fase1_analise' / 'manifesto_execucao.json'
        ]
        seeds_documented = any(Path(f).exists() for f in seed_files)
        
        if seeds_documented:
            score += 10
            checks.append("✓ Seeds documented (10/10)")
        else:
            checks.append("✗ Seeds not documented (0/10)")
        
        # Check 1.3: Pipeline executable (10 pts)
        print("  Checking execution scripts...")
        exec_scripts = [
            'framework_investigativo_completo.py',
            'executar_framework.sh',
            'README.md'
        ]
        exec_found = any(Path(f).exists() for f in exec_scripts)
        
        if exec_found:
            score += 10
            checks.append("✓ Execution pipeline documented (10/10)")
        else:
            checks.append("✗ No execution pipeline found (0/10)")
        
        self.scores['reproducibility'] = score
        self.max_scores['reproducibility'] = max_score
        self.checks['reproducibility'] = checks
        
        print(f"  Score: {score}/{max_score}")
        return score, max_score
    
    def check_traceability(self) -> Tuple[int, int]:
        """
        Category 2: Traceability (30 points)
        
        Returns:
            Tuple of (score, max_score)
        """
        print("\n📋 Category 2: Traceability")
        score = 0
        max_score = 30
        checks = []
        
        # Check 2.1: Traceability table (15 pts)
        print("  Checking traceability table...")
        trace_file = self.manuscript_dir / 'fase6_consolidacao' / 'rastreabilidade_completa.md'
        
        if trace_file.exists():
            with open(trace_file, 'r', encoding='utf-8') as f:
                content = f.read()
                # Check if table has essential columns
                has_table = '|' in content and 'Seção' in content and 'Evidência' in content
                if has_table:
                    score += 15
                    checks.append("✓ Traceability table complete (15/15)")
                else:
                    score += 5
                    checks.append("⚠ Traceability table incomplete (5/15)")
        else:
            checks.append("✗ Traceability table not found (0/15)")
        
        # Check 2.2: Code→Method mapping (15 pts)
        print("  Checking code-method mapping...")
        mapping_file = self.manuscript_dir / 'fase6_consolidacao' / 'tabela_codigo_metodo.md'
        
        if mapping_file.exists():
            with open(mapping_file, 'r', encoding='utf-8') as f:
                content = f.read()
                has_mapping = '|' in content and 'Componente' in content
                if has_mapping:
                    score += 15
                    checks.append("✓ Code→Method mapping complete (15/15)")
                else:
                    score += 5
                    checks.append("⚠ Code→Method mapping incomplete (5/15)")
        else:
            checks.append("✗ Code→Method mapping not found (0/15)")
        
        self.scores['traceability'] = score
        self.max_scores['traceability'] = max_score
        self.checks['traceability'] = checks
        
        print(f"  Score: {score}/{max_score}")
        return score, max_score
    
    def check_statistical_rigor(self) -> Tuple[int, int]:
        """
        Category 3: Statistical Rigor (20 points)
        
        Returns:
            Tuple of (score, max_score)
        """
        print("\n📋 Category 3: Statistical Rigor")
        score = 0
        max_score = 20
        checks = []
        
        # Check 3.1: Appropriate tests (5 pts)
        print("  Checking statistical tests...")
        methods_file = self.manuscript_dir / 'fase4_secoes' / 'metodologia_completa.md'
        
        if methods_file.exists():
            with open(methods_file, 'r', encoding='utf-8') as f:
                content = f.read().lower()
                has_tests = any(test in content for test in ['anova', 't-test', 'mann-whitney', 'kruskal'])
                if has_tests:
                    score += 5
                    checks.append("✓ Statistical tests documented (5/5)")
                else:
                    checks.append("✗ No statistical tests mentioned (0/5)")
        else:
            checks.append("✗ Methods file not found (0/5)")
        
        # Check 3.2: Multiple comparison correction (5 pts)
        print("  Checking multiple comparison correction...")
        if methods_file.exists():
            with open(methods_file, 'r', encoding='utf-8') as f:
                content = f.read().lower()
                has_correction = any(corr in content for corr in ['bonferroni', 'holm', 'fdr', 'benjamini'])
                if has_correction:
                    score += 5
                    checks.append("✓ Multiple comparison correction applied (5/5)")
                else:
                    checks.append("✗ No correction for multiple comparisons (0/5)")
        
        # Check 3.3: Confidence intervals (5 pts)
        print("  Checking confidence intervals...")
        results_file = self.manuscript_dir / 'fase4_secoes' / 'resultados_completo.md'
        
        if results_file.exists():
            with open(results_file, 'r', encoding='utf-8') as f:
                content = f.read().lower()
                has_ci = any(ci in content for ci in ['confidence interval', '95% ci', 'intervalo de confiança'])
                if has_ci:
                    score += 5
                    checks.append("✓ Confidence intervals reported (5/5)")
                else:
                    checks.append("✗ No confidence intervals found (0/5)")
        else:
            checks.append("✗ Results file not found (0/5)")
        
        # Check 3.4: Effect sizes (5 pts)
        print("  Checking effect sizes...")
        if results_file.exists():
            with open(results_file, 'r', encoding='utf-8') as f:
                content = f.read().lower()
                has_effect = any(eff in content for eff in ["cohen's d", "hedges' g", "effect size", "eta squared"])
                if has_effect:
                    score += 5
                    checks.append("✓ Effect sizes reported (5/5)")
                else:
                    checks.append("✗ No effect sizes found (0/5)")
        
        self.scores['statistical'] = score
        self.max_scores['statistical'] = max_score
        self.checks['statistical'] = checks
        
        print(f"  Score: {score}/{max_score}")
        return score, max_score
    
    def check_transparency(self) -> Tuple[int, int]:
        """
        Category 4: Transparency (20 points)
        
        Returns:
            Tuple of (score, max_score)
        """
        print("\n📋 Category 4: Transparency")
        score = 0
        max_score = 20
        checks = []
        
        # Check 4.1: Code publicly available (10 pts)
        print("  Checking code availability...")
        readme = Path('README.md')
        
        if readme.exists():
            with open(readme, 'r', encoding='utf-8') as f:
                content = f.read().lower()
                has_repo = 'github.com' in content or 'gitlab' in content or 'zenodo' in content
                if has_repo:
                    score += 10
                    checks.append("✓ Code publicly available (10/10)")
                else:
                    checks.append("✗ No public repository mentioned (0/10)")
        else:
            checks.append("✗ README not found (0/10)")
        
        # Check 4.2: Data availability (5 pts)
        print("  Checking data availability...")
        if readme.exists():
            with open(readme, 'r', encoding='utf-8') as f:
                content = f.read().lower()
                has_data = any(word in content for word in ['dataset', 'data availability', 'zenodo', 'figshare'])
                if has_data:
                    score += 5
                    checks.append("✓ Data availability documented (5/5)")
                else:
                    checks.append("⚠ Data availability not clearly stated (0/5)")
        
        # Check 4.3: Limitations discussed (5 pts)
        print("  Checking threats to validity...")
        discussion_file = self.manuscript_dir / 'fase4_secoes' / 'discussao_completa.md'
        
        if discussion_file.exists():
            with open(discussion_file, 'r', encoding='utf-8') as f:
                content = f.read().lower()
                has_limitations = any(lim in content for lim in ['limitation', 'threat', 'validity', 'limitação'])
                if has_limitations:
                    score += 5
                    checks.append("✓ Limitations and threats to validity discussed (5/5)")
                else:
                    checks.append("✗ No limitations section found (0/5)")
        else:
            checks.append("✗ Discussion file not found (0/5)")
        
        self.scores['transparency'] = score
        self.max_scores['transparency'] = max_score
        self.checks['transparency'] = checks
        
        print(f"  Score: {score}/{max_score}")
        return score, max_score
    
    def generate_report(self, output_path: str):
        """Generate audit report in Markdown format."""
        print("\n📝 Generating audit report...")
        
        # Run all checks
        self.check_reproducibility()
        self.check_traceability()
        self.check_statistical_rigor()
        self.check_transparency()
        
        # Calculate total score
        total_score = sum(self.scores.values())
        total_max = sum(self.max_scores.values())
        percentage = (total_score / total_max * 100) if total_max > 0 else 0
        
        # Generate report
        report = []
        report.append("# Checklist de Auditoria - 0-100 Pontos")
        report.append("")
        report.append(f"**Data**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report.append(f"**Configuração**: `{self.config_path}`")
        report.append(f"**Manuscrito**: `{self.manuscript_dir}`")
        report.append("")
        
        # Overall score
        report.append("## 🎯 Pontuação Final")
        report.append("")
        report.append(f"### **{total_score}/{total_max} pontos ({percentage:.1f}%)**")
        report.append("")
        
        if percentage >= 90:
            report.append("✅ **EXCELENTE** - Pronto para submissão a periódicos Qualis A1")
        elif percentage >= 80:
            report.append("⚠️  **BOM** - Pequenos ajustes recomendados antes da submissão")
        elif percentage >= 70:
            report.append("⚠️  **ACEITÁVEL** - Requer melhorias significativas")
        else:
            report.append("❌ **INSUFICIENTE** - Revisão completa necessária")
        
        report.append("")
        
        # Category breakdown
        report.append("## 📊 Detalhamento por Categoria")
        report.append("")
        
        categories = [
            ('reproducibility', 'Reprodutibilidade'),
            ('traceability', 'Rastreabilidade'),
            ('statistical', 'Rigor Estatístico'),
            ('transparency', 'Transparência')
        ]
        
        for key, name in categories:
            score = self.scores.get(key, 0)
            max_score = self.max_scores.get(key, 0)
            pct = (score / max_score * 100) if max_score > 0 else 0
            
            report.append(f"### {name}: {score}/{max_score} pts ({pct:.0f}%)")
            report.append("")
            
            for check in self.checks.get(key, []):
                report.append(f"- {check}")
            report.append("")
        
        # Recommendations
        report.append("## 💡 Recomendações")
        report.append("")
        
        if percentage >= 90:
            report.append("Excelente trabalho! O manuscrito atende aos mais altos padrões:")
            report.append("1. Revisar uma última vez antes da submissão")
            report.append("2. Verificar formatação específica do periódico-alvo")
            report.append("3. Preparar carta de apresentação (cover letter)")
        elif percentage >= 80:
            report.append("Bom progresso! Ajustes recomendados:")
            report.append("1. Completar itens marcados com ✗")
            report.append("2. Revisar itens marcados com ⚠")
            report.append("3. Executar novamente após correções")
        else:
            report.append("Ação necessária:")
            report.append("1. Focar em categorias com pontuação < 70%")
            report.append("2. Revisar FAQ e guias de troubleshooting")
            report.append("3. Considerar template de exemplo fornecido")
        
        report.append("")
        report.append("## 📋 Checklist Detalhado")
        report.append("")
        report.append("### Reprodutibilidade (30 pts)")
        report.append("- [ ] Ambiente documentado (10 pts)")
        report.append("- [ ] Seeds fixas e reportadas (10 pts)")
        report.append("- [ ] Pipeline executável (10 pts)")
        report.append("")
        report.append("### Rastreabilidade (30 pts)")
        report.append("- [ ] Tabela de rastreabilidade completa (15 pts)")
        report.append("- [ ] Mapa código→método completo (15 pts)")
        report.append("")
        report.append("### Rigor Estatístico (20 pts)")
        report.append("- [ ] Testes apropriados (5 pts)")
        report.append("- [ ] Correção para múltiplas comparações (5 pts)")
        report.append("- [ ] Intervalos de confiança (5 pts)")
        report.append("- [ ] Tamanhos de efeito (5 pts)")
        report.append("")
        report.append("### Transparência (20 pts)")
        report.append("- [ ] Código disponível publicamente (10 pts)")
        report.append("- [ ] Dados disponíveis publicamente (5 pts)")
        report.append("- [ ] Limitações e ameaças à validade discutidas (5 pts)")
        report.append("")
        
        report.append("---")
        report.append("*Gerado automaticamente por MegaPrompt v2.0 Audit Checklist*")
        
        # Write report
        os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(report))
        
        print(f"✓ Report generated: {output_path}")
        print(f"\n🎯 Final Score: {total_score}/{total_max} ({percentage:.1f}%)")
        
        return total_score, total_max


def main():
    parser = argparse.ArgumentParser(
        description='Run comprehensive audit checklist for Qualis A1 compliance'
    )
    parser.add_argument(
        '--config',
        type=str,
        default='config.json',
        help='Path to config.json (default: config.json)'
    )
    parser.add_argument(
        '--manuscript',
        type=str,
        default='artigo_cientifico/',
        help='Path to manuscript directory (default: artigo_cientifico/)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='artigo_cientifico/fase6_consolidacao/checklist_auditoria_100pts.md',
        help='Output report path'
    )
    
    args = parser.parse_args()
    
    print("="*60)
    print("AUDIT CHECKLIST - MEGAPROMPT V2.0")
    print("="*60)
    
    # Run audit
    auditor = AuditChecklist(args.config, args.manuscript)
    score, max_score = auditor.generate_report(args.output)
    
    print("\n✅ Audit complete!")
    print(f"   Report: {args.output}")
    print(f"   Score: {score}/{max_score} ({score/max_score*100:.1f}%)")
    
    # Exit code based on score
    if score >= 90:
        sys.exit(0)
    elif score >= 70:
        sys.exit(1)
    else:
        sys.exit(2)


if __name__ == '__main__':
    main()
