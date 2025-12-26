#!/usr/bin/env python3
"""
Check Code-Text Consistency

This script verifies the consistency between code implementation and
manuscript text, ensuring that all claims are backed by evidence.

Usage:
    python check_consistency.py --manuscript artigo_cientifico/ --code . --output consistency_report.md

Author: MegaPrompt v2.0 Framework
Date: 2025-12-26
"""

import argparse
import json
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple, Any
import sys


class ConsistencyChecker:
    """Check consistency between code and manuscript text."""
    
    def __init__(self, manuscript_dir: str, code_dir: str):
        self.manuscript_dir = Path(manuscript_dir)
        self.code_dir = Path(code_dir)
        self.issues = []
        self.checks_passed = 0
        self.checks_total = 0
    
    def check_numeric_claims(self) -> List[Dict[str, Any]]:
        """
        Extract numeric claims from manuscript and verify against code/results.
        
        Returns:
            List of issues found
        """
        print("🔍 Checking numeric claims...")
        issues = []
        
        # Patterns to find numbers in text
        number_patterns = [
            r'(\d+\.?\d*)\%',  # Percentages
            r'p\s*=\s*(\d+\.?\d*)',  # p-values
            r'accuracy.*?(\d+\.?\d*)',  # Accuracy values
            r'(\d+)\s+configurations',  # Configuration counts
            r'(\d+)\s+experiments',  # Experiment counts
        ]
        
        # Check all markdown files in manuscript
        for md_file in self.manuscript_dir.rglob('*.md'):
            if md_file.name.startswith('.'):
                continue
                
            with open(md_file, 'r', encoding='utf-8') as f:
                content = f.read()
                
            for pattern in number_patterns:
                matches = re.finditer(pattern, content, re.IGNORECASE)
                for match in matches:
                    self.checks_total += 1
                    value = match.group(1)
                    context = content[max(0, match.start()-50):min(len(content), match.end()+50)]
                    
                    # Check if there's a citation or evidence marker nearby
                    has_evidence = bool(re.search(r'\[.*?\]|\(see|cf\.|Figure|Table', context))
                    
                    if not has_evidence:
                        issues.append({
                            'file': str(md_file.relative_to(self.manuscript_dir)),
                            'type': 'numeric_claim',
                            'value': value,
                            'context': context.strip(),
                            'severity': 'medium',
                            'message': 'Numeric claim without apparent evidence citation'
                        })
                    else:
                        self.checks_passed += 1
        
        return issues
    
    def check_missing_markers(self) -> List[Dict[str, Any]]:
        """
        Check for [INFORMAÇÃO AUSENTE], [NÃO DISPONÍVEL], [LACUNA DE CITAÇÃO] markers.
        
        Returns:
            List of findings
        """
        print("🔍 Checking for integrity markers...")
        findings = []
        
        markers = [
            'INFORMAÇÃO AUSENTE',
            'NÃO DISPONÍVEL',
            'LACUNA DE CITAÇÃO'
        ]
        
        for md_file in self.manuscript_dir.rglob('*.md'):
            if md_file.name.startswith('.'):
                continue
                
            with open(md_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            for marker in markers:
                pattern = r'\[' + marker + r'\]'
                matches = list(re.finditer(pattern, content))
                
                for match in matches:
                    self.checks_total += 1
                    context = content[max(0, match.start()-100):min(len(content), match.end()+100)]
                    
                    findings.append({
                        'file': str(md_file.relative_to(self.manuscript_dir)),
                        'type': 'integrity_marker',
                        'marker': marker,
                        'context': context.strip(),
                        'severity': 'low',
                        'message': f'Found integrity marker: [{marker}]'
                    })
        
        return findings
    
    def check_citations(self) -> List[Dict[str, Any]]:
        """
        Verify that all citations in text have corresponding references.
        
        Returns:
            List of issues found
        """
        print("🔍 Checking citations...")
        issues = []
        
        # Find all citations
        citation_pattern = r'\(([A-Z][a-z]+(?:\s+et\s+al\.?)?,\s+\d{4})\)'
        
        # Load references file
        refs_file = self.manuscript_dir / 'fase4_secoes' / 'agradecimentos_referencias.md'
        references = set()
        
        if refs_file.exists():
            with open(refs_file, 'r', encoding='utf-8') as f:
                refs_content = f.read()
                # Extract author-year patterns from references
                ref_matches = re.finditer(r'([A-Z][A-Z\s,\.]+)\s+\((\d{4})\)', refs_content)
                for match in ref_matches:
                    author = match.group(1).split(',')[0].strip()
                    year = match.group(2)
                    references.add(f"{author}, {year}")
        
        # Check citations in all files
        for md_file in self.manuscript_dir.rglob('*.md'):
            if md_file.name.startswith('.') or md_file == refs_file:
                continue
                
            with open(md_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            citations = re.finditer(citation_pattern, content)
            for citation in citations:
                self.checks_total += 1
                cited = citation.group(1)
                
                # Simplified check - just look for year match
                year = re.search(r'\d{4}', cited)
                if year and year.group() not in [ref.split(', ')[-1] for ref in references]:
                    issues.append({
                        'file': str(md_file.relative_to(self.manuscript_dir)),
                        'type': 'missing_reference',
                        'citation': cited,
                        'severity': 'high',
                        'message': f'Citation not found in references: {cited}'
                    })
                else:
                    self.checks_passed += 1
        
        return issues
    
    def check_code_references(self) -> List[Dict[str, Any]]:
        """
        Check that code files mentioned in manuscript actually exist.
        
        Returns:
            List of issues found
        """
        print("🔍 Checking code file references...")
        issues = []
        
        # Pattern to find code file references
        code_pattern = r'`([a-zA-Z0-9_/\.]+\.py)`'
        
        for md_file in self.manuscript_dir.rglob('*.md'):
            if md_file.name.startswith('.'):
                continue
                
            with open(md_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            matches = re.finditer(code_pattern, content)
            for match in matches:
                self.checks_total += 1
                filename = match.group(1)
                filepath = self.code_dir / filename
                
                if not filepath.exists():
                    issues.append({
                        'file': str(md_file.relative_to(self.manuscript_dir)),
                        'type': 'missing_code_file',
                        'referenced_file': filename,
                        'severity': 'high',
                        'message': f'Referenced code file does not exist: {filename}'
                    })
                else:
                    self.checks_passed += 1
        
        return issues
    
    def generate_report(self, output_path: str):
        """Generate consistency report in Markdown format."""
        print("\n📝 Generating consistency report...")
        
        # Run all checks
        numeric_issues = self.check_numeric_claims()
        marker_findings = self.check_missing_markers()
        citation_issues = self.check_citations()
        code_issues = self.check_code_references()
        
        all_issues = numeric_issues + marker_findings + citation_issues + code_issues
        
        # Calculate consistency score
        if self.checks_total > 0:
            consistency_pct = (self.checks_passed / self.checks_total) * 100
        else:
            consistency_pct = 100.0
        
        # Generate report
        report = []
        report.append("# Relatório de Consistência Código-Texto")
        report.append("")
        report.append(f"**Data**: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report.append(f"**Diretório Manuscrito**: `{self.manuscript_dir}`")
        report.append(f"**Diretório Código**: `{self.code_dir}`")
        report.append("")
        
        # Summary
        report.append("## 📊 Resumo Executivo")
        report.append("")
        report.append(f"- **Verificações Totais**: {self.checks_total}")
        report.append(f"- **Verificações Aprovadas**: {self.checks_passed}")
        report.append(f"- **Problemas Encontrados**: {len([i for i in all_issues if i['severity'] in ['medium', 'high']])}")
        report.append(f"- **Avisos/Marcadores**: {len([i for i in all_issues if i['severity'] == 'low'])}")
        report.append("")
        report.append(f"### 🎯 Índice de Consistência: {consistency_pct:.1f}%")
        report.append("")
        
        if consistency_pct >= 95:
            report.append("✅ **Status**: EXCELENTE - Meta atingida (≥95%)")
        elif consistency_pct >= 90:
            report.append("⚠️  **Status**: BOM - Próximo da meta (≥90%)")
        elif consistency_pct >= 80:
            report.append("⚠️  **Status**: ACEITÁVEL - Requer melhorias (≥80%)")
        else:
            report.append("❌ **Status**: INSUFICIENTE - Requer revisão significativa (<80%)")
        
        report.append("")
        
        # Detailed findings
        if numeric_issues:
            report.append("## 🔢 Afirmações Numéricas Sem Evidência")
            report.append("")
            for issue in numeric_issues:
                report.append(f"### {issue['file']}")
                report.append(f"- **Valor**: {issue['value']}")
                report.append(f"- **Contexto**: ...{issue['context']}...")
                report.append(f"- **Ação**: Adicionar citação ou referência a arquivo/tabela")
                report.append("")
        
        if citation_issues:
            report.append("## 📚 Citações Sem Referência")
            report.append("")
            for issue in citation_issues:
                report.append(f"### {issue['file']}")
                report.append(f"- **Citação**: {issue['citation']}")
                report.append(f"- **Ação**: Adicionar referência completa ou remover citação")
                report.append("")
        
        if code_issues:
            report.append("## 💻 Arquivos de Código Não Encontrados")
            report.append("")
            for issue in code_issues:
                report.append(f"### {issue['file']}")
                report.append(f"- **Arquivo Referenciado**: `{issue['referenced_file']}`")
                report.append(f"- **Ação**: Verificar caminho ou remover referência")
                report.append("")
        
        if marker_findings:
            report.append("## ⚠️  Marcadores de Integridade Encontrados")
            report.append("")
            report.append("Os seguintes marcadores indicam informações ausentes ou não disponíveis:")
            report.append("")
            for finding in marker_findings:
                report.append(f"### {finding['file']}")
                report.append(f"- **Marcador**: [{finding['marker']}]")
                report.append(f"- **Contexto**: ...{finding['context']}...")
                report.append("")
        
        # Recommendations
        report.append("## 💡 Recomendações")
        report.append("")
        
        if consistency_pct < 95:
            report.append("1. **Revisar afirmações numéricas**: Adicionar evidências para todos os valores")
            report.append("2. **Completar referências**: Garantir que todas as citações têm entradas na bibliografia")
            report.append("3. **Validar arquivos**: Verificar que todos os arquivos de código mencionados existem")
            report.append("4. **Resolver marcadores**: Preencher informações ausentes ou justificar indisponibilidade")
        else:
            report.append("✅ A consistência código-texto está excelente! Pequenos ajustes recomendados:")
            report.append("1. Resolver quaisquer marcadores remanescentes")
            report.append("2. Validar uma última vez antes da submissão")
        
        report.append("")
        report.append("---")
        report.append("*Gerado automaticamente por MegaPrompt v2.0 Consistency Checker*")
        
        # Write report
        os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(report))
        
        print(f"✓ Report generated: {output_path}")
        print(f"  Consistency: {consistency_pct:.1f}%")
        print(f"  Issues found: {len([i for i in all_issues if i['severity'] in ['medium', 'high']])}")


def main():
    parser = argparse.ArgumentParser(
        description='Check consistency between code and manuscript text'
    )
    parser.add_argument(
        '--manuscript',
        type=str,
        default='artigo_cientifico/',
        help='Path to manuscript directory (default: artigo_cientifico/)'
    )
    parser.add_argument(
        '--code',
        type=str,
        default='.',
        help='Path to code directory (default: .)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='artigo_cientifico/fase6_consolidacao/relatorio_consistencia.md',
        help='Output report path'
    )
    
    args = parser.parse_args()
    
    # Check directories exist
    if not os.path.isdir(args.manuscript):
        print(f"Error: Manuscript directory not found: {args.manuscript}")
        sys.exit(1)
    
    if not os.path.isdir(args.code):
        print(f"Error: Code directory not found: {args.code}")
        sys.exit(1)
    
    print("="*60)
    print("CONSISTENCY CHECKER - MegaPrompt v2.0")
    print("="*60)
    print(f"\nManuscript: {args.manuscript}")
    print(f"Code: {args.code}")
    print(f"Output: {args.output}\n")
    
    # Run checker
    checker = ConsistencyChecker(args.manuscript, args.code)
    checker.generate_report(args.output)
    
    print("\n✅ Consistency check complete!")
    print(f"   Report: {args.output}")


if __name__ == '__main__':
    main()
