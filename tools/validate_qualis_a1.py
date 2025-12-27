#!/usr/bin/env python3
"""
Validador de Conformidade QUALIS A1
Verifica se o artigo científico atende aos critérios QUALIS A1 estabelecidos.

Uso:
    python tools/validate_qualis_a1.py --article artigo_cientifico/ --report validation_report.md
"""

import argparse
import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass
class ValidationResult:
    """Resultado de uma validação."""
    criterion: str
    expected: str
    achieved: str
    percentage: float
    status: str  # "✅", "⚠️", "❌"
    details: str = ""


@dataclass
class QualisA1Criteria:
    """Critérios QUALIS A1."""
    # Estruturais
    min_references: int = 35
    max_references: int = 50
    min_doi_coverage: float = 0.80  # 80%
    min_hypotheses: int = 3
    min_objectives: int = 3
    
    # Extensão (palavras)
    abstract_min: int = 250
    abstract_max: int = 300
    introduction_min: int = 3000
    introduction_max: int = 4000
    literature_min: int = 4000
    literature_max: int = 5000
    methodology_min: int = 4000
    methodology_max: int = 5000
    results_min: int = 3000
    results_max: int = 4000
    discussion_min: int = 4000
    discussion_max: int = 5000
    conclusion_min: int = 1000
    conclusion_max: int = 1500
    
    # Qualidade
    min_tables: int = 5
    min_equations: int = 10
    min_congruence: float = 0.95  # 95%


class QualisA1Validator:
    """Validador de conformidade QUALIS A1."""
    
    def __init__(self, article_path: Path):
        self.article_path = Path(article_path)
        self.criteria = QualisA1Criteria()
        self.results: List[ValidationResult] = []
        self.score: float = 0.0
        
    def validate_all(self) -> Tuple[float, List[ValidationResult]]:
        """Executa todas as validações e retorna pontuação e resultados."""
        print("🔍 Iniciando validação QUALIS A1...")
        print("=" * 80)
        
        # Validações estruturais
        self._validate_references()
        self._validate_hypotheses()
        self._validate_objectives()
        
        # Validações de extensão
        self._validate_abstract()
        self._validate_introduction()
        self._validate_literature()
        self._validate_methodology()
        self._validate_results()
        self._validate_discussion()
        self._validate_conclusion()
        
        # Validações de qualidade
        self._validate_tables()
        self._validate_equations()
        self._validate_congruence()
        
        # Calcular pontuação final
        self._calculate_score()
        
        return self.score, self.results
    
    def _validate_references(self):
        """Valida número e qualidade de referências."""
        refs_file = self.article_path / "fase4_secoes" / "agradecimentos_referencias.md"
        
        if not refs_file.exists():
            self.results.append(ValidationResult(
                criterion="Referências",
                expected=f"{self.criteria.min_references}-{self.criteria.max_references}",
                achieved="0 (arquivo não encontrado)",
                percentage=0.0,
                status="❌",
                details="Arquivo de referências não encontrado"
            ))
            return
        
        content = refs_file.read_text(encoding='utf-8')
        
        # Contar referências (linhas começando com número ou padrão de citação)
        ref_patterns = [
            r'^\d+\.\s+[A-Z]',  # Formato: 1. AUTOR
            r'^\[[^\]]+\]',     # Formato: [AUTOR, ano]
        ]
        
        refs_found = 0
        for pattern in ref_patterns:
            refs_found += len(re.findall(pattern, content, re.MULTILINE))
        
        # Remover duplicatas aproximadas
        refs_found = min(refs_found, len(set(re.findall(r'^[^\n]{20,}', content, re.MULTILINE))))
        
        # Contar DOIs
        dois = re.findall(r'(?:doi:|DOI:|https://doi\.org/)[^\s\n]+', content, re.IGNORECASE)
        doi_coverage = len(dois) / refs_found if refs_found > 0 else 0.0
        
        # Avaliar
        in_range = self.criteria.min_references <= refs_found <= self.criteria.max_references
        doi_ok = doi_coverage >= self.criteria.min_doi_coverage
        
        percentage = 100.0 if in_range and doi_ok else (
            50.0 if in_range or doi_ok else 0.0
        )
        
        status = "✅" if in_range and doi_ok else "⚠️" if in_range or doi_ok else "❌"
        
        self.results.append(ValidationResult(
            criterion="Referências",
            expected=f"{self.criteria.min_references}-{self.criteria.max_references} refs, ≥{self.criteria.min_doi_coverage*100}% DOI",
            achieved=f"{refs_found} refs, {doi_coverage*100:.1f}% DOI",
            percentage=percentage,
            status=status,
            details=f"Encontradas {refs_found} referências com {len(dois)} DOIs"
        ))
    
    def _validate_hypotheses(self):
        """Valida número de hipóteses."""
        hyp_file = self.article_path / "fase3_estrutura" / "hipoteses_objetivos.md"
        
        if not hyp_file.exists():
            self.results.append(ValidationResult(
                criterion="Hipóteses",
                expected=f"≥{self.criteria.min_hypotheses}",
                achieved="0 (arquivo não encontrado)",
                percentage=0.0,
                status="❌"
            ))
            return
        
        content = hyp_file.read_text(encoding='utf-8')
        
        # Contar hipóteses (H₀, H₁, H₂, etc.)
        hypotheses = re.findall(r'H[₀₁₂₃₄₅₆₇₈₉\d]+[:：]', content)
        num_hypotheses = len(set(hypotheses))
        
        percentage = min(100.0, (num_hypotheses / self.criteria.min_hypotheses) * 100.0)
        status = "✅" if num_hypotheses >= self.criteria.min_hypotheses else "❌"
        
        self.results.append(ValidationResult(
            criterion="Hipóteses",
            expected=f"≥{self.criteria.min_hypotheses}",
            achieved=str(num_hypotheses),
            percentage=percentage,
            status=status,
            details=f"Hipóteses identificadas: {', '.join(set(hypotheses[:5]))}"
        ))
    
    def _validate_objectives(self):
        """Valida número de objetivos SMART."""
        obj_file = self.article_path / "fase3_estrutura" / "hipoteses_objetivos.md"
        
        if not obj_file.exists():
            self.results.append(ValidationResult(
                criterion="Objetivos SMART",
                expected=f"≥{self.criteria.min_objectives}",
                achieved="0 (arquivo não encontrado)",
                percentage=0.0,
                status="❌"
            ))
            return
        
        content = obj_file.read_text(encoding='utf-8')
        
        # Contar objetivos
        objectives = re.findall(r'(?:Objetivo|OE)[^\n]{0,20}[₀₁₂₃₄₅₆₇₈₉\d]+[:：]', content, re.IGNORECASE)
        num_objectives = len(set(objectives))
        
        # Se não encontrar padrão específico, procurar seção de objetivos
        if num_objectives == 0:
            obj_section = re.search(r'##\s*Objetivos.*?(?=##|$)', content, re.DOTALL | re.IGNORECASE)
            if obj_section:
                num_objectives = len(re.findall(r'^\s*[-*•]\s+\*\*', obj_section.group(), re.MULTILINE))
        
        percentage = min(100.0, (num_objectives / self.criteria.min_objectives) * 100.0)
        status = "✅" if num_objectives >= self.criteria.min_objectives else "❌"
        
        self.results.append(ValidationResult(
            criterion="Objetivos SMART",
            expected=f"≥{self.criteria.min_objectives}",
            achieved=str(num_objectives),
            percentage=percentage,
            status=status,
            details=f"Objetivos específicos identificados: {num_objectives}"
        ))
    
    def _validate_section(self, section_name: str, file_name: str, min_words: int, max_words: int):
        """Valida extensão de uma seção."""
        section_file = self.article_path / "fase4_secoes" / file_name
        
        if not section_file.exists():
            self.results.append(ValidationResult(
                criterion=section_name,
                expected=f"{min_words}-{max_words} palavras",
                achieved="0 (arquivo não encontrado)",
                percentage=0.0,
                status="❌"
            ))
            return
        
        content = section_file.read_text(encoding='utf-8')
        
        # Remover metadados e contar palavras
        content = re.sub(r'^---.*?---', '', content, flags=re.DOTALL)
        content = re.sub(r'```.*?```', '', content, flags=re.DOTALL)
        words = len(re.findall(r'\b\w+\b', content))
        
        # Avaliar
        in_range = min_words <= words <= max_words
        percentage = 100.0 if in_range else (
            80.0 if abs(words - min_words) < 500 or abs(words - max_words) < 500 else
            50.0 if words > 0 else 0.0
        )
        
        status = "✅" if in_range else "⚠️" if words > 0 else "❌"
        
        self.results.append(ValidationResult(
            criterion=section_name,
            expected=f"{min_words}-{max_words} palavras",
            achieved=f"{words} palavras",
            percentage=percentage,
            status=status,
            details=f"{'Dentro do intervalo' if in_range else 'Fora do intervalo esperado'}"
        ))
    
    def _validate_abstract(self):
        """Valida extensão do abstract."""
        self._validate_section(
            "Abstract/Resumo",
            "resumo_abstract.md",
            self.criteria.abstract_min,
            self.criteria.abstract_max
        )
    
    def _validate_introduction(self):
        """Valida extensão da introdução."""
        self._validate_section(
            "Introdução",
            "introducao_completa.md",
            self.criteria.introduction_min,
            self.criteria.introduction_max
        )
    
    def _validate_literature(self):
        """Valida extensão da revisão de literatura."""
        self._validate_section(
            "Revisão de Literatura",
            "revisao_literatura_completa.md",
            self.criteria.literature_min,
            self.criteria.literature_max
        )
    
    def _validate_methodology(self):
        """Valida extensão da metodologia."""
        self._validate_section(
            "Metodologia",
            "metodologia_completa.md",
            self.criteria.methodology_min,
            self.criteria.methodology_max
        )
    
    def _validate_results(self):
        """Valida extensão dos resultados."""
        self._validate_section(
            "Resultados",
            "resultados_completo.md",
            self.criteria.results_min,
            self.criteria.results_max
        )
    
    def _validate_discussion(self):
        """Valida extensão da discussão."""
        self._validate_section(
            "Discussão",
            "discussao_completa.md",
            self.criteria.discussion_min,
            self.criteria.discussion_max
        )
    
    def _validate_conclusion(self):
        """Valida extensão da conclusão."""
        self._validate_section(
            "Conclusão",
            "conclusao_completa.md",
            self.criteria.conclusion_min,
            self.criteria.conclusion_max
        )
    
    def _validate_tables(self):
        """Valida número de tabelas."""
        results_file = self.article_path / "fase4_secoes" / "resultados_completo.md"
        supp_file = self.article_path / "fase5_suplementar" / "tabelas_suplementares.md"
        
        total_tables = 0
        
        if results_file.exists():
            content = results_file.read_text(encoding='utf-8')
            total_tables += len(re.findall(r'(?:Tabela|Table)\s+\d+', content, re.IGNORECASE))
        
        if supp_file.exists():
            content = supp_file.read_text(encoding='utf-8')
            total_tables += len(re.findall(r'(?:Tabela|Table)\s+S\d+', content, re.IGNORECASE))
        
        percentage = min(100.0, (total_tables / self.criteria.min_tables) * 100.0)
        status = "✅" if total_tables >= self.criteria.min_tables else "⚠️" if total_tables > 0 else "❌"
        
        self.results.append(ValidationResult(
            criterion="Tabelas",
            expected=f"≥{self.criteria.min_tables}",
            achieved=str(total_tables),
            percentage=percentage,
            status=status,
            details=f"Tabelas principais e suplementares"
        ))
    
    def _validate_equations(self):
        """Valida número de equações LaTeX."""
        method_file = self.article_path / "fase4_secoes" / "metodologia_completa.md"
        
        if not method_file.exists():
            self.results.append(ValidationResult(
                criterion="Equações LaTeX",
                expected=f"≥{self.criteria.min_equations}",
                achieved="0 (arquivo não encontrado)",
                percentage=0.0,
                status="❌"
            ))
            return
        
        content = method_file.read_text(encoding='utf-8')
        
        # Contar equações (diversos formatos LaTeX)
        equations = len(re.findall(r'(?:\$\$.*?\$\$|\\\[.*?\\\]|\\begin\{equation\})', content, re.DOTALL))
        
        percentage = min(100.0, (equations / self.criteria.min_equations) * 100.0)
        status = "✅" if equations >= self.criteria.min_equations else "⚠️" if equations > 0 else "❌"
        
        self.results.append(ValidationResult(
            criterion="Equações LaTeX",
            expected=f"≥{self.criteria.min_equations}",
            achieved=str(equations),
            percentage=percentage,
            status=status,
            details="Equações matemáticas formalizadas"
        ))
    
    def _validate_congruence(self):
        """Valida conivência código-texto."""
        cong_file = self.article_path / "fase6_consolidacao" / "relatorio_conivencia.md"
        
        if not cong_file.exists():
            self.results.append(ValidationResult(
                criterion="Conivência Código-Texto",
                expected=f"≥{self.criteria.min_congruence*100}%",
                achieved="N/A (arquivo não encontrado)",
                percentage=0.0,
                status="❌"
            ))
            return
        
        content = cong_file.read_text(encoding='utf-8')
        
        # Procurar percentual de conivência
        congruence_match = re.search(r'(\d+(?:\.\d+)?)\s*%.*?(?:conivência|congruence|consistência)', 
                                      content, re.IGNORECASE)
        
        if congruence_match:
            congruence = float(congruence_match.group(1)) / 100.0
        else:
            congruence = 0.0
        
        percentage = (congruence / self.criteria.min_congruence) * 100.0
        status = "✅" if congruence >= self.criteria.min_congruence else "⚠️" if congruence > 0 else "❌"
        
        self.results.append(ValidationResult(
            criterion="Conivência Código-Texto",
            expected=f"≥{self.criteria.min_congruence*100}%",
            achieved=f"{congruence*100:.1f}%",
            percentage=percentage,
            status=status,
            details="Rastreabilidade código-texto verificada"
        ))
    
    def _calculate_score(self):
        """Calcula pontuação final."""
        if not self.results:
            self.score = 0.0
            return
        
        total_percentage = sum(r.percentage for r in self.results)
        self.score = total_percentage / len(self.results)
    
    def generate_report(self, output_file: Path):
        """Gera relatório de validação."""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("# Relatório de Validação QUALIS A1\n\n")
            f.write(f"**Data:** {Path(output_file).stat().st_mtime}\n")
            f.write(f"**Artigo:** {self.article_path}\n\n")
            
            # Pontuação final
            f.write("## 🏆 Pontuação Final\n\n")
            f.write(f"**{self.score:.1f}/100** ")
            
            if self.score >= 90:
                f.write("(🥇 EXCELENTE - Pronto para Nature/Science/Quantum)\n")
            elif self.score >= 80:
                f.write("(🥈 BOM - Pronto para periódicos Qualis A1)\n")
            elif self.score >= 70:
                f.write("(🥉 REGULAR - Necessita ajustes antes da submissão)\n")
            else:
                f.write("(❌ INSUFICIENTE - Revisão substancial necessária)\n")
            
            f.write("\n---\n\n")
            
            # Resultados detalhados
            f.write("## 📊 Resultados Detalhados\n\n")
            f.write("| Critério | Esperado | Alcançado | % | Status |\n")
            f.write("|----------|----------|-----------|---|--------|\n")
            
            for result in self.results:
                f.write(f"| {result.criterion} | {result.expected} | "
                       f"{result.achieved} | {result.percentage:.1f}% | {result.status} |\n")
            
            f.write("\n---\n\n")
            
            # Detalhes
            f.write("## 📝 Detalhes das Validações\n\n")
            for result in self.results:
                f.write(f"### {result.criterion} {result.status}\n\n")
                f.write(f"- **Esperado:** {result.expected}\n")
                f.write(f"- **Alcançado:** {result.achieved}\n")
                f.write(f"- **Conformidade:** {result.percentage:.1f}%\n")
                if result.details:
                    f.write(f"- **Detalhes:** {result.details}\n")
                f.write("\n")
            
            # Recomendações
            f.write("## 💡 Recomendações\n\n")
            
            failed = [r for r in self.results if r.status == "❌"]
            warning = [r for r in self.results if r.status == "⚠️"]
            
            if not failed and not warning:
                f.write("✅ **Parabéns!** O artigo atende a todos os critérios QUALIS A1.\n\n")
                f.write("**Próximos passos:**\n")
                f.write("1. Revisão final por pares\n")
                f.write("2. Formatação LaTeX para periódico-alvo\n")
                f.write("3. Geração de figuras finais\n")
                f.write("4. Submissão ao periódico\n")
            else:
                if failed:
                    f.write("### ❌ Critérios Críticos (Ação Obrigatória)\n\n")
                    for r in failed:
                        f.write(f"- **{r.criterion}:** {r.details or 'Revisar e corrigir'}\n")
                    f.write("\n")
                
                if warning:
                    f.write("### ⚠️ Critérios de Atenção (Ação Recomendada)\n\n")
                    for r in warning:
                        f.write(f"- **{r.criterion}:** {r.details or 'Considerar ajustes'}\n")
                    f.write("\n")


def main():
    """Função principal."""
    parser = argparse.ArgumentParser(
        description="Validador de Conformidade QUALIS A1"
    )
    parser.add_argument(
        "--article",
        type=str,
        required=True,
        help="Caminho para o diretório do artigo"
    )
    parser.add_argument(
        "--report",
        type=str,
        default="validation_report.md",
        help="Arquivo de saída para o relatório"
    )
    
    args = parser.parse_args()
    
    # Executar validação
    validator = QualisA1Validator(Path(args.article))
    score, results = validator.validate_all()
    
    # Gerar relatório
    output_file = Path(args.report)
    validator.generate_report(output_file)
    
    # Imprimir resumo
    print("\n" + "=" * 80)
    print(f"🏆 Pontuação Final: {score:.1f}/100")
    print(f"📄 Relatório salvo em: {output_file}")
    
    passed = len([r for r in results if r.status == "✅"])
    warning = len([r for r in results if r.status == "⚠️"])
    failed = len([r for r in results if r.status == "❌"])
    
    print(f"✅ Aprovados: {passed}/{len(results)}")
    print(f"⚠️  Atenção: {warning}/{len(results)}")
    print(f"❌ Falhou: {failed}/{len(results)}")
    print("=" * 80)


if __name__ == "__main__":
    main()
