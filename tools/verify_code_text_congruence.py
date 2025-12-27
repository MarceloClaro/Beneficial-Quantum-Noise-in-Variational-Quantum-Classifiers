#!/usr/bin/env python3
"""
Verificador de Conivência Código-Texto
Verifica correspondência entre código-fonte e artigo científico.

Uso:
    python tools/verify_code_text_congruence.py \
        --code framework_investigativo_completo.py \
        --article artigo_cientifico/ \
        --output congruence_report.md
"""

import argparse
import ast
import re
from pathlib import Path
from typing import Dict, List, Tuple
from dataclasses import dataclass, field


@dataclass
class CodeAnalysis:
    """Análise do código-fonte."""
    total_lines: int = 0
    num_classes: int = 0
    num_functions: int = 0
    class_names: List[str] = field(default_factory=list)
    datasets: List[str] = field(default_factory=list)
    noise_models: List[str] = field(default_factory=list)
    ansatze: List[str] = field(default_factory=list)
    metrics: List[str] = field(default_factory=list)
    libraries: Dict[str, str] = field(default_factory=dict)
    seeds: List[int] = field(default_factory=list)


@dataclass
class ArticleAnalysis:
    """Análise do artigo."""
    num_classes_mentioned: int = 0
    datasets_mentioned: List[str] = field(default_factory=list)
    noise_models_mentioned: List[str] = field(default_factory=list)
    ansatze_mentioned: List[str] = field(default_factory=list)
    metrics_mentioned: List[str] = field(default_factory=list)
    libraries_mentioned: Dict[str, str] = field(default_factory=dict)
    seeds_mentioned: List[int] = field(default_factory=list)


@dataclass
class CongruenceResult:
    """Resultado de verificação de congruência."""
    component: str
    code_value: str
    article_value: str
    match: bool
    percentage: float
    details: str = ""


class CodeTextVerifier:
    """Verificador de conivência código-texto."""
    
    def __init__(self, code_path: Path, article_path: Path):
        self.code_path = Path(code_path)
        self.article_path = Path(article_path)
        self.code_analysis = CodeAnalysis()
        self.article_analysis = ArticleAnalysis()
        self.results: List[CongruenceResult] = []
        
    def analyze_code(self):
        """Analisa o código-fonte."""
        print("🔍 Analisando código-fonte...")
        
        if not self.code_path.exists():
            print(f"❌ Arquivo não encontrado: {self.code_path}")
            return
        
        content = self.code_path.read_text(encoding='utf-8')
        
        # Contar linhas
        self.code_analysis.total_lines = len(content.splitlines())
        
        # Analisar AST (Abstract Syntax Tree)
        try:
            tree = ast.parse(content)
            
            # Contar classes e funções
            for node in ast.walk(tree):
                if isinstance(node, ast.ClassDef):
                    self.code_analysis.num_classes += 1
                    self.code_analysis.class_names.append(node.name)
                elif isinstance(node, ast.FunctionDef):
                    self.code_analysis.num_functions += 1
        except SyntaxError:
            print("⚠️  Erro ao analisar AST, usando regex")
        
        # Usar regex como fallback
        if self.code_analysis.num_classes == 0:
            classes = re.findall(r'^class\s+(\w+)', content, re.MULTILINE)
            self.code_analysis.num_classes = len(classes)
            self.code_analysis.class_names = classes
        
        if self.code_analysis.num_functions == 0:
            functions = re.findall(r'^def\s+(\w+)', content, re.MULTILINE)
            self.code_analysis.num_functions = len(functions)
        
        # Extrair datasets
        dataset_patterns = [
            r'datasets\.make_(\w+)',
            r'load_(\w+)',
            r'["\'](\w+)["\'].*dataset',
        ]
        for pattern in dataset_patterns:
            matches = re.findall(pattern, content, re.IGNORECASE)
            self.code_analysis.datasets.extend([m.lower() for m in matches])
        
        self.code_analysis.datasets = list(set(self.code_analysis.datasets))
        
        # Extrair modelos de ruído
        noise_keywords = ['depolarizing', 'amplitude', 'phase', 'damping', 'bitflip', 'phaseflip']
        for keyword in noise_keywords:
            if keyword.lower() in content.lower():
                self.code_analysis.noise_models.append(keyword)
        
        # Extrair ansätze
        ansatz_keywords = ['BasicEntangling', 'StronglyEntangling', 'SimplifiedTwo', 
                          'RandomLayers', 'Hardware', 'IQP', 'HardwareEfficient']
        for keyword in ansatz_keywords:
            if keyword in content:
                self.code_analysis.ansatze.append(keyword)
        
        # Extrair métricas
        metric_keywords = ['accuracy', 'precision', 'recall', 'f1', 'roc_auc', 
                          'confusion_matrix', 'classification_report']
        for keyword in metric_keywords:
            if keyword in content:
                self.code_analysis.metrics.append(keyword)
        
        # Extrair bibliotecas e versões
        import_pattern = r'^import\s+(\w+)|^from\s+(\w+)'
        imports = re.findall(import_pattern, content, re.MULTILINE)
        for imp in imports:
            lib = imp[0] or imp[1]
            if lib:
                self.code_analysis.libraries[lib] = "unknown"
        
        # Procurar versões em comentários
        version_pattern = r'(\w+)\s*[=:]\s*["\']?(\d+\.\d+\.\d+)'
        versions = re.findall(version_pattern, content)
        for lib, version in versions:
            if lib.lower() in [l.lower() for l in self.code_analysis.libraries]:
                self.code_analysis.libraries[lib] = version
        
        # Extrair seeds
        seed_pattern = r'(?:seed|random_state)\s*=\s*(\d+)'
        seeds = re.findall(seed_pattern, content, re.IGNORECASE)
        self.code_analysis.seeds = [int(s) for s in seeds]
        
        print(f"  ✓ Linhas: {self.code_analysis.total_lines}")
        print(f"  ✓ Classes: {self.code_analysis.num_classes}")
        print(f"  ✓ Funções: {self.code_analysis.num_functions}")
        print(f"  ✓ Datasets: {len(self.code_analysis.datasets)}")
        print(f"  ✓ Modelos de ruído: {len(self.code_analysis.noise_models)}")
    
    def analyze_article(self):
        """Analisa o artigo."""
        print("📄 Analisando artigo...")
        
        if not self.article_path.exists():
            print(f"❌ Diretório não encontrado: {self.article_path}")
            return
        
        # Ler todos os arquivos markdown
        article_content = ""
        for md_file in self.article_path.rglob("*.md"):
            article_content += md_file.read_text(encoding='utf-8') + "\n"
        
        # Contar menções de classes
        for class_name in self.code_analysis.class_names:
            if class_name in article_content:
                self.article_analysis.num_classes_mentioned += 1
        
        # Procurar datasets
        for dataset in ['moons', 'circles', 'iris', 'wine', 'digits', 'breast']:
            if dataset in article_content.lower():
                self.article_analysis.datasets_mentioned.append(dataset)
        
        # Procurar modelos de ruído
        noise_keywords = ['depolarizing', 'amplitude', 'phase', 'damping', 'bitflip', 'phaseflip']
        for keyword in noise_keywords:
            if keyword.lower() in article_content.lower():
                self.article_analysis.noise_models_mentioned.append(keyword)
        
        # Procurar ansätze
        ansatz_keywords = ['BasicEntangling', 'StronglyEntangling', 'SimplifiedTwo', 
                          'RandomLayers', 'Hardware', 'IQP', 'HardwareEfficient', 'ansatz', 'ansätze']
        for keyword in ansatz_keywords:
            if keyword.lower() in article_content.lower():
                self.article_analysis.ansatze_mentioned.append(keyword)
        
        # Procurar métricas
        metric_keywords = ['accuracy', 'precision', 'recall', 'f1', 'roc', 'auc']
        for keyword in metric_keywords:
            if keyword.lower() in article_content.lower():
                self.article_analysis.metrics_mentioned.append(keyword)
        
        # Procurar bibliotecas
        for lib in self.code_analysis.libraries:
            if lib.lower() in article_content.lower():
                # Procurar versão
                version_match = re.search(f'{lib}[^\n]*?(\\d+\\.\\d+\\.\\d+)', article_content, re.IGNORECASE)
                if version_match:
                    self.article_analysis.libraries_mentioned[lib] = version_match.group(1)
                else:
                    self.article_analysis.libraries_mentioned[lib] = "mentioned"
        
        # Procurar seeds
        seed_pattern = r'(?:seed|random)\s*[:\[]\s*(\d+)'
        seeds = re.findall(seed_pattern, article_content, re.IGNORECASE)
        self.article_analysis.seeds_mentioned = [int(s) for s in seeds]
        
        print(f"  ✓ Classes mencionadas: {self.article_analysis.num_classes_mentioned}")
        print(f"  ✓ Datasets mencionados: {len(self.article_analysis.datasets_mentioned)}")
        print(f"  ✓ Modelos de ruído mencionados: {len(self.article_analysis.noise_models_mentioned)}")
    
    def verify_congruence(self):
        """Verifica congruência entre código e artigo."""
        print("🔍 Verificando congruência...")
        
        # Verificar linhas de código
        # (Não é crítico que o artigo mencione o número exato, mas deve ser próximo)
        
        # Verificar classes
        if self.code_analysis.num_classes > 0:
            class_percentage = (self.article_analysis.num_classes_mentioned / 
                               self.code_analysis.num_classes) * 100.0
            self.results.append(CongruenceResult(
                component="Classes Implementadas",
                code_value=f"{self.code_analysis.num_classes} classes",
                article_value=f"{self.article_analysis.num_classes_mentioned} mencionadas",
                match=class_percentage >= 50.0,
                percentage=class_percentage,
                details=f"Classes: {', '.join(self.code_analysis.class_names[:5])}"
            ))
        
        # Verificar datasets
        code_datasets = set(d.lower() for d in self.code_analysis.datasets)
        article_datasets = set(d.lower() for d in self.article_analysis.datasets_mentioned)
        
        if code_datasets:
            dataset_match = code_datasets.intersection(article_datasets)
            dataset_percentage = (len(dataset_match) / len(code_datasets)) * 100.0
            self.results.append(CongruenceResult(
                component="Datasets",
                code_value=f"{len(code_datasets)}: {', '.join(sorted(code_datasets))}",
                article_value=f"{len(article_datasets)}: {', '.join(sorted(article_datasets))}",
                match=dataset_percentage >= 80.0,
                percentage=dataset_percentage,
                details=f"Match: {', '.join(sorted(dataset_match))}"
            ))
        
        # Verificar modelos de ruído
        code_noise = set(n.lower() for n in self.code_analysis.noise_models)
        article_noise = set(n.lower() for n in self.article_analysis.noise_models_mentioned)
        
        if code_noise:
            noise_match = code_noise.intersection(article_noise)
            noise_percentage = (len(noise_match) / len(code_noise)) * 100.0
            self.results.append(CongruenceResult(
                component="Modelos de Ruído",
                code_value=f"{len(code_noise)}: {', '.join(sorted(code_noise))}",
                article_value=f"{len(article_noise)}: {', '.join(sorted(article_noise))}",
                match=noise_percentage >= 80.0,
                percentage=noise_percentage,
                details=f"Match: {', '.join(sorted(noise_match))}"
            ))
        
        # Verificar ansätze
        code_ansatz = set(a.lower() for a in self.code_analysis.ansatze)
        article_ansatz = set(a.lower() for a in self.article_analysis.ansatze_mentioned)
        
        if code_ansatz:
            ansatz_match = code_ansatz.intersection(article_ansatz)
            ansatz_percentage = (len(ansatz_match) / len(code_ansatz)) * 100.0 if code_ansatz else 100.0
            self.results.append(CongruenceResult(
                component="Ansätze Quânticos",
                code_value=f"{len(code_ansatz)} tipos",
                article_value=f"{len(article_ansatz)} mencionados",
                match=ansatz_percentage >= 50.0 or len(article_ansatz) > 0,
                percentage=ansatz_percentage if code_ansatz else 100.0,
                details=f"Implementados: {len(code_ansatz)}, Mencionados: {len(article_ansatz)}"
            ))
        
        # Verificar métricas
        code_metrics = set(m.lower() for m in self.code_analysis.metrics)
        article_metrics = set(m.lower() for m in self.article_analysis.metrics_mentioned)
        
        if code_metrics:
            metrics_match = code_metrics.intersection(article_metrics)
            metrics_percentage = (len(metrics_match) / len(code_metrics)) * 100.0
            self.results.append(CongruenceResult(
                component="Métricas de Avaliação",
                code_value=f"{len(code_metrics)}: {', '.join(sorted(code_metrics))}",
                article_value=f"{len(article_metrics)}: {', '.join(sorted(article_metrics))}",
                match=metrics_percentage >= 70.0,
                percentage=metrics_percentage,
                details=f"Match: {', '.join(sorted(metrics_match))}"
            ))
        
        # Verificar bibliotecas
        code_libs = set(self.code_analysis.libraries.keys())
        article_libs = set(self.article_analysis.libraries_mentioned.keys())
        
        # Focar em bibliotecas principais
        main_libs = {'pennylane', 'qiskit', 'numpy', 'sklearn', 'tensorflow', 'cirq'}
        code_main = code_libs.intersection(main_libs)
        article_main = article_libs.intersection(main_libs)
        
        if code_main:
            lib_match = code_main.intersection(article_main)
            lib_percentage = (len(lib_match) / len(code_main)) * 100.0
            self.results.append(CongruenceResult(
                component="Bibliotecas Principais",
                code_value=f"{len(code_main)}: {', '.join(sorted(code_main))}",
                article_value=f"{len(article_main)}: {', '.join(sorted(article_main))}",
                match=lib_percentage >= 80.0,
                percentage=lib_percentage,
                details=f"Match: {', '.join(sorted(lib_match))}"
            ))
        
        # Verificar seeds de reprodutibilidade
        if self.code_analysis.seeds:
            seeds_match = set(self.code_analysis.seeds).intersection(
                set(self.article_analysis.seeds_mentioned)
            )
            seeds_percentage = (len(seeds_match) / len(set(self.code_analysis.seeds))) * 100.0
            self.results.append(CongruenceResult(
                component="Seeds de Reprodutibilidade",
                code_value=f"{sorted(set(self.code_analysis.seeds))}",
                article_value=f"{sorted(set(self.article_analysis.seeds_mentioned))}",
                match=seeds_percentage >= 80.0,
                percentage=seeds_percentage,
                details=f"Seeds críticos para reprodutibilidade"
            ))
    
    def calculate_overall_congruence(self) -> float:
        """Calcula congruência geral."""
        if not self.results:
            return 0.0
        
        total = sum(r.percentage for r in self.results)
        return total / len(self.results)
    
    def generate_report(self, output_file: Path):
        """Gera relatório de congruência."""
        overall = self.calculate_overall_congruence()
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("# Relatório de Conivência Código-Texto\n\n")
            f.write(f"**Código:** {self.code_path}\n")
            f.write(f"**Artigo:** {self.article_path}\n\n")
            
            # Congruência geral
            f.write("## 🎯 Congruência Geral\n\n")
            f.write(f"**{overall:.1f}%** ")
            
            if overall >= 95.0:
                f.write("(✅ EXCELENTE - Conivência total)\n")
            elif overall >= 85.0:
                f.write("(✅ BOA - Conivência adequada)\n")
            elif overall >= 70.0:
                f.write("(⚠️ REGULAR - Algumas inconsistências)\n")
            else:
                f.write("(❌ INSUFICIENTE - Revisão necessária)\n")
            
            f.write("\n---\n\n")
            
            # Análise do código
            f.write("## 💻 Análise do Código\n\n")
            f.write(f"- **Linhas de Código:** {self.code_analysis.total_lines}\n")
            f.write(f"- **Classes:** {self.code_analysis.num_classes}\n")
            f.write(f"- **Funções:** {self.code_analysis.num_functions}\n")
            f.write(f"- **Datasets:** {len(self.code_analysis.datasets)}\n")
            f.write(f"- **Modelos de Ruído:** {len(self.code_analysis.noise_models)}\n")
            f.write(f"- **Ansätze:** {len(self.code_analysis.ansatze)}\n")
            f.write(f"- **Métricas:** {len(self.code_analysis.metrics)}\n")
            f.write(f"- **Seeds:** {sorted(set(self.code_analysis.seeds))}\n")
            
            f.write("\n---\n\n")
            
            # Resultados detalhados
            f.write("## 📊 Verificação Componente por Componente\n\n")
            f.write("| Componente | Código | Artigo | Congruência | Status |\n")
            f.write("|------------|--------|--------|-------------|--------|\n")
            
            for result in self.results:
                status = "✅" if result.match else "❌"
                f.write(f"| {result.component} | {result.code_value} | "
                       f"{result.article_value} | {result.percentage:.1f}% | {status} |\n")
            
            f.write("\n---\n\n")
            
            # Detalhes
            f.write("## 📝 Detalhes das Verificações\n\n")
            for result in self.results:
                status = "✅" if result.match else "❌"
                f.write(f"### {result.component} {status}\n\n")
                f.write(f"- **Código:** {result.code_value}\n")
                f.write(f"- **Artigo:** {result.article_value}\n")
                f.write(f"- **Congruência:** {result.percentage:.1f}%\n")
                if result.details:
                    f.write(f"- **Detalhes:** {result.details}\n")
                f.write("\n")
            
            # Recomendações
            f.write("## 💡 Recomendações\n\n")
            
            failed = [r for r in self.results if not r.match]
            
            if not failed:
                f.write("✅ **Excelente!** O artigo está perfeitamente alinhado com o código.\n\n")
                f.write("A reprodutibilidade está garantida.\n")
            else:
                f.write("### ⚠️ Inconsistências Detectadas\n\n")
                for r in failed:
                    f.write(f"- **{r.component}:** Revisar e alinhar código e texto\n")
                f.write("\n")
                f.write("**Ação Recomendada:** Atualizar seções do artigo para refletir "
                       "com precisão os componentes implementados no código.\n")


def main():
    """Função principal."""
    parser = argparse.ArgumentParser(
        description="Verificador de Conivência Código-Texto"
    )
    parser.add_argument(
        "--code",
        type=str,
        required=True,
        help="Caminho para o arquivo de código-fonte"
    )
    parser.add_argument(
        "--article",
        type=str,
        required=True,
        help="Caminho para o diretório do artigo"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="congruence_report.md",
        help="Arquivo de saída para o relatório"
    )
    
    args = parser.parse_args()
    
    # Executar verificação
    verifier = CodeTextVerifier(Path(args.code), Path(args.article))
    verifier.analyze_code()
    verifier.analyze_article()
    verifier.verify_congruence()
    
    # Calcular congruência geral
    overall = verifier.calculate_overall_congruence()
    
    # Gerar relatório
    output_file = Path(args.output)
    verifier.generate_report(output_file)
    
    # Imprimir resumo
    print("\n" + "=" * 80)
    print(f"🎯 Congruência Geral: {overall:.1f}%")
    print(f"📄 Relatório salvo em: {output_file}")
    
    passed = len([r for r in verifier.results if r.match])
    failed = len([r for r in verifier.results if not r.match])
    
    print(f"✅ Componentes alinhados: {passed}/{len(verifier.results)}")
    print(f"❌ Inconsistências: {failed}/{len(verifier.results)}")
    print("=" * 80)


if __name__ == "__main__":
    main()
