"""
Testes para validar a estrutura e conteúdo dos notebooks.

Garante que os notebooks foram criados corretamente e contêm
os elementos necessários para reproduzir o framework completo.
"""

import json
import os
from pathlib import Path
import pytest


def get_notebook_path(notebook_name):
    """Retorna o caminho para um notebook."""
    repo_root = Path(__file__).parent.parent
    return repo_root / "notebooks" / notebook_name


def load_notebook(notebook_name):
    """Carrega um notebook como JSON."""
    nb_path = get_notebook_path(notebook_name)
    if not nb_path.exists():
        pytest.skip(f"Notebook {notebook_name} não encontrado")
    
    with open(nb_path, 'r', encoding='utf-8') as f:
        return json.load(f)


class TestNotebookStructure:
    """Testa a estrutura básica dos notebooks."""
    
    @pytest.mark.parametrize("notebook_name", [
        "01_introducao_vqc.ipynb",
        "02_beneficial_noise_demo.ipynb",
        "03_reproducao_experimentos.ipynb"
    ])
    def test_notebook_exists(self, notebook_name):
        """Verifica se o notebook existe."""
        nb_path = get_notebook_path(notebook_name)
        assert nb_path.exists(), f"Notebook {notebook_name} não encontrado"
    
    @pytest.mark.parametrize("notebook_name", [
        "01_introducao_vqc.ipynb",
        "02_beneficial_noise_demo.ipynb",
        "03_reproducao_experimentos.ipynb"
    ])
    def test_notebook_valid_json(self, notebook_name):
        """Verifica se o notebook é um JSON válido."""
        nb = load_notebook(notebook_name)
        assert isinstance(nb, dict), "Notebook deve ser um dicionário"
        assert "cells" in nb, "Notebook deve ter campo 'cells'"
        assert "metadata" in nb, "Notebook deve ter campo 'metadata'"
    
    @pytest.mark.parametrize("notebook_name,min_cells", [
        ("01_introducao_vqc.ipynb", 15),
        ("02_beneficial_noise_demo.ipynb", 15),
        ("03_reproducao_experimentos.ipynb", 20)
    ])
    def test_notebook_has_cells(self, notebook_name, min_cells):
        """Verifica se o notebook tem número mínimo de células."""
        nb = load_notebook(notebook_name)
        assert len(nb["cells"]) >= min_cells, \
            f"Notebook {notebook_name} deve ter pelo menos {min_cells} células"


class TestNotebookContent:
    """Testa o conteúdo dos notebooks."""
    
    @pytest.mark.parametrize("notebook_name", [
        "01_introducao_vqc.ipynb",
        "02_beneficial_noise_demo.ipynb",
        "03_reproducao_experimentos.ipynb"
    ])
    def test_has_colab_badge(self, notebook_name):
        """Verifica se o notebook tem badge 'Open in Colab'."""
        nb = load_notebook(notebook_name)
        first_cell_source = ''.join(nb["cells"][0]["source"]) if nb["cells"] else ""
        assert "colab" in first_cell_source.lower(), \
            f"Notebook {notebook_name} deve ter badge 'Open in Colab'"
    
    @pytest.mark.parametrize("notebook_name", [
        "01_introducao_vqc.ipynb",
        "02_beneficial_noise_demo.ipynb",
        "03_reproducao_experimentos.ipynb"
    ])
    def test_has_markdown_and_code_cells(self, notebook_name):
        """Verifica se o notebook tem células markdown e code."""
        nb = load_notebook(notebook_name)
        cell_types = [cell["cell_type"] for cell in nb["cells"]]
        assert "markdown" in cell_types, \
            f"Notebook {notebook_name} deve ter células markdown"
        assert "code" in cell_types, \
            f"Notebook {notebook_name} deve ter células de código"
    
    @pytest.mark.parametrize("notebook_name", [
        "01_introducao_vqc.ipynb",
        "02_beneficial_noise_demo.ipynb",
        "03_reproducao_experimentos.ipynb"
    ])
    def test_has_dual_audience_content(self, notebook_name):
        """Verifica se o notebook tem conteúdo para iniciantes e especialistas."""
        nb = load_notebook(notebook_name)
        all_text = ' '.join(''.join(cell.get("source", [])) 
                           for cell in nb["cells"])
        
        # Verificar marcadores de conteúdo para iniciantes
        has_beginner = ('iniciante' in all_text.lower() or 
                       '💡' in all_text or
                       'para iniciantes' in all_text.lower())
        
        # Verificar marcadores de conteúdo para especialistas
        has_expert = ('especialista' in all_text.lower() or 
                     '🎓' in all_text or
                     'para especialistas' in all_text.lower())
        
        assert has_beginner, \
            f"Notebook {notebook_name} deve ter conteúdo para iniciantes"
        assert has_expert, \
            f"Notebook {notebook_name} deve ter conteúdo para especialistas"


class TestNotebook03Framework:
    """Testa se o notebook 03 contém as funções do framework."""
    
    def test_has_constantes_fundamentais(self):
        """Verifica se notebook 03 tem classe ConstantesFundamentais."""
        nb = load_notebook("03_reproducao_experimentos.ipynb")
        all_code = ' '.join(''.join(cell.get("source", [])) 
                           for cell in nb["cells"] 
                           if cell["cell_type"] == "code")
        assert "class ConstantesFundamentais" in all_code, \
            "Notebook 03 deve ter classe ConstantesFundamentais"
    
    def test_has_modelo_ruido(self):
        """Verifica se notebook 03 tem classes de modelo de ruído."""
        nb = load_notebook("03_reproducao_experimentos.ipynb")
        all_code = ' '.join(''.join(cell.get("source", [])) 
                           for cell in nb["cells"] 
                           if cell["cell_type"] == "code")
        assert "class ModeloRuido" in all_code, \
            "Notebook 03 deve ter classe ModeloRuido"
    
    def test_has_classificador_vqc(self):
        """Verifica se notebook 03 tem classe ClassificadorVQC."""
        nb = load_notebook("03_reproducao_experimentos.ipynb")
        all_code = ' '.join(''.join(cell.get("source", [])) 
                           for cell in nb["cells"] 
                           if cell["cell_type"] == "code")
        assert "class ClassificadorVQC" in all_code, \
            "Notebook 03 deve ter classe ClassificadorVQC"
    
    def test_has_carregar_datasets(self):
        """Verifica se notebook 03 tem função carregar_datasets."""
        nb = load_notebook("03_reproducao_experimentos.ipynb")
        all_code = ' '.join(''.join(cell.get("source", [])) 
                           for cell in nb["cells"] 
                           if cell["cell_type"] == "code")
        assert "def carregar_datasets" in all_code, \
            "Notebook 03 deve ter função carregar_datasets"
    
    def test_has_executar_grid_search(self):
        """Verifica se notebook 03 tem função executar_grid_search."""
        nb = load_notebook("03_reproducao_experimentos.ipynb")
        all_code = ' '.join(''.join(cell.get("source", [])) 
                           for cell in nb["cells"] 
                           if cell["cell_type"] == "code")
        assert "def executar_grid_search" in all_code, \
            "Notebook 03 deve ter função executar_grid_search"
    
    def test_has_executar_analises_estatisticas(self):
        """Verifica se notebook 03 tem função executar_analises_estatisticas."""
        nb = load_notebook("03_reproducao_experimentos.ipynb")
        all_code = ' '.join(''.join(cell.get("source", [])) 
                           for cell in nb["cells"] 
                           if cell["cell_type"] == "code")
        assert "def executar_analises_estatisticas" in all_code, \
            "Notebook 03 deve ter função executar_analises_estatisticas"
    
    def test_has_gerar_visualizacoes(self):
        """Verifica se notebook 03 tem função gerar_visualizacoes."""
        nb = load_notebook("03_reproducao_experimentos.ipynb")
        all_code = ' '.join(''.join(cell.get("source", [])) 
                           for cell in nb["cells"] 
                           if cell["cell_type"] == "code")
        assert "def gerar_visualizacoes" in all_code, \
            "Notebook 03 deve ter função gerar_visualizacoes"
    
    def test_has_references(self):
        """Verifica se notebook 03 tem referências científicas."""
        nb = load_notebook("03_reproducao_experimentos.ipynb")
        all_text = ' '.join(''.join(cell.get("source", [])) 
                           for cell in nb["cells"])
        
        # Verificar referências chave
        references = [
            "Nielsen & Chuang",
            "Preskill",
            "Cerezo"
        ]
        
        for ref in references:
            assert ref in all_text, \
                f"Notebook 03 deve ter referência a {ref}"


class TestNotebookQUALISA1Compliance:
    """Testa conformidade com padrões QUALIS A1."""
    
    def test_notebook03_has_latex_equations(self):
        """Verifica se notebook 03 tem equações LaTeX."""
        nb = load_notebook("03_reproducao_experimentos.ipynb")
        all_text = ' '.join(''.join(cell.get("source", [])) 
                           for cell in nb["cells"] 
                           if cell["cell_type"] == "markdown")
        
        # Procurar por delimitadores LaTeX
        has_latex = ("$$" in all_text or 
                    "\\[" in all_text or 
                    "$" in all_text)
        
        assert has_latex, \
            "Notebook 03 deve ter equações LaTeX para rigor QUALIS A1"
    
    @pytest.mark.parametrize("notebook_name", [
        "03_reproducao_experimentos.ipynb"
    ])
    def test_has_scientific_methodology(self, notebook_name):
        """Verifica se notebook tem seções de metodologia científica."""
        nb = load_notebook(notebook_name)
        all_text = ' '.join(''.join(cell.get("source", [])) 
                           for cell in nb["cells"]).lower()
        
        # Palavras-chave de metodologia científica
        keywords = ["método", "experiment", "análise", "estatística", 
                   "reprodut", "rigor"]
        
        found_keywords = sum(1 for kw in keywords if kw in all_text)
        
        assert found_keywords >= 3, \
            f"Notebook {notebook_name} deve ter pelo menos 3 palavras-chave " \
            f"de metodologia científica (encontrado: {found_keywords})"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
