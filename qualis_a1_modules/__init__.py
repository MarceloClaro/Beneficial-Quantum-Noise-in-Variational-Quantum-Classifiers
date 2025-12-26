"""
Qualis A1 Modules for Enhanced Scientific Rigor.

This package provides modules for improving the reproducibility,
mathematical rigor, and auditability of the Beneficial Quantum Noise framework.

Modules:
--------
- validation: Mathematical validation of quantum operators
- reproducibility: Tools for ensuring reproducible experiments
- auditing: Code traceability and audit trail generation
- visualization: Circuit diagram generation
- statistical_extensions: Advanced statistical analysis with Bonferroni correction and power analysis

Version: 8.0-QAI
"""

__version__ = "8.0-QAI"

# Import main functions from modules
try:
    from .validation import validar_operadores_kraus
    from .reproducibility import configurar_seeds_reprodutiveis, gerar_manifesto_execucao
    from .auditing import gerar_tabela_codigo_metodo
    from .visualization import gerar_diagrama_circuito, gerar_todos_diagramas_ansatze
    from .statistical_extensions import (
        testes_post_hoc_com_correcao,
        analise_poder_estatistico,
        calcular_tamanho_amostral_necessario
    )
    
    __all__ = [
        'validar_operadores_kraus',
        'configurar_seeds_reprodutiveis',
        'gerar_manifesto_execucao',
        'gerar_tabela_codigo_metodo',
        'gerar_diagrama_circuito',
        'gerar_todos_diagramas_ansatze',
        'testes_post_hoc_com_correcao',
        'analise_poder_estatistico',
        'calcular_tamanho_amostral_necessario',
    ]
except ImportError as e:
    # Allow package to be imported even if some dependencies are missing
    import warnings
    warnings.warn(f"Some qualis_a1_modules could not be imported: {e}")
    __all__ = []
