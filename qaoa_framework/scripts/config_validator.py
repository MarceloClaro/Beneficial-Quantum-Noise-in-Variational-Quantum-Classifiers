"""
Módulo de validação de configurações YAML para QAOA Framework.

Este módulo implementa validação de esquemas JSON para garantir que as
configurações YAML sejam válidas antes da execução.

Task 10: Adaptar Reprodutibilidade - Schema Validation
"""

import yaml
import json
from typing import Dict, Any, List, Tuple
from pathlib import Path


# Schema JSON para experiment_qaoa.yaml
EXPERIMENT_SCHEMA = {
    "type": "object",
    "required": ["run", "reproducibility", "problem", "model", "optimization", "frameworks"],
    "properties": {
        "run": {
            "type": "object",
            "required": ["run_id"],
            "properties": {
                "run_id": {"type": "string"}
            }
        },
        "reproducibility": {
            "type": "object",
            "required": ["seeds"],
            "properties": {
                "seeds": {
                    "type": "array",
                    "items": {"type": "integer"},
                    "minItems": 1
                }
            }
        },
        "problem": {
            "type": "object",
            "required": ["type", "n_nodes", "graph_type"],
            "properties": {
                "type": {"type": "string", "enum": ["maxcut"]},
                "n_nodes": {"type": "integer", "minimum": 2, "maximum": 100},
                "graph_type": {"type": "string", "enum": ["erdos_renyi", "regular", "complete"]},
                "edge_probability": {"type": "number", "minimum": 0.0, "maximum": 1.0},
                "degree": {"type": "integer", "minimum": 1}
            }
        },
        "model": {
            "type": "object",
            "required": ["algorithm", "p_layers"],
            "properties": {
                "algorithm": {"type": "string", "enum": ["qaoa"]},
                "p_layers": {"type": "integer", "minimum": 1, "maximum": 10}
            }
        },
        "noise": {
            "type": "object",
            "properties": {
                "enabled": {"type": "boolean"},
                "model": {"type": "string", "enum": ["depolarizing", "amplitude_damping", "phase_damping", "thermal", "pauli"]},
                "schedule": {"type": "string", "enum": ["constant", "linear", "exponential"]},
                "params": {
                    "type": "object",
                    "properties": {
                        "p": {
                            "type": "array",
                            "items": {"type": "number", "minimum": 0.0, "maximum": 0.1}
                        }
                    }
                }
            }
        },
        "optimization": {
            "type": "object",
            "required": ["optimizer", "maxiter"],
            "properties": {
                "optimizer": {"type": "string", "enum": ["COBYLA", "SLSQP", "Powell", "Nelder-Mead", "L-BFGS-B"]},
                "maxiter": {"type": "integer", "minimum": 1, "maximum": 1000}
            }
        },
        "frameworks": {
            "type": "object",
            "required": ["qiskit"],
            "properties": {
                "qiskit": {
                    "type": "object",
                    "required": ["backend_name"],
                    "properties": {
                        "backend_name": {"type": "string"},
                        "method": {"type": "string", "enum": ["statevector", "matrix_product_state", "stabilizer", "automatic"]},
                        "shots": {"type": "integer", "minimum": 1, "maximum": 100000},
                        "optimization_level": {"type": "integer", "minimum": 0, "maximum": 3},
                        "seed_simulator": {"type": "integer"}
                    }
                }
            }
        }
    }
}


def validar_tipo(valor: Any, tipo_esperado: str, caminho: str) -> Tuple[bool, List[str]]:
    """
    Valida o tipo de um valor.
    
    Args:
        valor: Valor a ser validado
        tipo_esperado: Tipo esperado ('string', 'integer', 'number', 'boolean', 'array', 'object')
        caminho: Caminho no documento para mensagens de erro
        
    Returns:
        Tupla (válido, lista_de_erros)
    """
    erros = []
    
    tipo_mapping = {
        'string': str,
        'integer': int,
        'number': (int, float),
        'boolean': bool,
        'array': list,
        'object': dict
    }
    
    tipo_python = tipo_mapping.get(tipo_esperado)
    if tipo_python is None:
        erros.append(f"{caminho}: Tipo desconhecido '{tipo_esperado}'")
        return False, erros
    
    if not isinstance(valor, tipo_python):
        erros.append(f"{caminho}: Esperado tipo '{tipo_esperado}', recebido '{type(valor).__name__}'")
        return False, erros
    
    return True, []


def validar_enum(valor: Any, opcoes: List, caminho: str) -> Tuple[bool, List[str]]:
    """Valida se o valor está em uma lista de opções."""
    erros = []
    if valor not in opcoes:
        erros.append(f"{caminho}: Valor '{valor}' não é válido. Opções: {opcoes}")
        return False, erros
    return True, []


def validar_range(valor: float, minimo: float = None, maximo: float = None, caminho: str = "") -> Tuple[bool, List[str]]:
    """Valida se o valor está dentro de um range."""
    erros = []
    if minimo is not None and valor < minimo:
        erros.append(f"{caminho}: Valor {valor} menor que o mínimo {minimo}")
        return False, erros
    if maximo is not None and valor > maximo:
        erros.append(f"{caminho}: Valor {valor} maior que o máximo {maximo}")
        return False, erros
    return True, []


def validar_propriedade(valor: Any, schema_prop: Dict, caminho: str) -> Tuple[bool, List[str]]:
    """
    Valida uma propriedade contra seu schema.
    
    Args:
        valor: Valor da propriedade
        schema_prop: Schema da propriedade
        caminho: Caminho no documento
        
    Returns:
        Tupla (válido, lista_de_erros)
    """
    erros = []
    
    # Validar tipo
    if 'type' in schema_prop:
        valido, erros_tipo = validar_tipo(valor, schema_prop['type'], caminho)
        if not valido:
            return False, erros_tipo
    
    # Validar enum
    if 'enum' in schema_prop:
        valido, erros_enum = validar_enum(valor, schema_prop['enum'], caminho)
        if not valido:
            erros.extend(erros_enum)
    
    # Validar range numérico
    if isinstance(valor, (int, float)):
        valido, erros_range = validar_range(
            valor,
            schema_prop.get('minimum'),
            schema_prop.get('maximum'),
            caminho
        )
        if not valido:
            erros.extend(erros_range)
    
    # Validar array
    if schema_prop.get('type') == 'array' and isinstance(valor, list):
        if 'minItems' in schema_prop and len(valor) < schema_prop['minItems']:
            erros.append(f"{caminho}: Array deve ter pelo menos {schema_prop['minItems']} itens")
        
        if 'items' in schema_prop:
            for i, item in enumerate(valor):
                valido_item, erros_item = validar_propriedade(item, schema_prop['items'], f"{caminho}[{i}]")
                erros.extend(erros_item)
    
    # Validar object
    if schema_prop.get('type') == 'object' and isinstance(valor, dict):
        valido_obj, erros_obj = validar_objeto(valor, schema_prop, caminho)
        erros.extend(erros_obj)
    
    return len(erros) == 0, erros


def validar_objeto(obj: Dict, schema: Dict, caminho: str = "root") -> Tuple[bool, List[str]]:
    """
    Valida um objeto contra um schema.
    
    Args:
        obj: Objeto a ser validado
        schema: Schema JSON
        caminho: Caminho no documento
        
    Returns:
        Tupla (válido, lista_de_erros)
    """
    erros = []
    
    # Validar campos obrigatórios
    if 'required' in schema:
        for campo in schema['required']:
            if campo not in obj:
                erros.append(f"{caminho}: Campo obrigatório '{campo}' não encontrado")
    
    # Validar propriedades
    if 'properties' in schema:
        for campo, valor in obj.items():
            if campo in schema['properties']:
                valido, erros_prop = validar_propriedade(
                    valor,
                    schema['properties'][campo],
                    f"{caminho}.{campo}"
                )
                erros.extend(erros_prop)
    
    return len(erros) == 0, erros


def validar_config_qaoa(config: Dict) -> Tuple[bool, List[str]]:
    """
    Valida uma configuração QAOA contra o schema.
    
    Args:
        config: Dicionário de configuração
        
    Returns:
        Tupla (válido, lista_de_erros)
    
    Exemplo:
        >>> config = yaml.safe_load(open('experiment_qaoa.yaml'))
        >>> valido, erros = validar_config_qaoa(config)
        >>> if not valido:
        ...     for erro in erros:
        ...         print(f"❌ {erro}")
    """
    return validar_objeto(config, EXPERIMENT_SCHEMA)


def validar_config_arquivo(caminho_yaml: str) -> Tuple[bool, List[str]]:
    """
    Valida um arquivo YAML de configuração.
    
    Args:
        caminho_yaml: Caminho para o arquivo YAML
        
    Returns:
        Tupla (válido, lista_de_erros)
        
    Exemplo:
        >>> valido, erros = validar_config_arquivo('configs/experiment_qaoa.yaml')
        >>> if valido:
        ...     print("✅ Configuração válida!")
    """
    try:
        with open(caminho_yaml, 'r') as f:
            config = yaml.safe_load(f)
        
        return validar_config_qaoa(config)
    
    except FileNotFoundError:
        return False, [f"Arquivo não encontrado: {caminho_yaml}"]
    except yaml.YAMLError as e:
        return False, [f"Erro ao parsear YAML: {str(e)}"]
    except Exception as e:
        return False, [f"Erro inesperado: {str(e)}"]


def gerar_relatorio_validacao(caminho_yaml: str, output_file: str = None) -> str:
    """
    Gera um relatório de validação detalhado.
    
    Args:
        caminho_yaml: Caminho para o arquivo YAML
        output_file: Arquivo de saída (opcional)
        
    Returns:
        String com o relatório
        
    Exemplo:
        >>> relatorio = gerar_relatorio_validacao('configs/experiment_qaoa.yaml', 'validacao.txt')
        >>> print(relatorio)
    """
    valido, erros = validar_config_arquivo(caminho_yaml)
    
    relatorio = []
    relatorio.append("=" * 80)
    relatorio.append("RELATÓRIO DE VALIDAÇÃO DE CONFIGURAÇÃO QAOA")
    relatorio.append("=" * 80)
    relatorio.append(f"\nArquivo: {caminho_yaml}\n")
    
    if valido:
        relatorio.append("✅ CONFIGURAÇÃO VÁLIDA\n")
        relatorio.append("Todos os campos obrigatórios estão presentes e válidos.")
        relatorio.append("A configuração está pronta para uso.")
    else:
        relatorio.append("❌ CONFIGURAÇÃO INVÁLIDA\n")
        relatorio.append(f"Encontrados {len(erros)} erro(s):\n")
        for i, erro in enumerate(erros, 1):
            relatorio.append(f"  {i}. {erro}")
    
    relatorio.append("\n" + "=" * 80)
    
    relatorio_str = "\n".join(relatorio)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(relatorio_str)
    
    return relatorio_str


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Uso: python config_validator.py <caminho_yaml>")
        sys.exit(1)
    
    caminho = sys.argv[1]
    relatorio = gerar_relatorio_validacao(caminho)
    print(relatorio)
