#!/bin/bash
# ============================================================================
# Script para executar o Framework Investigativo Completo v7.2
# ============================================================================

set -e  # Exit on error

# Cores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}============================================================================${NC}"
echo -e "${BLUE} FRAMEWORK INVESTIGATIVO COMPLETO v7.2${NC}"
echo -e "${BLUE} Beneficial Quantum Noise in Variational Quantum Classifiers${NC}"
echo -e "${BLUE}============================================================================${NC}"
echo ""

# Verificar se Python está instalado
if ! command -v python &> /dev/null; then
    echo -e "${RED}❌ Python não encontrado! Por favor, instale Python 3.9+${NC}"
    exit 1
fi

PYTHON_VERSION=$(python --version 2>&1 | awk '{print $2}')
echo -e "${GREEN}✓ Python $PYTHON_VERSION detectado${NC}"

# Verificar se pip está instalado
if ! command -v pip &> /dev/null; then
    echo -e "${RED}❌ pip não encontrado! Por favor, instale pip${NC}"
    exit 1
fi

# Instalar dependências se necessário
echo -e "\n${YELLOW}Verificando dependências...${NC}"
if ! python -c "import pennylane" &> /dev/null; then
    echo -e "${YELLOW}📦 Instalando dependências (isso pode levar alguns minutos)...${NC}"
    pip install -q -r requirements.txt
    echo -e "${GREEN}✓ Dependências instaladas${NC}"
else
    echo -e "${GREEN}✓ Dependências já instaladas${NC}"
fi

# Mostrar menu de opções
echo -e "\n${BLUE}Escolha o modo de execução:${NC}"
echo -e "  ${YELLOW}1)${NC} Modo Rápido Bayesiano (recomendado para teste - ~15 minutos)"
echo -e "  ${YELLOW}2)${NC} Modo Bayesiano Completo (~1-2 horas)"
echo -e "  ${YELLOW}3)${NC} Modo Grid Search Rápido (~5-6 horas)"
echo -e "  ${YELLOW}4)${NC} Modo Grid Search Completo (~15-20 horas)"
echo -e "  ${YELLOW}5)${NC} Modo Híbrido - Grid + Bayesiano (~20-25 horas)"
echo -e "  ${YELLOW}6)${NC} Personalizado (você define os parâmetros)"
echo ""

read -p "Digite sua escolha [1-6]: " choice

case $choice in
    1)
        echo -e "\n${GREEN}Executando Modo Rápido Bayesiano...${NC}"
        export VQC_QUICK=1
        export VQC_BAYESIAN=1
        python framework_investigativo_completo.py --bayes --trials 10 --dataset-bayes moons || exit 1
        ;;
    2)
        echo -e "\n${GREEN}Executando Modo Bayesiano Completo...${NC}"
        export VQC_BAYESIAN=1
        python framework_investigativo_completo.py --bayes --trials 200 --dataset-bayes all || exit 1
        ;;
    3)
        echo -e "\n${GREEN}Executando Modo Grid Search Rápido...${NC}"
        export VQC_QUICK=1
        python framework_investigativo_completo.py || exit 1
        ;;
    4)
        echo -e "\n${GREEN}Executando Modo Grid Search Completo...${NC}"
        python framework_investigativo_completo.py || exit 1
        ;;
    5)
        echo -e "\n${GREEN}Executando Modo Híbrido...${NC}"
        python framework_investigativo_completo.py --run-both || exit 1
        ;;
    6)
        echo -e "\n${YELLOW}Modo Personalizado${NC}"
        read -p "Usar modo rápido? (s/n): " quick
        read -p "Usar otimização Bayesiana? (s/n): " bayes
        
        if [ "$bayes" == "s" ]; then
            read -p "Número de trials (padrão: 100): " trials
            trials=${trials:-100}
            read -p "Dataset (moons/circles/iris/breast_cancer/wine/all): " dataset
            dataset=${dataset:-moons}
            
            if [ "$quick" == "s" ]; then
                export VQC_QUICK=1
            fi
            export VQC_BAYESIAN=1
            python framework_investigativo_completo.py --bayes --trials $trials --dataset-bayes $dataset || exit 1
        else
            if [ "$quick" == "s" ]; then
                export VQC_QUICK=1
            fi
            python framework_investigativo_completo.py || exit 1
        fi
        ;;
    *)
        echo -e "${RED}❌ Opção inválida!${NC}"
        exit 1
        ;;
esac

# Verificar se a execução foi bem-sucedida
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo -e "\n${GREEN}============================================================================${NC}"
    echo -e "${GREEN}✅ FRAMEWORK EXECUTADO COM SUCESSO!${NC}"
    echo -e "${GREEN}============================================================================${NC}"
    
    # Encontrar o diretório de resultados mais recente
    RESULTS_DIR=$(ls -td resultados_* 2>/dev/null | head -1)
    if [ -n "$RESULTS_DIR" ]; then
        echo -e "\n${BLUE}📁 Resultados salvos em: ${NC}$RESULTS_DIR"
        echo -e "\n${BLUE}Arquivos gerados:${NC}"
        echo -e "  • Visualizações interativas (HTML)"
        echo -e "  • Análises estatísticas (CSV)"
        echo -e "  • Otimização Bayesiana (JSON)"
        echo -e "  • Circuitos quânticos (PNG)"
        echo -e "  • Metadados completos (JSON)"
        
        # Contar arquivos
        HTML_COUNT=$(find "$RESULTS_DIR" -name "*.html" 2>/dev/null | wc -l)
        CSV_COUNT=$(find "$RESULTS_DIR" -name "*.csv" 2>/dev/null | wc -l)
        JSON_COUNT=$(find "$RESULTS_DIR" -name "*.json" 2>/dev/null | wc -l)
        
        echo -e "\n${YELLOW}Resumo:${NC}"
        echo -e "  • $HTML_COUNT visualizações HTML"
        echo -e "  • $CSV_COUNT arquivos CSV"
        echo -e "  • $JSON_COUNT arquivos JSON"
    fi
else
    echo -e "\n${RED}============================================================================${NC}"
    echo -e "${RED}❌ ERRO NA EXECUÇÃO DO FRAMEWORK${NC}"
    echo -e "${RED}============================================================================${NC}"
    echo -e "\n${YELLOW}Verifique os logs para mais detalhes.${NC}"
    exit 1
fi
