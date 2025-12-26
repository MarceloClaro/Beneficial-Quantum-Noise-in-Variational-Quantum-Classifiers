#!/bin/bash
# Wrapper script to run the methodological consultant tool

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONSULTOR="$SCRIPT_DIR/consultor_metodologico.py"

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}================================${NC}"
echo -e "${BLUE}Consultor Metodológico Qualis A1${NC}"
echo -e "${BLUE}================================${NC}"
echo ""

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo -e "${YELLOW}⚠️  Python 3 não encontrado. Tentando 'python'...${NC}"
    if ! command -v python &> /dev/null; then
        echo "❌ Erro: Python não encontrado no sistema."
        echo "   Por favor, instale Python 3.9+ e tente novamente."
        exit 1
    fi
    PYTHON_CMD="python"
else
    PYTHON_CMD="python3"
fi

echo -e "${GREEN}✅ Python encontrado: $(which $PYTHON_CMD)${NC}"
echo ""

# Check if consultor script exists
if [ ! -f "$CONSULTOR" ]; then
    echo "❌ Erro: Script consultor_metodologico.py não encontrado."
    echo "   Esperado em: $CONSULTOR"
    exit 1
fi

# Make script executable if needed
chmod +x "$CONSULTOR"

# Run the consultant
echo -e "${BLUE}🚀 Executando consultor metodológico...${NC}"
echo ""

$PYTHON_CMD "$CONSULTOR" "$@"

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}✅ Análise concluída com sucesso!${NC}"
else
    echo -e "${YELLOW}⚠️  Análise concluída com avisos (código: $EXIT_CODE)${NC}"
fi

exit $EXIT_CODE
