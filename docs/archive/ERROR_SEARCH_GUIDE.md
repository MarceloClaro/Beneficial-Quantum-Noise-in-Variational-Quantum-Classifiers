# Error Search Framework - Quick Start Guide

This document explains how to use the Error Search Framework to find and fix errors in the codebase.

## Overview

The Error Search Framework (`error_search_framework.py`) is a comprehensive tool that performs automated error detection across the entire codebase, checking for:

1. **Test failures** - Runs all unit tests and reports failures
2. **Syntax errors** - Validates Python syntax in all `.py` files
3. **Import errors** - Checks if all required dependencies are installed
4. **Code style violations** - Uses Ruff to identify style issues


## Quick Start

### Basic Usage

Run the error search framework without any arguments to perform all checks:

```bash
python error_search_framework.py

```text

### With Auto-Fix

Automatically fix auto-fixable issues (like whitespace):

```bash
python error_search_framework.py --fix

```text

### Detailed Report

Generate a detailed report with full context:

```bash
python error_search_framework.py --detailed

```text

## Output Files

The framework generates two main output files:

1. **ERROR_SEARCH_REPORT.md** - Human-readable markdown report with detailed analysis
2. **error_search_results.json** - Machine-readable JSON with all results


## Understanding the Results

### Overall Status

- **EXCELLENT** - No issues found
- **GOOD** - Minor warnings only
- **NEEDS ATTENTION** - Multiple warnings (>50)
- **CRITICAL** - Critical issues found (syntax errors, missing dependencies, test failures)


### Quality Score

The framework calculates a quality score out of 100 based on:

- Tests: -30 points if failing
- Syntax: -30 points if errors found
- Imports: -20 points if missing dependencies
- Code style: -20 points if violations found


### Example Output

```

================================================================================
ERROR SEARCH FRAMEWORK - FINAL REPORT
================================================================================

Overall Status: GOOD
Quality Score: 80/100

Checks Passed: 3/4
Critical Issues: 0
Warnings: 69

Recommendations:

  1. Consider running ruff with --fix to auto-fix issues


================================================================================
Detailed report saved to: ERROR_SEARCH_REPORT.md
JSON results saved to: error_search_results.json
================================================================================

```text

## Common Workflows

### Before Committing Code

Run the framework to ensure your changes don't introduce errors:

```bash

# Run with auto-fix to clean up minor issues
python error_search_framework.py --fix

# Check if all tests pass
python -m pytest tests/ -v

```text

### Before Creating a Pull Request

Generate a detailed report to include in your PR:

```bash

# Generate detailed report
python error_search_framework.py --detailed

# Review the report
cat ERROR_SEARCH_REPORT.md

```text

### Continuous Integration

Add the framework to your CI pipeline:

```yaml

# .github/workflows/error-check.yml
- name: Run Error Search Framework

  run: python error_search_framework.py

```text

## Error Types

### 1. Test Failures

**Example:**

```

❌ Tests FAILED (3 failures)

```text

#### Resolution:
- Review failed test output in `error_search_results.json`
- Fix the underlying code issues
- Re-run tests to verify fixes


### 2. Syntax Errors

**Example:**

```

❌ Syntax error in framework_investigativo_completo.py: invalid syntax

```text

#### Resolution:
- Open the file and navigate to the reported line
- Fix the syntax error
- Re-run the framework


### 3. Import Errors

**Example:**

```

❌ Missing package: pennylane

```text

**Resolution:**

```bash
pip install -r requirements.txt

```text

### 4. Code Style Violations

**Example:**

```

⚠️  Found 69 linting errors

   - E501: 69 errors

```text

**Resolution:**

```bash

# Auto-fix fixable issues
python -m ruff check . --fix

# Manual fixes for remaining issues
# (See ERROR_SEARCH_REPORT.md for details)

```text

## Integration with Existing Tools

### Ruff

The framework uses Ruff for linting. You can also run Ruff directly:

```bash

# Check for errors
python -m ruff check .

# Auto-fix errors
python -m ruff check . --fix

# Format code
python -m ruff format .

```text

### Pytest

The framework runs pytest tests. You can also run pytest directly:

```bash

# Run all tests
python -m pytest tests/ -v

# Run specific test file
python -m pytest tests/test_repo_smoke.py -v

# Run with coverage
python -m pytest tests/ --cov=. --cov-report=html

```text

## Troubleshooting

### Framework Fails to Run

**Issue:** `ModuleNotFoundError: No module named 'pytest'`


**Solution:**

```bash
pip install pytest ruff

```text

### All Tests Pass But Still Shows Critical Status

**Issue:** Missing dependencies even though tests pass


**Solution:**

```bash

# Ensure all dependencies are installed
pip install -r requirements.txt

```text

### Line Length Violations

**Issue:** Many E501 errors (line too long)


**Solution:**

These are cosmetic issues and don't affect functionality. You can:

1. Manually break long lines
2. Update `.ruff.toml` to increase line length limit
3. Add specific line ignores with `# noqa: E501`


## Best Practices

1. **Run Before Every Commit**

   ```bash
   python error_search_framework.py --fix
   ```text

2. **Include Report in PRs**
   - Attach `ERROR_SEARCH_REPORT.md` to your PR description
   - Shows reviewers you've checked for errors


3. **Monitor Quality Score**
   - Aim for 80+ quality score
   - Address critical issues immediately
   - Clean up warnings over time


4. **Use in Pre-commit Hooks**

   ```bash

   # .git/hooks/pre-commit
   #!/bin/bash
   python error_search_framework.py
   ```

## Configuration

The framework uses configuration from:

- **`.ruff.toml`** - Ruff linting rules
- **`pytest.ini`** or **`pyproject.toml`** - Pytest configuration
- **`requirements.txt`** - List of required dependencies


## Support

For issues or questions about the Error Search Framework:

1. Check the detailed report in `ERROR_SEARCH_REPORT.md`
2. Review the JSON output in `error_search_results.json`
3. Open an issue on GitHub with the error details


## Version History

- **v1.0** (2025-12-23) - Initial release
  - Test execution
  - Syntax checking
  - Import validation
  - Ruff linting integration
  - JSON and Markdown reports
  - Auto-fix support


---


**Last Updated:** 2025-12-23  
**Maintainer:** Framework Team

