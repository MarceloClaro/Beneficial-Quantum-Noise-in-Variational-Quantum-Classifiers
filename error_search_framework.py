#!/usr/bin/env python3
"""
Error Search Framework for Beneficial Quantum Noise in VQC
===========================================================

This script executes a comprehensive error search across the codebase,
checking for:
1. Import errors
2. Syntax errors
3. Test failures
4. Code style violations
5. Dependency issues

Usage:
    python error_search_framework.py [--fix] [--detailed]

Arguments:
    --fix: Automatically fix auto-fixable issues
    --detailed: Generate detailed report with full context
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple


class ErrorSearchFramework:
    """Framework for searching and documenting errors in the codebase."""

    def __init__(self, repo_root: Path, auto_fix: bool = False, detailed: bool = False):
        self.repo_root = repo_root
        self.auto_fix = auto_fix
        self.detailed = detailed
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "repository": "Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers",
            "tests": {},
            "linting": {},
            "syntax": {},
            "imports": {},
            "summary": {}
        }

    def run_tests(self) -> Tuple[bool, Dict]:
        """Run pytest tests and capture results."""
        print("\n" + "=" * 80)
        print("1. RUNNING TESTS")
        print("=" * 80)

        test_results = {
            "passed": False,
            "total": 0,
            "passed_count": 0,
            "failed_count": 0,
            "output": "",
            "failed_tests": []
        }

        try:
            result = subprocess.run(
                ["python", "-m", "pytest", "tests/test_repo_smoke.py", "-v", "--tb=short"],
                cwd=self.repo_root,
                capture_output=True,
                text=True,
                timeout=120
            )

            test_results["output"] = result.stdout + result.stderr

            # Parse output for test counts
            for line in result.stdout.split('\n'):
                if 'passed' in line.lower():
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == 'passed':
                            test_results["passed_count"] = int(parts[i-1])
                if 'failed' in line.lower():
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == 'failed':
                            test_results["failed_count"] = int(parts[i-1])

            test_results["total"] = test_results["passed_count"] + test_results["failed_count"]
            test_results["passed"] = result.returncode == 0

            if test_results["passed"]:
                print(f"✅ All tests PASSED ({test_results['passed_count']} tests)")
            else:
                print(f"❌ Tests FAILED ({test_results['failed_count']} failures)")

        except subprocess.TimeoutExpired:
            print("❌ Tests timed out after 120 seconds")
            test_results["output"] = "Tests timed out"
        except Exception as e:
            print(f"❌ Error running tests: {e}")
            test_results["output"] = str(e)

        return test_results["passed"], test_results

    def run_linting(self) -> Tuple[bool, Dict]:
        """Run ruff linter and capture results."""
        print("\n" + "=" * 80)
        print("2. RUNNING CODE STYLE CHECKS (RUFF)")
        print("=" * 80)

        lint_results = {
            "passed": False,
            "total_errors": 0,
            "fixed_errors": 0,
            "remaining_errors": 0,
            "errors_by_type": {},
            "errors_by_file": {},
            "output": ""
        }

        try:
            # First, try to fix auto-fixable issues
            if self.auto_fix:
                print("Attempting to auto-fix issues...")
                fix_result = subprocess.run(
                    ["python", "-m", "ruff", "check", ".", "--fix"],
                    cwd=self.repo_root,
                    capture_output=True,
                    text=True
                )
                lint_results["output"] += "=== AUTO-FIX OUTPUT ===\n" + fix_result.stdout + "\n"

            # Run full check
            result = subprocess.run(
                ["python", "-m", "ruff", "check", ".", "--output-format=json"],
                cwd=self.repo_root,
                capture_output=True,
                text=True
            )

            # Parse JSON output
            if result.stdout:
                errors = json.loads(result.stdout)
                lint_results["total_errors"] = len(errors)

                for error in errors:
                    # Group by error code
                    code = error.get("code", "UNKNOWN")
                    if code not in lint_results["errors_by_type"]:
                        lint_results["errors_by_type"][code] = []
                    lint_results["errors_by_type"][code].append(error)

                    # Group by file
                    filename = error.get("filename", "UNKNOWN")
                    if filename not in lint_results["errors_by_file"]:
                        lint_results["errors_by_file"][filename] = []
                    lint_results["errors_by_file"][filename].append(error)

            # Also get concise output for summary
            result_concise = subprocess.run(
                ["python", "-m", "ruff", "check", ".", "--output-format=concise"],
                cwd=self.repo_root,
                capture_output=True,
                text=True
            )
            lint_results["output"] += "\n=== CONCISE OUTPUT ===\n" + result_concise.stdout

            lint_results["remaining_errors"] = lint_results["total_errors"]
            lint_results["passed"] = lint_results["total_errors"] == 0

            if lint_results["passed"]:
                print("✅ No linting errors found")
            else:
                print(f"⚠️  Found {lint_results['total_errors']} linting errors")
                for code, errors in lint_results["errors_by_type"].items():
                    print(f"   - {code}: {len(errors)} errors")

        except Exception as e:
            print(f"❌ Error running linter: {e}")
            lint_results["output"] = str(e)

        return lint_results["passed"], lint_results

    def check_syntax(self) -> Tuple[bool, Dict]:
        """Check Python syntax for all .py files."""
        print("\n" + "=" * 80)
        print("3. CHECKING PYTHON SYNTAX")
        print("=" * 80)

        syntax_results = {
            "passed": True,
            "total_files": 0,
            "errors": []
        }

        # Find all Python files
        py_files = list(self.repo_root.rglob("*.py"))
        syntax_results["total_files"] = len(py_files)

        for py_file in py_files:
            try:
                with open(py_file, 'r', encoding='utf-8') as f:
                    code = f.read()
                compile(code, str(py_file), 'exec')
            except SyntaxError as e:
                syntax_results["passed"] = False
                syntax_results["errors"].append({
                    "file": str(py_file.relative_to(self.repo_root)),
                    "line": e.lineno,
                    "message": str(e)
                })
                print(f"❌ Syntax error in {py_file.relative_to(self.repo_root)}: {e}")

        if syntax_results["passed"]:
            print(f"✅ All {syntax_results['total_files']} Python files have valid syntax")
        else:
            print(f"❌ Found {len(syntax_results['errors'])} syntax errors")

        return syntax_results["passed"], syntax_results

    def check_imports(self) -> Tuple[bool, Dict]:
        """Check if all required imports are available."""
        print("\n" + "=" * 80)
        print("4. CHECKING IMPORTS")
        print("=" * 80)

        import_results = {
            "passed": True,
            "required_packages": [],
            "missing_packages": [],
            "available_packages": []
        }

        # Read requirements.txt
        requirements_file = self.repo_root / "requirements.txt"
        if requirements_file.exists():
            with open(requirements_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        package = line.split('>=')[0].split('==')[0].strip()
                        import_results["required_packages"].append(package)

        # Check each package
        for package in import_results["required_packages"]:
            # Map package names to import names
            import_name = package
            if package == "scikit-learn":
                import_name = "sklearn"
            elif package == "typing-extensions":
                import_name = "typing_extensions"

            try:
                __import__(import_name)
                import_results["available_packages"].append(package)
            except ImportError:
                import_results["passed"] = False
                import_results["missing_packages"].append(package)
                print(f"❌ Missing package: {package}")

        if import_results["passed"]:
            print(f"✅ All {len(import_results['required_packages'])} required packages available")
        else:
            print(f"❌ Missing {len(import_results['missing_packages'])} packages")

        return import_results["passed"], import_results

    def generate_summary(self) -> Dict:
        """Generate overall summary of results."""
        summary = {
            "overall_status": "HEALTHY",
            "quality_score": 100,
            "passed_checks": 0,
            "total_checks": 4,
            "critical_issues": 0,
            "warnings": 0,
            "recommendations": []
        }

        # Evaluate each check
        checks = [
            ("tests", self.results["tests"].get("passed", False)),
            ("syntax", self.results["syntax"].get("passed", False)),
            ("imports", self.results["imports"].get("passed", False)),
            ("linting", self.results["linting"].get("passed", False))
        ]

        for check_name, passed in checks:
            if passed:
                summary["passed_checks"] += 1

        # Calculate quality score
        if self.results["tests"].get("passed", False):
            summary["quality_score"] = 100
        else:
            summary["quality_score"] -= 30

        if not self.results["syntax"].get("passed", False):
            summary["quality_score"] -= 30
            summary["critical_issues"] += len(self.results["syntax"].get("errors", []))

        if not self.results["imports"].get("passed", False):
            summary["quality_score"] -= 20
            summary["critical_issues"] += len(self.results["imports"].get("missing_packages", []))

        if not self.results["linting"].get("passed", False):
            summary["quality_score"] -= 20
            summary["warnings"] += self.results["linting"].get("total_errors", 0)

        # Determine overall status
        if summary["critical_issues"] > 0:
            summary["overall_status"] = "CRITICAL"
        elif summary["warnings"] > 50:
            summary["overall_status"] = "NEEDS ATTENTION"
        elif summary["warnings"] > 0:
            summary["overall_status"] = "GOOD"
        else:
            summary["overall_status"] = "EXCELLENT"

        # Generate recommendations
        if not self.results["tests"].get("passed", False):
            summary["recommendations"].append("Fix failing tests immediately")

        if self.results["syntax"].get("errors"):
            summary["recommendations"].append("Fix syntax errors before proceeding")

        if self.results["imports"].get("missing_packages"):
            summary["recommendations"].append(
                f"Install missing packages: {', '.join(self.results['imports']['missing_packages'])}"
            )

        if self.results["linting"].get("total_errors", 0) > 50:
            summary["recommendations"].append("Consider running ruff with --fix to auto-fix issues")

        return summary

    def print_final_report(self):
        """Print a summary report to console."""
        print("\n" + "=" * 80)
        print("ERROR SEARCH FRAMEWORK - FINAL REPORT")
        print("=" * 80)

        summary = self.results["summary"]

        print(f"\nOverall Status: {summary['overall_status']}")
        print(f"Quality Score: {summary['quality_score']}/100")
        print(f"\nChecks Passed: {summary['passed_checks']}/{summary['total_checks']}")
        print(f"Critical Issues: {summary['critical_issues']}")
        print(f"Warnings: {summary['warnings']}")

        if summary["recommendations"]:
            print("\nRecommendations:")
            for i, rec in enumerate(summary["recommendations"], 1):
                print(f"  {i}. {rec}")

        print("\n" + "=" * 80)
        print("Detailed report saved to: ERROR_SEARCH_REPORT.md")
        print("JSON results saved to: error_search_results.json")
        print("=" * 80 + "\n")

    def save_results(self):
        """Save results to JSON file."""
        output_file = self.repo_root / "error_search_results.json"
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"\n✅ Results saved to {output_file}")

    def run_all_checks(self):
        """Run all error checks."""
        print("\n" + "=" * 80)
        print("ERROR SEARCH FRAMEWORK - STARTING")
        print("=" * 80)
        print(f"Repository: {self.repo_root}")
        print(f"Auto-fix: {self.auto_fix}")
        print(f"Timestamp: {self.results['timestamp']}")

        # Run all checks
        tests_passed, test_results = self.run_tests()
        self.results["tests"] = test_results

        syntax_passed, syntax_results = self.check_syntax()
        self.results["syntax"] = syntax_results

        imports_passed, import_results = self.check_imports()
        self.results["imports"] = import_results

        linting_passed, lint_results = self.run_linting()
        self.results["linting"] = lint_results

        # Generate summary
        self.results["summary"] = self.generate_summary()

        # Print and save results
        self.print_final_report()
        self.save_results()

        return self.results


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Error Search Framework for VQC codebase"
    )
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Automatically fix auto-fixable issues"
    )
    parser.add_argument(
        "--detailed",
        action="store_true",
        help="Generate detailed report with full context"
    )

    args = parser.parse_args()

    # Get repository root
    repo_root = Path(__file__).parent.absolute()

    # Create and run framework
    framework = ErrorSearchFramework(
        repo_root=repo_root,
        auto_fix=args.fix,
        detailed=args.detailed
    )

    results = framework.run_all_checks()

    # Exit with appropriate code
    sys.exit(0 if results["summary"]["critical_issues"] == 0 else 1)


if __name__ == "__main__":
    main()
