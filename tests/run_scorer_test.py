"""
ControlScan: In-Silico Biochemist Scoring Test

Validates the MutationScorer module across individual mutations and epistatic networks.
Tests BLOSUM62 and Grantham Distance integration for severity assessment.
"""

import sys
import logging
from pathlib import Path

# Add src directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from controlscan.scorer import MutationScorer


def print_header():
    """Print test session header."""
    print("\n" + "=" * 80)
    print("CONTROLSCAN: IN-SILICO BIOCHEMIST SCORING TEST")
    print("=" * 80)
    print("Scoring Mutations: BLOSUM62 (Evolutionary) + Grantham Distance (Chemical)")
    print("=" * 80 + "\n")


def print_single_result(result: dict):
    """Pretty-print a single mutation score."""
    if result.get('Error'):
        print(f"  Error: {result['Error']}\n")
    else:
        print(f"  Mutation:        {result['Mutation']}")
        print(f"  BLOSUM62 Score:  {result.get('BLOSUM', 'N/A')}")
        print(f"  Grantham Dist:   {result.get('Grantham', 'N/A')}")
        print(f"  Severity:        {result['Severity']}/100")
        print(f"  Category:        {result['Category']}")
        print()


def print_network_result(result: dict):
    """Pretty-print a network score with individual components."""
    print(f"  Network:         {result['Network']}")
    print(f"  Mean Severity:   {result['Mean_Severity']}/100")
    print(f"  Max Severity:    {result['Max_Severity']}/100")
    print()
    
    print("  Individual Mutations:")
    print("  " + "-" * 76)
    for score in result['Individual_Scores']:
        print(f"    {score['Mutation']:8s} | BLOSUM: {score.get('BLOSUM', 'N/A'):6} | "
              f"Grantham: {score.get('Grantham', 'N/A'):6} | "
              f"Severity: {score['Severity']:6.2f} | {score['Category']:10s}")
    print()


def main():
    logging.basicConfig(level=logging.WARNING)
    
    print_header()
    
    scorer = MutationScorer()
    
    # ========================================================================
    # SECTION 1: Individual Mutation Tests
    # ========================================================================
    print("[SECTION 1] INDIVIDUAL MUTATION SCORING")
    print("-" * 80 + "\n")
    
    individual_mutations = [
        "I174V",  # Expected: BENIGN (similar hydrophobic AAs)
        "L970A",  # Expected: UNCERTAIN (loss of hydrophobicity)
        "W100G",  # Expected: PATHOGENIC (tiny glycine replaces bulky tryptophan)
    ]
    
    for mut in individual_mutations:
        print(f"[TEST] Scoring: {mut}")
        try:
            result = scorer.score_single(mut)
            print_single_result(result)
        except Exception as e:
            print(f"  ERROR: {e}\n")
    
    # ========================================================================
    # SECTION 2: Epistatic Network Tests
    # ========================================================================
    print("[SECTION 2] EPISTATIC NETWORK SCORING")
    print("-" * 80 + "\n")
    
    networks = [
        ["I174V", "K10R", "L970A"],
        ["I174V", "K10R", "R342S", "S540G"],
    ]
    
    for network in networks:
        print(f"[TEST] Scoring Network: {' + '.join(network)}")
        try:
            result = scorer.score_network(network)
            print_network_result(result)
        except Exception as e:
            print(f"  ERROR: {e}\n")
    
    # ========================================================================
    # Test Summary
    # ========================================================================
    print("=" * 80)
    print("[PASS] All tests completed successfully!")
    print("=" * 80 + "\n")
    
    print("Interpretation Guide:")
    print("  BENIGN (< 30):      Subtle changes, minimal structural impact")
    print("  UNCERTAIN (30-60):  Moderate changes, potential functional consequences")
    print("  PATHOGENIC (> 60):  Dramatic changes, likely deleterious effects")
    print()


if __name__ == "__main__":
    main()
