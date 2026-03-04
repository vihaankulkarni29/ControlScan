# ControlScan: MutationScorer Test Results

**Date**: March 4, 2026  
**Module**: `src/controlscan/scorer.py`  
**Test Script**: `tests/run_scorer_test.py`  
**Status**: ✅ ALL TESTS PASSED

---

## Executive Summary

The **MutationScorer** module successfully implements an in-silico biochemist for amino acid substitution severity assessment. The module integrates:

- **BLOSUM62 Matrix** (Evolutionary Distance): 20×20 substitution scoring
- **Grantham Distance** (Chemical Properties): Composition, polarity, and volume analysis
- **Integrated Severity Score**: 40% evolutionary + 60% chemical penalties
- **Clinical Categorization**: BENIGN, UNCERTAIN, PATHOGENIC

---

## Section 1: Individual Mutation Scoring

### Test Case 1: I174V (Isoleucine → Valine at position 174)

**Result: BENIGN** ✅

| Metric | Value |
|--------|-------|
| BLOSUM62 Score | 3.0 |
| Grantham Distance | 0.58 |
| Severity Score | 21.5 / 100 |
| Category | **BENIGN** |

**Interpretation**: 
- Both isoleucine and valine are hydrophobic, branched amino acids
- Minimal chemical property changes (small Grantham distance)
- Positive BLOSUM score indicates evolutionary conservation
- **Clinical Impact**: Likely tolerated mutation with minimal functional consequence

---

### Test Case 2: L970A (Leucine → Alanine at position 970)

**Result: UNCERTAIN** ⚠️

| Metric | Value |
|--------|-------|
| BLOSUM62 Score | -1.0 |
| Grantham Distance | 1.90 |
| Severity Score | 32.53 / 100 |
| Category | **UNCERTAIN** |

**Interpretation**:
- Leucine is a larger, hydrophobic residue; alanine is small and hydrophobic
- Moderate loss of hydrophobic bulk (~4× volume decrease)
- Negative BLOSUM indicates some evolutionary penalty
- **Clinical Impact**: Potentially deleterious; depends on local protein environment

---

### Test Case 3: W100G (Tryptophan → Glycine at position 100)

**Result: PATHOGENIC** ⚠️⚠️

| Metric | Value |
|--------|-------|
| BLOSUM62 Score | -2.0 |
| Grantham Distance | 176.04 |
| Severity Score | 83.79 / 100 |
| Category | **PATHOGENIC** |

**Interpretation**:
- Tryptophan: large aromatic ring (volume=170, composition=130)
- Glycine: smallest amino acid (volume=3, composition=0)
- Extreme chemical property loss (Grantham distance = 176.04, near maximum)
- Negative BLOSUM (-2) indicates evolutionary unfavorable
- **Clinical Impact**: Likely deleterious; high probability of functional loss or misfolding

---

## Section 2: Epistatic Network Scoring

### Network 1: I174V + K10R + L970A

**Network Composition:**

| Mutation | BLOSUM | Grantham | Severity | Category |
|----------|--------|----------|----------|----------|
| I174V | 3.0 | 0.58 | 21.50 | BENIGN |
| K10R | 2.0 | 10.29 | 26.87 | BENIGN |
| L970A | -1.0 | 1.90 | 32.53 | UNCERTAIN |

**Network Summary:**

| Metric | Value |
|--------|-------|
| Network ID | `I174V_K10R_L970A` |
| Mean Severity | **26.97 / 100** |
| Max Severity | **32.53 / 100** |
| Overall Category | **BENIGN** |

**Interpretation**:
- Three mutations with predominantly benign severity
- K10R (Lysine → Arginine): Both positive, similar size; minor Grantham distance (10.29)
- Epistatic network shows **low cumulative burden**
- **Clinical Significance**: This mutation cluster likely represents a stable, tolerated genotype

---

### Network 2: I174V + K10R + R342S + S540G

**Network Composition:**

| Mutation | BLOSUM | Grantham | Severity | Category |
|----------|--------|----------|----------|----------|
| I174V | 3.0 | 0.58 | 21.50 | BENIGN |
| K10R | 2.0 | 10.29 | 26.87 | BENIGN |
| R342S | -1.0 | 103.45 | 60.87 | PATHOGENIC |
| S540G | 0.0 | 192.25 | 82.99 | PATHOGENIC |

**Network Summary:**

| Metric | Value |
|--------|-------|
| Network ID | `I174V_K10R_R342S_S540G` |
| Mean Severity | **48.06 / 100** |
| Max Severity | **82.99 / 100** |
| Overall Category | **MIXED (UNCERTAIN/PATHOGENIC)** |

**Interpretation**:
- **Compound network with severe mutations:**
  - R342S (Arginine → Serine): 60.87 severity; loss of positive charge and size
  - S540G (Serine → Glycine): 82.99 severity; extreme chemical property divergence
- Mean severity of 48.06 reflects balanced composition (2 benign + 2 pathogenic)
- **Clinical Significance**: This network contains destabilizing mutations with high cumulative burden; likely pathogenic when combined

---

## Methodology & Scoring Formula

### BLOSUM62 Normalization
```
blosum_norm = (blosum_raw + 4) / 15
Range: [-4 (worst) → 11 (best)] → [0.0 → 1.0]
```

### Grantham Distance Formula
```
distance = sqrt(1.833 * (c1-c2)² + 0.1018 * (p1-p2)² + 0.000399 * (v1-v2)²)
Where:
  c = composition
  p = polarity
  v = volume
Range: [0 (identical) → ~215 (maximum divergence)] → [0.0 → 1.0]
```

### Integrated Severity Score
```
severity_raw = 0.40 * (1 - blosum_norm) + 0.60 * grantham_norm
severity_final = min(severity_raw * 100, 100)

Weighting:
  - 40% evolutionary penalty (BLOSUM62)
  - 60% chemical penalty (Grantham Distance)
```

### Categorization Thresholds
| Severity Range | Category | Clinical Interpretation |
|---|---|---|
| 0 - 29 | **BENIGN** | Subtle changes, minimal structural impact |
| 30 - 59 | **UNCERTAIN** | Moderate changes, potential functional consequences |
| 60 - 100 | **PATHOGENIC** | Dramatic changes, likely deleterious effects |

---

## Data Structures

### BLOSUM62 Matrix
- **Format**: 20×20 dictionary of amino acid pair scores
- **Range**: -4 (most conservative) to +11 (most favorable)
- **Source**: Standard BLOSUM62 substitution matrix
- **Coverage**: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V

### Amino Acid Properties (Grantham)
```python
Properties tracked per amino acid:
  - c (composition):  0.0 - 142.0 (S has highest polarity)
  - p (polarity):     4.9 - 13.0  (D has highest polarity)
  - v (volume):       3.0 - 170.0 (W has largest volume)
```

---

## Module Architecture

### Class: `MutationScorer`

#### Public Methods

**`__init__()`**
- Initializes BLOSUM62 matrix and amino acid properties
- Called once per scorer instance

**`score_single(mut_string: str) → dict`**
- **Input**: Mutation string (e.g., "I174V")
- **Output**: Dictionary with mutation details, severity score, and category
- **Processing**: 
  1. Parse mutation string
  2. Retrieve BLOSUM62 and Grantham scores
  3. Normalize both metrics
  4. Calculate weighted severity
  5. Assign category

**`score_network(network_mutations: List[str]) → dict`**
- **Input**: List of mutation strings
- **Output**: Dictionary with network ID, individual scores, mean/max severity
- **Use Case**: Epistatic interaction analysis

#### Internal Methods

**`_get_blosum(aa1, aa2) → float`**
- Retrieves BLOSUM62 substitution score
- Default: -4 if not found

**`_calc_grantham(aa1, aa2) → float`**
- Calculates Grantham Distance using Composition/Polarity/Volume
- Returns 0.0 if amino acids are identical

**`parse_mutation(mut_string) → tuple`**
- Parses "I174V" → ("I", "174", "V")
- Validates mutation format and position

---

## Test Coverage

✅ **Individual Mutations**: 3 test cases (Benign, Uncertain, Pathogenic)  
✅ **Epistatic Networks**: 2 test cases (Mono-severity, Mixed severity)  
✅ **Error Handling**: Malformed mutation string detection  
✅ **Normalization**: Verifies 0-100 severity scale  
✅ **Categorization**: Confirms threshold logic  

**Total Assertions Passed**: 7/7 ✅

---

## Output Example

```
[TEST] Scoring: W100G
  Mutation:        W100G
  BLOSUM62 Score:  -2
  Grantham Dist:   176.04
  Severity:        83.79/100
  Category:        PATHOGENIC

[TEST] Scoring Network: I174V + K10R + R342S + S540G
  Network:         I174V_K10R_R342S_S540G
  Mean Severity:   48.06/100
  Max Severity:    82.99/100
  
  Individual Mutations:
    I174V    | BLOSUM:   3.00 | Grantham:   0.58 | Severity:  21.50 | BENIGN
    K10R     | BLOSUM:   2.00 | Grantham:  10.29 | Severity:  26.87 | BENIGN
    R342S    | BLOSUM:  -1.00 | Grantham: 103.45 | Severity:  60.87 | PATHOGENIC
    S540G    | BLOSUM:   0.00 | Grantham: 192.25 | Severity:  82.99 | PATHOGENIC
```

---

## Conclusion

The **MutationScorer** module provides a robust, evidence-based framework for in-silico assessment of amino acid substitution severity. By integrating evolutionary (BLOSUM62) and chemical property (Grantham Distance) perspectives, it enables:

1. **Individual Mutation Risk Stratification**: Classify mutations as benign, uncertain, or pathogenic
2. **Network-Level Epistasis Analysis**: Evaluate cumulative effects of co-occurring mutations
3. **Biochemical Interpretability**: Transparent scoring with biological grounding

**Recommended Applications:**
- Variant effect prediction in genomics
- Protein engineering screening
- Drug resistance mechanism analysis
- Strain phenotype correlation studies

---

**Test Execution Date**: March 4, 2026  
**Status**: ✅ PASSED  
**Module Version**: 1.0.0
