# ControlScan: Systems Biology Microservice for Epistasis Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**ControlScan** is a standalone Python microservice designed for computational systems biology. It processes transcriptomic networks and calculates genomic co-mutation probabilities (epistasis analysis).

The tool is built for:
- **Transcriptomic analysis**: Expression normalization and correlation-based network construction
- **Mutation epistasis**: Co-occurrence analysis and mutual exclusivity/inclusivity scoring
- **Data integration**: Seamless integration with NCBI Gene Expression Omnibus (GEO) databases

## Features

### 🧬 Core Modules

1. **Epistasis Engine** (`epistasis.py`)
   - Co-occurrence matrix construction from mutation samples
   - Jaccard Index and conditional probability calculations
   - Mutual exclusivity/inclusivity scoring
   - Mutation frequency analysis

2. **Data Fetcher** (`fetcher.py`)
   - Download expression matrices from NCBI GEO
   - Parse SOFT files and count matrices
   - Sample metadata extraction

3. **Expression Normalizer** (`normalizer.py`)
   - TPM (Transcripts Per Million) normalization
   - Log transformation and z-score normalization
   - Quantile normalization for batch correction
   - Comprehensive normalization statistics

4. **Co-expression Network** (`network.py`)
   - Pearson, Spearman, and Kendall correlation analysis
   - Complete network construction with statistical filtering
   - Regulatory candidate identification
   - Gene correlation heatmaps

## Installation

### Requirements
- Python 3.8+
- NumPy, Pandas, SciPy
- GEOparse (for GEO data download)

### Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/ControlScan.git
cd ControlScan

# Install dependencies
pip install -r requirements.txt

# Install the package (development mode)
pip install -e .
```

## Quick Start

### 1. Download Expression Data from GEO

```python
from controlscan import GEOFetcher, ExpressionNormalizer

# Fetch data
fetcher = GEOFetcher()
matrix = fetcher.download_matrix('GSE12345')

# Normalize to TPM
normalizer = ExpressionNormalizer(matrix)
tpm = normalizer.normalize_tpm(log_transform=True)
```

### 2. Build Co-expression Network

```python
from controlscan import NetworkBuilder

# Create network
builder = NetworkBuilder(tpm)
correlations = builder.calculate_correlation('BRCA1')
print(correlations.head(10))

# Get regulatory candidates
regulators = builder.get_regulatory_candidates('BRCA1', top_n=15)
```

### 3. Epistasis Analysis

```python
from controlscan import MutationMatrix

# Load mutation data
samples = [
    {'Sample_ID': 'S1', 'Mutations': ['acrA_T104A', 'acrB_N596H']},
    {'Sample_ID': 'S2', 'Mutations': ['acrA_T104A', 'rpoB_H526Y']},
    {'Sample_ID': 'S3', 'Mutations': ['acrB_N596H', 'rpoB_H526Y']},
]

# Build co-occurrence matrix
mutation_analyzer = MutationMatrix(samples)
co_matrix = mutation_analyzer.build_co_occurrence_matrix()

# Calculate exclusivity scores
score = mutation_analyzer.calculate_exclusivity('acrA_T104A', 'acrB_N596H')
print(f"Jaccard Index: {score:.3f}")  # 1.0 = mutually inclusive, 0.0 = mutually exclusive

# Conditional probability
prob = mutation_analyzer.calculate_conditional_probability('acrA_T104A', 'acrB_N596H')
print(f"P(acrA_T104A | acrB_N596H) = {prob:.3f}")
```

## API Documentation

### MutationMatrix

**Methods:**

- `build_co_occurrence_matrix()` → pd.DataFrame
  - Builds N×N co-occurrence matrix
  - Returns symmetric DataFrame of mutation co-occurrence counts

- `calculate_exclusivity(mutation_a, mutation_b)` → float
  - Calculates Jaccard Index (mutual exclusivity score)
  - Returns score in [0, 1]

- `calculate_conditional_probability(mutation_given, mutation_conditional)` → float
  - Calculates P(A | B)
  - Returns probability in [0, 1]

- `get_mutation_summary()` → pd.DataFrame
  - Returns mutation frequency statistics

### ExpressionNormalizer

**Methods:**

- `normalize_tpm(expression_matrix, gene_lengths, log_transform)` → pd.DataFrame
  - Normalize to Transcripts Per Million
  - Accounts for gene length and sequencing depth

- `normalize_z_score()` → pd.DataFrame
  - Z-score normalization per gene
  - Produces mean=0, std=1 distribution

- `normalize_quantile()` → pd.DataFrame
  - Quantile normalization across samples
  - Removes technical batch effects

### NetworkBuilder

**Methods:**

- `calculate_correlation(target_gene, method, min_samples)` → pd.DataFrame
  - Correlate target gene against all genes
  - Methods: 'pearson', 'spearman', 'kendall'
  - Returns correlation, p-value, and significance

- `build_correlation_network(correlation_threshold, pvalue_threshold)` → pd.DataFrame
  - Complete pairwise correlation network
  - Returns edge list

- `get_regulatory_candidates(target_gene, top_n)` → pd.DataFrame
  - Top correlated genes (potential regulators)

## Data Formats

### Mutation Data
```python
# List of dictionaries
samples = [
    {'Sample_ID': 'S1', 'Mutations': ['mutA', 'mutB']},
    {'Sample_ID': 'S2', 'Mutations': ['mutA', 'mutC']},
]

# Or pandas DataFrame
df = pd.DataFrame({
    'Sample_ID': ['S1', 'S2'],
    'Mutations': [['mutA', 'mutB'], ['mutA', 'mutC']]
})
```

### Expression Data
```python
# Genes × Samples DataFrame
expression_matrix = pd.DataFrame(
    data=[[1.2, 3.4, 2.1], [4.5, 2.3, 5.6]],
    index=['GENE1', 'GENE2'],      # Gene identifiers (rows)
    columns=['Sample1', 'Sample2', 'Sample3']  # Sample IDs (columns)
)
```

## Testing

Run the test suite:

```bash
# Run all tests
python -m pytest tests/

# Run specific test file
python -m pytest tests/test_epistasis.py -v

# Run with coverage
python -m pytest tests/ --cov=controlscan
```

## Project Structure

```
ControlScan/
├── src/
│   └── controlscan/
│       ├── __init__.py           # Package initialization
│       ├── epistasis.py          # Co-mutation analysis
│       ├── fetcher.py            # GEO data fetching
│       ├── normalizer.py         # Expression normalization
│       └── network.py            # Co-expression networks
├── tests/
│   ├── test_epistasis.py         # Epistasis tests
│   └── test_fetcher.py           # Fetcher tests
├── data/
│   ├── raw/                      # Raw input data
│   └── processed/                # Processed outputs
├── requirements.txt              # Python dependencies
└── README.md                     # This file
```

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| NumPy | ≥1.20.0 | Numerical computing |
| Pandas | ≥1.3.0 | Data manipulation |
| SciPy | ≥1.7.0 | Statistical functions |
| scikit-learn | ≥0.24.0 | Machine learning utilities |
| GEOparse | ≥2.0.0 | NCBI GEO data parsing |

## Logging

All modules include comprehensive logging. Configure logging in your application:

```python
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
```

## Performance Notes

- **Co-occurrence matrices**: O(n_samples × m_mutations²) - scales well up to 100K samples
- **Correlation calculations**: O(n_genes × m_samples) - typical genomic data (~20K genes, 100+ samples) processes in seconds
- **Network construction**: Memory-intensive for >50K genes; use gene filtering for large datasets

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit changes (`git commit -m 'Add AmazingFeature'`)
4. Push to branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

MIT License - See LICENSE file for details

## Citation

If you use ControlScan in your research, please cite:

```bibtex
@software{controlscan2025,
  title={ControlScan: Systems Biology Microservice for Epistasis Analysis},
  author={Your Name},
  year={2025},
  url={https://github.com/yourusername/ControlScan}
}
```

## References

- **GEO Database**: https://www.ncbi.nlm.nih.gov/geo/
- **GEOparse Documentation**: https://geoparse.readthedocs.io/
- **Pearson Correlation**: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
- **Jaccard Index**: https://en.wikipedia.org/wiki/Jaccard_index
- **TPM Normalization**: Wagner et al. (2012) Theory in Biosciences

## Authors

- Senior Computational Biologist & Python Software Architect

## Acknowledgments

Built with standardized bioinformatics practices for reproducible, scalable systems biology research.

---

**Last Updated**: February 2025
**Version**: 0.1.0
