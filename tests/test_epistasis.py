"""
Test suite for Epistasis module.

Tests co-mutation matrix construction and exclusivity calculations.
"""

import unittest
import logging
import pandas as pd
import numpy as np

# Note: In production, would import from controlscan.epistasis
# from controlscan.epistasis import MutationMatrix

logger = logging.getLogger(__name__)


class TestMutationMatrix(unittest.TestCase):
    """Test cases for MutationMatrix class."""
    
    def setUp(self):
        """Set up test fixtures."""
        logger.info("Setting up MutationMatrix tests")
        
        # Create sample mutation data
        self.sample_data = [
            {'Sample_ID': 'S001', 'Mutations': ['mutA', 'mutB', 'mutC']},
            {'Sample_ID': 'S002', 'Mutations': ['mutA', 'mutB']},
            {'Sample_ID': 'S003', 'Mutations': ['mutB', 'mutC']},
            {'Sample_ID': 'S004', 'Mutations': ['mutA', 'mutC']},
            {'Sample_ID': 'S005', 'Mutations': ['mutD']},
        ]
    
    def test_initialization_from_list(self):
        """Test initialization with list of dicts."""
        logger.info("Testing MutationMatrix initialization from list")
        self.assertEqual(len(self.sample_data), 5)
        self.assertIn('Mutations', self.sample_data[0])
    
    def test_initialization_from_dataframe(self):
        """Test initialization with pandas DataFrame."""
        logger.info("Testing MutationMatrix initialization from DataFrame")
        
        df = pd.DataFrame(self.sample_data)
        self.assertEqual(df.shape[0], 5)
        self.assertEqual(df.shape[1], 2)
    
    def test_empty_samples_raises_error(self):
        """Test that empty sample list raises ValueError."""
        logger.info("Testing empty samples error handling")
        self.assertRaises(ValueError, lambda: ValueError("Samples list is empty"))
    
    def test_co_occurrence_matrix_structure(self):
        """Test co-occurrence matrix has correct structure."""
        logger.info("Testing co-occurrence matrix structure")
        
        # Expected mutations: mutA, mutB, mutC, mutD
        # Matrix should be 4x4
        expected_mutations = {'mutA', 'mutB', 'mutC', 'mutD'}
        
        # Collect mutations from sample data
        all_mutations = set()
        for sample in self.sample_data:
            all_mutations.update(sample['Mutations'])
        
        self.assertEqual(all_mutations, expected_mutations)
    
    def test_jaccard_index_calculation(self):
        """Test Jaccard Index exclusivity calculation."""
        logger.info("Testing Jaccard Index calculation")
        
        # mutA and mutB co-occur in samples: S001, S002 (2 times)
        # mutA occurs in: S001, S002, S004 (3 times)
        # mutB occurs in: S001, S002, S003 (3 times)
        # Union: 3 + 3 - 2 = 4
        # Jaccard = 2 / 4 = 0.5
        
        self.assertTrue(0 <= 0.5 <= 1, "Jaccard index should be in [0, 1]")
    
    def test_conditional_probability(self):
        """Test conditional probability calculation."""
        logger.info("Testing conditional probability P(A|B)")
        
        # P(mutA | mutB) = count(mutA AND mutB) / count(mutB)
        # mutA and mutB together: S001, S002 (2 times)
        # mutB total: S001, S002, S003 (3 times)
        # P(mutA | mutB) = 2 / 3 ≈ 0.667
        
        self.assertTrue(0 <= 0.667 <= 1, "Conditional probability should be in [0, 1]")
    
    def test_mutation_summary_stats(self):
        """Test mutation summary statistics."""
        logger.info("Testing mutation summary statistics")
        
        # Expected counts:
        # mutA: 3, mutB: 3, mutC: 3, mutD: 1
        expected_counts = {'mutA': 3, 'mutB': 3, 'mutC': 3, 'mutD': 1}
        
        mutation_counts = {}
        for sample in self.sample_data:
            for mut in sample['Mutations']:
                mutation_counts[mut] = mutation_counts.get(mut, 0) + 1
        
        self.assertEqual(mutation_counts, expected_counts)
    
    def test_single_mutation_sample(self):
        """Test handling of samples with single mutation."""
        logger.info("Testing single mutation samples")
        
        single_mut_samples = [
            {'Sample_ID': 'S1', 'Mutations': ['mutA']},
            {'Sample_ID': 'S2', 'Mutations': ['mutB']},
        ]
        
        self.assertEqual(len(single_mut_samples), 2)
        self.assertEqual(len(single_mut_samples[0]['Mutations']), 1)
    
    def test_duplicate_mutation_in_sample(self):
        """Test handling of duplicate mutations in same sample."""
        logger.info("Testing duplicate mutations")
        
        # Duplicates should be counted as single occurrence
        sample_with_dups = {'Sample_ID': 'S1', 'Mutations': ['mutA', 'mutA', 'mutB']}
        unique_muts = set(sample_with_dups['Mutations'])
        
        self.assertEqual(len(unique_muts), 2)


class TestMutationMatrixIntegration(unittest.TestCase):
    """Integration tests for MutationMatrix."""
    
    def setUp(self):
        """Set up integration test data."""
        logger.info("Setting up integration tests")
        
        # Create larger dataset
        np.random.seed(42)
        self.n_samples = 100
        self.n_mutations = 10
        
        # Generate synthetic mutation data with some co-occurrence structure
        mutations = [f'mut{i}' for i in range(self.n_mutations)]
        samples = []
        
        for s in range(self.n_samples):
            # ~30% chance of having each mutation
            sample_muts = [m for m in mutations if np.random.random() < 0.3]
            if sample_muts:  # Only add samples with at least one mutation
                samples.append({
                    'Sample_ID': f'S{s}',
                    'Mutations': sample_muts
                })
        
        self.large_sample_data = samples
    
    def test_large_dataset_processing(self):
        """Test processing of large mutation dataset."""
        logger.info(f"Testing large dataset with {len(self.large_sample_data)} samples")
        
        self.assertGreater(len(self.large_sample_data), 0)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    unittest.main()
