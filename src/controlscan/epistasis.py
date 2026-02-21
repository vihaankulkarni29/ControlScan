"""
Epistasis Module: Co-mutation matrix and probability calculations.

This module provides tools for analyzing co-occurrence patterns of mutations
across samples and calculating mutual exclusivity/inclusivity metrics.
"""

import logging
from typing import Dict, List, Union, Tuple
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class MutationMatrix:
    """
    Analyzes co-occurrence patterns of mutations across samples.
    
    Attributes:
        samples (List[Dict]): List of sample dictionaries with mutation data.
        co_occurrence_matrix (pd.DataFrame): N x N matrix of co-occurrence counts.
        unique_mutations (List[str]): Sorted list of unique mutations.
    """
    
    def __init__(self, samples: Union[List[Dict[str, List[str]]], pd.DataFrame]):
        """
        Initialize MutationMatrix with sample data.
        
        Args:
            samples: List of dicts with structure {'Sample_ID': str, 'Mutations': List[str]}
                    or a pandas DataFrame with 'Sample_ID' and 'Mutations' columns.
        
        Raises:
            ValueError: If input is empty or malformed.
        """
        logger.info("Initializing MutationMatrix")
        
        if isinstance(samples, pd.DataFrame):
            if samples.empty:
                raise ValueError("Input DataFrame is empty")
            self.samples = samples.to_dict('records')
        elif isinstance(samples, list):
            if not samples:
                raise ValueError("Samples list is empty")
            self.samples = samples
        else:
            raise ValueError("Samples must be a list of dicts or a pandas DataFrame")
        
        self.unique_mutations: List[str] = []
        self.co_occurrence_matrix: pd.DataFrame = pd.DataFrame()
        
        logger.debug(f"Loaded {len(self.samples)} samples")
    
    def build_co_occurrence_matrix(self) -> pd.DataFrame:
        """
        Build an N x N co-occurrence matrix for mutations.
        
        Each cell (i, j) contains the count of samples where Mutation_i 
        and Mutation_j occur together.
        
        Returns:
            pd.DataFrame: Symmetric co-occurrence matrix with mutations as index/columns.
        """
        logger.info("Building co-occurrence matrix")
        
        # Collect all unique mutations
        all_mutations = set()
        for sample in self.samples:
            mutations = sample.get('Mutations', [])
            if isinstance(mutations, str):
                mutations = [mutations]  # Handle single mutation
            all_mutations.update(mutations)
        
        self.unique_mutations = sorted(list(all_mutations))
        n_mutations = len(self.unique_mutations)
        
        logger.debug(f"Found {n_mutations} unique mutations")
        
        # Initialize co-occurrence matrix
        co_matrix = np.zeros((n_mutations, n_mutations), dtype=np.int32)
        mutation_to_idx = {mut: idx for idx, mut in enumerate(self.unique_mutations)}
        
        # Count co-occurrences
        for sample in self.samples:
            mutations = sample.get('Mutations', [])
            if isinstance(mutations, str):
                mutations = [mutations]
            
            mutation_indices = [mutation_to_idx[mut] for mut in mutations if mut in mutation_to_idx]
            
            # Mark co-occurrences (including self-occurrence on diagonal)
            for i in mutation_indices:
                for j in mutation_indices:
                    co_matrix[i, j] += 1
        
        # Convert to DataFrame
        self.co_occurrence_matrix = pd.DataFrame(
            co_matrix,
            index=self.unique_mutations,
            columns=self.unique_mutations
        )
        
        logger.info(f"Co-occurrence matrix built: {self.co_occurrence_matrix.shape}")
        return self.co_occurrence_matrix
    
    def calculate_exclusivity(self, mutation_a: str, mutation_b: str) -> float:
        """
        Calculate mutual exclusivity score using Jaccard Index.
        
        Returns a score from 0 to 1:
        - 1.0: Mutually Inclusive (always occur together)
        - 0.0: Mutually Exclusive (never occur together)
        - 0.5: Random association
        
        The score is calculated as: P(A ∩ B) / P(A ∪ B)
        
        Args:
            mutation_a (str): First mutation identifier.
            mutation_b (str): Second mutation identifier.
        
        Returns:
            float: Exclusivity score in range [0, 1].
        
        Raises:
            ValueError: If mutations not found in matrix.
        """
        logger.debug(f"Calculating exclusivity for {mutation_a} <-> {mutation_b}")
        
        if self.co_occurrence_matrix.empty:
            logger.warning("Co-occurrence matrix not built. Building now...")
            self.build_co_occurrence_matrix()
        
        if mutation_a not in self.unique_mutations or mutation_b not in self.unique_mutations:
            raise ValueError(f"One or both mutations not found in dataset: {mutation_a}, {mutation_b}")
        
        # Get counts from co-occurrence matrix
        co_occur = self.co_occurrence_matrix.loc[mutation_a, mutation_b]
        count_a = self.co_occurrence_matrix.loc[mutation_a, mutation_a]
        count_b = self.co_occurrence_matrix.loc[mutation_b, mutation_b]
        
        # Jaccard Index: |A ∩ B| / |A ∪ B|
        if count_a == 0 or count_b == 0:
            logger.warning(f"Mutation {mutation_a} or {mutation_b} has zero count")
            return 0.0
        
        union = count_a + count_b - co_occur
        jaccard_index = co_occur / union if union > 0 else 0.0
        
        logger.debug(f"Jaccard Index: {jaccard_index:.4f} (co-occur={co_occur}, union={union})")
        return float(jaccard_index)
    
    def calculate_conditional_probability(self, mutation_given: str, mutation_conditional: str) -> float:
        """
        Calculate conditional probability P(A | B).
        
        Returns probability that mutation_given occurs given that mutation_conditional is present.
        
        Args:
            mutation_given (str): Mutation whose probability we're calculating.
            mutation_conditional (str): Condition mutation.
        
        Returns:
            float: Conditional probability in range [0, 1].
        
        Raises:
            ValueError: If mutations not found in matrix.
        """
        logger.debug(f"Calculating P({mutation_given} | {mutation_conditional})")
        
        if self.co_occurrence_matrix.empty:
            self.build_co_occurrence_matrix()
        
        if mutation_given not in self.unique_mutations or mutation_conditional not in self.unique_mutations:
            raise ValueError(f"One or both mutations not found: {mutation_given}, {mutation_conditional}")
        
        co_occur = self.co_occurrence_matrix.loc[mutation_given, mutation_conditional]
        count_conditional = self.co_occurrence_matrix.loc[mutation_conditional, mutation_conditional]
        
        if count_conditional == 0:
            logger.warning(f"Condition mutation {mutation_conditional} has zero count")
            return 0.0
        
        prob = co_occur / count_conditional
        logger.debug(f"P({mutation_given} | {mutation_conditional}) = {prob:.4f}")
        return float(prob)
    
    def get_mutation_summary(self) -> pd.DataFrame:
        """
        Get summary statistics for each mutation.
        
        Returns:
            pd.DataFrame: Summary with mutation counts and frequencies.
        """
        logger.info("Generating mutation summary")
        
        if self.co_occurrence_matrix.empty:
            self.build_co_occurrence_matrix()
        
        diagonal = np.diag(self.co_occurrence_matrix)
        total_samples = len(self.samples)
        
        summary = pd.DataFrame({
            'Mutation': self.unique_mutations,
            'Sample_Count': diagonal,
            'Frequency': diagonal / total_samples
        })
        
        return summary.sort_values('Sample_Count', ascending=False)
