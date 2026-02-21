"""
Normalizer Module: RNA-seq expression normalization.

This module provides methods for normalizing gene expression data
into standard metrics like TPM (Transcripts Per Million).
"""

import logging
from typing import Optional
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class ExpressionNormalizer:
    """
    Normalizes RNA-seq count matrices to standard expression metrics.
    
    Supports multiple normalization strategies:
        - TPM (Transcripts Per Million): Accounts for gene length and sequencing depth
        - Log-transformed TPM: Stabilizes variance
        - Z-score normalization: Centers and scales by variance
    
    Attributes:
        raw_matrix (pd.DataFrame): Original count matrix.
        normalized_matrix (pd.DataFrame): Normalized expression values.
        normalization_type (str): Type of normalization applied.
    """
    
    def __init__(self, expression_matrix: Optional[pd.DataFrame] = None):
        """
        Initialize ExpressionNormalizer.
        
        Args:
            expression_matrix (Optional[pd.DataFrame]): Raw count matrix to normalize.
                                                       Genes as rows, samples as columns.
        """
        logger.info("Initializing ExpressionNormalizer")
        self.raw_matrix: Optional[pd.DataFrame] = expression_matrix.copy() if expression_matrix is not None else None
        self.normalized_matrix: pd.DataFrame = pd.DataFrame()
        self.normalization_type: str = ""
        
        if expression_matrix is not None:
            logger.debug(f"Loaded matrix: {expression_matrix.shape}")
    
    def normalize_tpm(
        self,
        expression_matrix: Optional[pd.DataFrame] = None,
        gene_lengths: Optional[pd.Series] = None,
        log_transform: bool = False,
        pseudocount: float = 1.0
    ) -> pd.DataFrame:
        """
        Normalize counts to Transcripts Per Million (TPM).
        
        TPM normalization accounts for both:
        1. Gene length bias (longer genes accumulate more reads)
        2. Sequencing depth bias (samples with different total reads)
        
        Formula:
            TPM = (read_count / gene_length) / sum(read_count / gene_length) * 1e6
        
        Args:
            expression_matrix (Optional[pd.DataFrame]): Raw count matrix. If None, uses initialized matrix.
                                                       Genes as rows, samples as columns.
            gene_lengths (Optional[pd.Series]): Gene lengths indexed by gene name.
                                              If None, assumes all genes have equal length (length normalization only).
            log_transform (bool): If True, returns log2(TPM + pseudocount). Default: False.
            pseudocount (float): Small value to avoid log(0). Default: 1.0.
        
        Returns:
            pd.DataFrame: TPM-normalized expression matrix with same shape as input.
        
        Raises:
            ValueError: If no matrix is available or genes/samples dimension mismatch.
        
        Examples:
            >>> normalizer = ExpressionNormalizer(raw_counts)
            >>> tpm_matrix = normalizer.normalize_tpm(log_transform=False)
            >>> log_tpm_matrix = normalizer.normalize_tpm(log_transform=True)
        """
        logger.info("Normalizing expression to TPM")
        
        # Get matrix to normalize
        if expression_matrix is not None:
            matrix = expression_matrix.copy()
            self.raw_matrix = matrix.copy()
        elif self.raw_matrix is not None:
            matrix = self.raw_matrix.copy()
        else:
            raise ValueError("No expression matrix provided or initialized")
        
        # Validate data
        if matrix.empty:
            raise ValueError("Expression matrix is empty")
        
        if matrix.isnull().any().any():
            logger.warning("Expression matrix contains NaN values. Replacing with 0.")
            matrix.fillna(0, inplace=True)
        
        if (matrix < 0).any().any():
            logger.warning("Expression matrix contains negative values. Setting to 0.")
            matrix[matrix < 0] = 0
        
        # Apply length normalization if gene lengths provided
        if gene_lengths is not None:
            logger.debug(f"Applying gene length normalization for {len(gene_lengths)} genes")
            
            # Align gene lengths with matrix
            aligned_lengths = gene_lengths.reindex(matrix.index)
            
            if aligned_lengths.isnull().any():
                logger.warning(
                    f"{aligned_lengths.isnull().sum()} genes missing length information. Using 1.0 for missing values."
                )
                aligned_lengths.fillna(1.0, inplace=True)
            
            # Divide by gene length (avoiding division by zero)
            aligned_lengths = aligned_lengths.replace(0, 1)
            length_normalized = matrix.div(aligned_lengths, axis=0)
        else:
            logger.debug("No gene length information provided. Skipping length normalization.")
            length_normalized = matrix
        
        # Calculate per-sample scaling factors (sum of normalized counts)
        scaling_factors = length_normalized.sum(axis=0)
        
        if (scaling_factors == 0).any():
            logger.warning("Some samples have zero total expression. Handling gracefully...")
            scaling_factors = scaling_factors.replace(0, 1)
        
        # Calculate TPM: (normalized / scaling_factor) * 1e6
        tpm_matrix = (length_normalized / scaling_factors) * 1e6
        
        # Apply log transformation if requested
        if log_transform:
            logger.debug(f"Applying log2 transformation with pseudocount={pseudocount}")
            tpm_matrix = np.log2(tpm_matrix + pseudocount)
            self.normalization_type = "log2(TPM)"
        else:
            self.normalization_type = "TPM"
        
        self.normalized_matrix = tpm_matrix
        logger.info(f"TPM normalization complete. Range: [{tpm_matrix.min().min():.2f}, {tpm_matrix.max().max():.2f}]")
        
        return tpm_matrix
    
    def normalize_z_score(self) -> pd.DataFrame:
        """
        Z-score normalize each gene across samples.
        
        Produces a matrix where each gene has:
        - Mean = 0
        - Standard Deviation = 1
        
        Formula: z = (x - mean) / std
        
        Returns:
            pd.DataFrame: Z-score normalized matrix.
        
        Raises:
            ValueError: If no normalized matrix exists.
        """
        logger.info("Applying Z-score normalization")
        
        if self.normalized_matrix.empty:
            raise ValueError("No normalized matrix available. Run normalize_tpm() first.")
        
        # Calculate mean and std for each gene (row)
        gene_means = self.normalized_matrix.mean(axis=1)
        gene_stds = self.normalized_matrix.std(axis=1)
        
        # Avoid division by zero
        gene_stds = gene_stds.replace(0, 1)
        
        # Z-score normalization
        z_matrix = (self.normalized_matrix - gene_means.values[:, np.newaxis]) / gene_stds.values[:, np.newaxis]
        
        self.normalization_type = "Z-score"
        self.normalized_matrix = z_matrix
        
        logger.info(f"Z-score normalization complete. Mean ≈ 0, Std ≈ 1")
        
        return z_matrix
    
    def normalize_quantile(self) -> pd.DataFrame:
        """
        Quantile normalization across samples.
        
        This method ensures all samples have the same distribution of expression values.
        Useful for removing technical batch effects.
        
        Returns:
            pd.DataFrame: Quantile-normalized matrix.
        
        Raises:
            ValueError: If no normalized matrix exists.
        """
        logger.info("Applying Quantile normalization")
        
        if self.normalized_matrix.empty:
            raise ValueError("No normalized matrix available. Run normalize_tpm() first.")
        
        matrix = self.normalized_matrix.copy()
        
        # Sort each sample
        sorted_matrix = matrix.apply(np.sort, axis=0)
        
        # Calculate mean for each rank position
        mean_sorted = sorted_matrix.mean(axis=1)
        
        # Create mapping from sorted values back to original
        quantile_matrix = pd.DataFrame(
            data=np.zeros_like(matrix.values),
            index=matrix.index,
            columns=matrix.columns
        )
        
        for col in matrix.columns:
            # Get ranks of original values
            ranks = matrix[col].argsort().argsort()
            # Map to mean sorted values
            quantile_matrix[col] = mean_sorted.values[ranks]
        
        self.normalization_type = "Quantile"
        self.normalized_matrix = quantile_matrix
        
        logger.info("Quantile normalization complete")
        
        return quantile_matrix
    
    def get_normalization_stats(self) -> pd.DataFrame:
        """
        Get descriptive statistics of the normalized matrix.
        
        Returns:
            pd.DataFrame: Summary statistics (mean, std, min, max, median) per sample.
        """
        logger.info(f"Computing statistics for {self.normalization_type} normalized data")
        
        if self.normalized_matrix.empty:
            raise ValueError("No normalized matrix available")
        
        stats = pd.DataFrame({
            'Sample': self.normalized_matrix.columns,
            'Mean': self.normalized_matrix.mean(axis=0).values,
            'Std': self.normalized_matrix.std(axis=0).values,
            'Min': self.normalized_matrix.min(axis=0).values,
            'Max': self.normalized_matrix.max(axis=0).values,
            'Median': self.normalized_matrix.median(axis=0).values,
        })
        
        return stats.set_index('Sample')
