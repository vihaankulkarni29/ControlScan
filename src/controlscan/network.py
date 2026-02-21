"""
Network Module: Co-expression network analysis.

This module constructs and analyzes gene co-expression networks
using correlation metrics and identifies regulatory relationships.
"""

import logging
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


class NetworkBuilder:
    """
    Constructs and analyzes gene co-expression networks.
    
    Provides methods for calculating correlation metrics between genes
    and identifying potential regulatory relationships based on
    expression co-variation.
    
    Attributes:
        expression_matrix (pd.DataFrame): Gene expression matrix (genes x samples).
        correlation_cache (Dict): Cached correlation results.
    """
    
    def __init__(self, expression_matrix: Optional[pd.DataFrame] = None):
        """
        Initialize NetworkBuilder.
        
        Args:
            expression_matrix (Optional[pd.DataFrame]): Expression matrix with genes as rows
                                                       and samples as columns.
        """
        logger.info("Initializing NetworkBuilder")
        self.expression_matrix: Optional[pd.DataFrame] = expression_matrix.copy() if expression_matrix is not None else None
        self.correlation_cache: Dict = {}
        self.pvalue_cache: Dict = {}
        
        if expression_matrix is not None:
            logger.debug(f"Loaded matrix: {expression_matrix.shape}")
    
    def calculate_correlation(
        self,
        target_gene: str,
        expression_matrix: Optional[pd.DataFrame] = None,
        method: str = 'pearson',
        min_samples: int = 3
    ) -> pd.DataFrame:
        """
        Calculate expression correlation between a target gene and all other genes.
        
        Uses Pearson, Spearman, or Kendall correlation metrics to identify genes
        that co-vary with the target gene across samples.
        
        Args:
            target_gene (str): Gene identifier to calculate correlations for.
            expression_matrix (Optional[pd.DataFrame]): Expression matrix. If None, uses initialized matrix.
            method (str): Correlation method ('pearson', 'spearman', 'kendall'). Default: 'pearson'.
            min_samples (int): Minimum number of samples required (raises error if fewer). Default: 3.
        
        Returns:
            pd.DataFrame: DataFrame with columns:
                - Gene: Gene identifier
                - Correlation: Correlation coefficient (-1 to 1)
                - P_Value: Two-tailed p-value
                - Significant: Boolean (p < 0.05)
        
        Raises:
            ValueError: If target gene not found, insufficient samples, or invalid method.
        
        Examples:
            >>> builder = NetworkBuilder(expression_matrix)
            >>> correlations = builder.calculate_correlation('acrB', method='pearson')
            >>> top_correlated = correlations.nlargest(10, 'Correlation')
        """
        logger.info(f"Calculating {method} correlations for target gene: {target_gene}")
        
        # Get matrix to use
        if expression_matrix is not None:
            matrix = expression_matrix.copy()
            self.expression_matrix = matrix.copy()
        elif self.expression_matrix is not None:
            matrix = self.expression_matrix.copy()
        else:
            raise ValueError("No expression matrix provided or initialized")
        
        # Validate inputs
        if matrix.empty:
            raise ValueError("Expression matrix is empty")
        
        if target_gene not in matrix.index:
            raise ValueError(f"Target gene '{target_gene}' not found in matrix. Available: {list(matrix.index)[:5]}...")
        
        if matrix.shape[1] < min_samples:
            raise ValueError(f"Insufficient samples ({matrix.shape[1]} < {min_samples})")
        
        if method not in ['pearson', 'spearman', 'kendall']:
            raise ValueError(f"Invalid method '{method}'. Must be 'pearson', 'spearman', or 'kendall'")
        
        # Get target gene expression
        target_expr = matrix.loc[target_gene].values.astype(float)
        
        # Check for validity
        if np.isnan(target_expr).all():
            raise ValueError(f"Target gene '{target_gene}' has all NaN values")
        
        # Calculate correlations
        correlations = []
        pvalues = []
        
        for gene in matrix.index:
            gene_expr = matrix.loc[gene].values.astype(float)
            
            # Skip if all NaN
            if np.isnan(gene_expr).all() or np.isnan(target_expr).all():
                correlations.append(np.nan)
                pvalues.append(np.nan)
                continue
            
            # Calculate correlation and p-value
            if method == 'pearson':
                corr, pval = stats.pearsonr(target_expr, gene_expr)
            elif method == 'spearman':
                corr, pval = stats.spearmanr(target_expr, gene_expr)
            elif method == 'kendall':
                corr, pval = stats.kendalltau(target_expr, gene_expr)
            
            correlations.append(corr)
            pvalues.append(pval)
        
        # Create results DataFrame
        results = pd.DataFrame({
            'Gene': matrix.index,
            'Correlation': correlations,
            'P_Value': pvalues,
        })
        
        # Add significance column (p < 0.05)
        results['Significant'] = results['P_Value'] < 0.05
        
        # Sort by absolute correlation
        results['Abs_Correlation'] = results['Correlation'].abs()
        results = results.sort_values('Abs_Correlation', ascending=False).drop('Abs_Correlation', axis=1)
        
        # Cache results
        cache_key = f"{target_gene}_{method}"
        self.correlation_cache[cache_key] = results.copy()
        
        logger.info(
            f"Calculated correlations for {target_gene}: "
            f"{(results['Significant']).sum()} significant at p<0.05"
        )
        
        return results
    
    def build_correlation_network(
        self,
        expression_matrix: Optional[pd.DataFrame] = None,
        correlation_threshold: float = 0.7,
        pvalue_threshold: float = 0.05,
        method: str = 'pearson'
    ) -> pd.DataFrame:
        """
        Build a complete gene co-expression network (all pairwise correlations).
        
        Args:
            expression_matrix (Optional[pd.DataFrame]): Expression matrix.
            correlation_threshold (float): Minimum absolute correlation to include edge.
            pvalue_threshold (float): Maximum p-value to include edge.
            method (str): Correlation method ('pearson', 'spearman', 'kendall').
        
        Returns:
            pd.DataFrame: Edge list with columns:
                - Source: First gene
                - Target: Second gene
                - Correlation: Correlation coefficient
                - P_Value: Statistical significance
        """
        logger.info(f"Building {method} co-expression network (|r| >= {correlation_threshold})")
        
        if expression_matrix is not None:
            matrix = expression_matrix.copy()
            self.expression_matrix = matrix.copy()
        elif self.expression_matrix is not None:
            matrix = self.expression_matrix.copy()
        else:
            raise ValueError("No expression matrix provided or initialized")
        
        # Calculate correlation matrix
        if method == 'pearson':
            corr_matrix, pval_matrix = self._correlation_matrix_pearson(matrix)
        elif method == 'spearman':
            corr_matrix, pval_matrix = self._correlation_matrix_spearman(matrix)
        elif method == 'kendall':
            corr_matrix, pval_matrix = self._correlation_matrix_kendall(matrix)
        else:
            raise ValueError(f"Invalid method: {method}")
        
        # Extract edges above threshold
        edges = []
        for i in range(len(matrix.index)):
            for j in range(i + 1, len(matrix.index)):
                corr = corr_matrix[i, j]
                pval = pval_matrix[i, j]
                
                if abs(corr) >= correlation_threshold and pval < pvalue_threshold:
                    edges.append({
                        'Source': matrix.index[i],
                        'Target': matrix.index[j],
                        'Correlation': corr,
                        'P_Value': pval
                    })
        
        edge_df = pd.DataFrame(edges)
        logger.info(f"Network contains {len(edge_df)} edges")
        
        return edge_df.sort_values('Correlation', key=abs, ascending=False)
    
    @staticmethod
    def _correlation_matrix_pearson(matrix: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate Pearson correlation matrix."""
        corr_matrix = np.corrcoef(matrix.values)
        
        # Calculate p-values
        n = matrix.shape[1]
        pval_matrix = np.zeros_like(corr_matrix)
        
        for i in range(len(matrix)):
            for j in range(i + 1, len(matrix)):
                _, pval = stats.pearsonr(matrix.iloc[i].values, matrix.iloc[j].values)
                pval_matrix[i, j] = pval
                pval_matrix[j, i] = pval
        
        return corr_matrix, pval_matrix
    
    @staticmethod
    def _correlation_matrix_spearman(matrix: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate Spearman correlation matrix."""
        corr_matrix = np.zeros((len(matrix), len(matrix)))
        pval_matrix = np.zeros((len(matrix), len(matrix)))
        
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                if i == j:
                    corr_matrix[i, j] = 1.0
                    pval_matrix[i, j] = 0.0
                else:
                    corr, pval = stats.spearmanr(matrix.iloc[i].values, matrix.iloc[j].values)
                    corr_matrix[i, j] = corr
                    pval_matrix[i, j] = pval
        
        return corr_matrix, pval_matrix
    
    @staticmethod
    def _correlation_matrix_kendall(matrix: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate Kendall correlation matrix."""
        corr_matrix = np.zeros((len(matrix), len(matrix)))
        pval_matrix = np.zeros((len(matrix), len(matrix)))
        
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                if i == j:
                    corr_matrix[i, j] = 1.0
                    pval_matrix[i, j] = 0.0
                else:
                    corr, pval = stats.kendalltau(matrix.iloc[i].values, matrix.iloc[j].values)
                    corr_matrix[i, j] = corr
                    pval_matrix[i, j] = pval
        
        return corr_matrix, pval_matrix
    
    def get_regulatory_candidates(
        self,
        target_gene: str,
        top_n: int = 10,
        expression_matrix: Optional[pd.DataFrame] = None
    ) -> pd.DataFrame:
        """
        Get top regulatory candidates for a target gene.
        
        Returns genes with strongest correlation to target, useful for
        identifying potential regulators or regulatory targets.
        
        Args:
            target_gene (str): Target gene identifier.
            top_n (int): Number of top candidates to return.
            expression_matrix (Optional[pd.DataFrame]): Expression matrix.
        
        Returns:
            pd.DataFrame: Top correlated genes with correlation and p-values.
        """
        logger.info(f"Identifying top {top_n} regulatory candidates for {target_gene}")
        
        results = self.calculate_correlation(target_gene, expression_matrix)
        
        # Remove self-correlation
        results = results[results['Gene'] != target_gene]
        
        # Get top by absolute correlation
        results['Abs_Correlation'] = results['Correlation'].abs()
        top_genes = results.nlargest(top_n, 'Abs_Correlation')
        
        return top_genes.drop('Abs_Correlation', axis=1)
