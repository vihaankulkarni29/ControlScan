"""
Fetcher Module: Data retrieval from GEO (Gene Expression Omnibus).

This module provides utilities for downloading and parsing gene expression data
from the NCBI GEO database using the GEOparse library.
"""

import logging
from typing import Optional, Tuple
import pandas as pd

logger = logging.getLogger(__name__)


class GEOFetcher:
    """
    Fetches gene expression data from NCBI Gene Expression Omnibus (GEO).
    
    Attributes:
        geo_id (str): GEO accession ID (e.g., 'GSE12345').
        expression_matrix (pd.DataFrame): Downloaded expression data.
    """
    
    def __init__(self):
        """Initialize GEOFetcher."""
        logger.info("Initializing GEOFetcher (GEOparse required)")
        self.geo_id: Optional[str] = None
        self.expression_matrix: pd.DataFrame = pd.DataFrame()
        self._geoparse_available = self._check_geoparse()
    
    @staticmethod
    def _check_geoparse() -> bool:
        """
        Check if GEOparse library is available.
        
        Returns:
            bool: True if GEOparse can be imported, False otherwise.
        """
        try:
            import GEOparse
            logger.debug("GEOparse library available")
            return True
        except ImportError:
            logger.warning("GEOparse not installed. Install with: pip install GEOparse")
            return False
    
    def download_matrix(
        self,
        geo_id: str,
        gsm_subset: Optional[list] = None,
        gene_subset: Optional[list] = None
    ) -> pd.DataFrame:
        """
        Download gene expression matrix from GEO.
        
        This method fetches SOFT files and extracts expression count matrices
        from a specified GEO Series (GSE) or GEO Sample (GSM) identifier.
        
        Args:
            geo_id (str): GEO accession ID (e.g., 'GSE12345', 'GSM987654').
            gsm_subset (Optional[list]): Subset of GSM IDs to include. If None, all samples included.
            gene_subset (Optional[list]): Subset of genes to include. If None, all genes included.
        
        Returns:
            pd.DataFrame: Expression matrix with genes as rows, samples as columns.
                         Contains normalized count values.
        
        Raises:
            ImportError: If GEOparse is not installed.
            ValueError: If GEO ID is invalid or data cannot be retrieved.
        
        Examples:
            >>> fetcher = GEOFetcher()
            >>> matrix = fetcher.download_matrix('GSE12345')
            >>> print(matrix.shape)  # (genes, samples)
        """
        logger.info(f"Downloading expression matrix for {geo_id}")
        
        if not self._geoparse_available:
            raise ImportError(
                "GEOparse is required. Install with: pip install GEOparse"
            )
        
        try:
            import GEOparse
        except ImportError:
            raise ImportError("GEOparse library not available")
        
        try:
            # Fetch GEO object
            logger.debug(f"Fetching GEO object: {geo_id}")
            geo = GEOparse.get_GEO(geo_id, silent=False)
            
            # Extract sample data
            expression_data = {}
            
            if hasattr(geo, 'phenotype_data'):
                # For GSE (Series)
                logger.debug("Parsing GSE (Series) object")
                samples = geo.gsms
                
                for gsm_id, gsm in samples.items():
                    # Filter by GSM subset if provided
                    if gsm_subset and gsm_id not in gsm_subset:
                        continue
                    
                    # Extract expression values
                    if hasattr(gsm, 'table') and not gsm.table.empty:
                        table = gsm.table
                        
                        # Attempt to find expression column
                        expr_col = None
                        for col in ['VALUE', 'SIGNAL', 'counts', 'expression']:
                            if col in table.columns:
                                expr_col = col
                                break
                        
                        if expr_col is None and len(table.columns) > 0:
                            expr_col = table.columns[-1]
                        
                        if expr_col:
                            expression_data[gsm_id] = pd.Series(
                                table[expr_col].values,
                                index=table['ID_REF'].values if 'ID_REF' in table.columns else range(len(table))
                            )
            else:
                # For individual GSM (Sample)
                logger.debug("Parsing GSM (Sample) object")
                if hasattr(geo, 'table') and not geo.table.empty:
                    table = geo.table
                    expr_col = None
                    for col in ['VALUE', 'SIGNAL', 'counts', 'expression']:
                        if col in table.columns:
                            expr_col = col
                            break
                    
                    if expr_col is None and len(table.columns) > 0:
                        expr_col = table.columns[-1]
                    
                    if expr_col:
                        expression_data[geo_id] = pd.Series(
                            table[expr_col].values,
                            index=table['ID_REF'].values if 'ID_REF' in table.columns else range(len(table))
                        )
            
            # Combine into DataFrame
            if expression_data:
                self.expression_matrix = pd.DataFrame(expression_data)
                
                # Apply gene subset filter
                if gene_subset:
                    self.expression_matrix = self.expression_matrix.loc[
                        self.expression_matrix.index.isin(gene_subset)
                    ]
                
                logger.info(
                    f"Successfully retrieved matrix: {self.expression_matrix.shape[0]} genes, "
                    f"{self.expression_matrix.shape[1]} samples"
                )
                self.geo_id = geo_id
                return self.expression_matrix
            else:
                raise ValueError(f"No expression data found for {geo_id}")
        
        except Exception as e:
            logger.error(f"Failed to download {geo_id}: {str(e)}")
            raise ValueError(f"Cannot retrieve expression matrix for {geo_id}: {str(e)}")
    
    def get_sample_metadata(self) -> Optional[pd.DataFrame]:
        """
        Get metadata for samples in the downloaded dataset.
        
        Returns:
            pd.DataFrame: Sample metadata with characteristics and descriptions.
        """
        logger.info(f"Extracting metadata for {self.geo_id}")
        
        if self.expression_matrix.empty:
            logger.warning("No expression matrix loaded")
            return None
        
        try:
            import GEOparse
            geo = GEOparse.get_GEO(self.geo_id, silent=True)
            
            metadata = []
            for gsm_id in self.expression_matrix.columns:
                if gsm_id in geo.gsms:
                    gsm = geo.gsms[gsm_id]
                    metadata.append({
                        'Sample_ID': gsm_id,
                        'Title': gsm.metadata.get('title', [''])[0] if hasattr(gsm, 'metadata') else '',
                        'Organism': gsm.metadata.get('organism_ch1', [''])[0] if hasattr(gsm, 'metadata') else '',
                    })
            
            return pd.DataFrame(metadata) if metadata else None
        
        except Exception as e:
            logger.warning(f"Could not retrieve metadata: {str(e)}")
            return None
