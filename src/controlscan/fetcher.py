"""
Fetcher Module: Data retrieval from GEO (Gene Expression Omnibus).

This module provides utilities for downloading and parsing gene expression data
from the NCBI GEO database using the GEOparse library.

The GEOFetcher class automates SOFT file download, expression matrix extraction,
and sample metadata retrieval for transcriptomic analysis workflows.
"""

import logging
import os
from typing import Optional
import pandas as pd

logger = logging.getLogger(__name__)


class GEOFetcher:
    """
    Fetches gene expression data from NCBI Gene Expression Omnibus (GEO).
    
    Automatically downloads SOFT files from GEO and extracts expression matrices
    and phenotype data for downstream analysis.
    
    Attributes:
        dest_dir (str): Local directory where downloaded files are stored.
        geo_id (str): Most recently accessed GEO accession ID.
        expression_matrix (pd.DataFrame): Last downloaded expression matrix.
        phenotype_data (pd.DataFrame): Last downloaded sample metadata.
    
    Example:
        >>> fetcher = GEOFetcher(dest_dir="./data/raw")
        >>> matrix = fetcher.download_matrix("GSE45243")
        >>> print(matrix.shape)  # (genes, samples)
        >>> metadata = fetcher.get_phenotype_data("GSE45243")
        >>> print(metadata.head())
    """
    
    def __init__(self, dest_dir: str = "./data/raw"):
        """
        Initialize GEOFetcher with download destination.
        
        Args:
            dest_dir (str): Directory path for storing downloaded files.
                           Created if it doesn't exist. Default: "./data/raw"
        
        Raises:
            OSError: If directory cannot be created.
        """
        logger.info(f"Initializing GEOFetcher with dest_dir={dest_dir}")
        
        self.dest_dir = dest_dir
        self.geo_id: Optional[str] = None
        self.expression_matrix: pd.DataFrame = pd.DataFrame()
        self.phenotype_data: pd.DataFrame = pd.DataFrame()
        
        # Create destination directory if it doesn't exist
        try:
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir, exist_ok=True)
                logger.info(f"Created directory: {dest_dir}")
            else:
                logger.debug(f"Directory already exists: {dest_dir}")
        except OSError as e:
            logger.error(f"Failed to create directory {dest_dir}: {e}")
            raise
        
        self._check_geoparse()
    
    @staticmethod
    def _check_geoparse() -> bool:
        """
        Check if GEOparse library is available.
        
        Returns:
            bool: True if GEOparse is installed, False otherwise.
        """
        try:
            import GEOparse
            logger.debug("GEOparse library is available")
            return True
        except ImportError:
            logger.warning("GEOparse not installed. Install with: pip install GEOparse")
            return False
    
    def download_matrix(self, geo_id: str) -> pd.DataFrame:
        """
        Download expression matrix from GEO.
        
        Fetches the SOFT file for a given GEO Series (GSE) ID and extracts
        the normalized expression matrix. Rows = genes/probes, Columns = samples.
        
        Args:
            geo_id (str): GEO accession ID (e.g., 'GSE45243', 'GSE12345').
        
        Returns:
            pd.DataFrame: Expression matrix with genes/probes as index and 
                         sample IDs as columns. Values are normalized expression
                         (typically log2 or TPM).
        
        Raises:
            ImportError: If GEOparse is not installed.
            ValueError: If download fails or matrix cannot be extracted.
            Exception: Any error raised by GEOparse during download.
        
        Examples:
            >>> fetcher = GEOFetcher()
            >>> matrix = fetcher.download_matrix("GSE45243")
            >>> print(f"Downloaded {matrix.shape[0]} genes x {matrix.shape[1]} samples")
        """
        logger.info(f"Downloading expression matrix for {geo_id}")
        
        try:
            import GEOparse
        except ImportError:
            logger.error("GEOparse is required but not installed")
            raise ImportError("Install GEOparse with: pip install GEOparse")
        
        try:
            # Download the GEO object
            logger.debug(f"Fetching GEO object: {geo_id} to {self.dest_dir}")
            gse = GEOparse.get_GEO(geo=geo_id, destdir=self.dest_dir, silent=False)
            
            logger.debug(f"Successfully retrieved GEO object. Type: {type(gse)}")
            
            # Try to extract expression matrix by pivoting samples
            logger.debug("Attempting to extract expression matrix...")
            try:
                matrix = gse.pivot_samples('VALUE')
            except (KeyError, ValueError) as e:
                logger.warning(f"pivot_samples('VALUE') failed ({e}). Trying alternative approach...")
                
                # Alternative approach: manually construct the matrix from samples
                expression_data = {}
                
                if hasattr(gse, 'gsms') and gse.gsms:
                    logger.debug(f"Found {len(gse.gsms)} samples. Constructing matrix manually...")
                    
                    for gsm_id, gsm_obj in gse.gsms.items():
                        if hasattr(gsm_obj, 'table') and not gsm_obj.table.empty:
                            table = gsm_obj.table
                            
                            # Find the VALUE column (or similar expression column)
                            expr_col = None
                            for col in ['VALUE', 'SIGNAL', 'Cy3', 'Cy5', 'Signal']:
                                if col in table.columns:
                                    expr_col = col
                                    break
                            
                            if expr_col is None and len(table.columns) > 1:
                                # Default to last column if VALUE not found
                                expr_col = table.columns[-1]
                            
                            if expr_col:
                                # Use ID_REF, ID, or Gene ID as index
                                index_col = None
                                for idx in ['ID_REF', 'ID', 'GeneID', 'ProbeID']:
                                    if idx in table.columns:
                                        index_col = idx
                                        break
                                
                                if index_col:
                                    expression_data[gsm_id] = pd.Series(
                                        table[expr_col].values,
                                        index=table[index_col].values
                                    )
                                else:
                                    # Use table index if no ID column found
                                    expression_data[gsm_id] = pd.Series(
                                        table[expr_col].values,
                                        index=range(len(table))
                                    )
                    
                    if expression_data:
                        matrix = pd.DataFrame(expression_data)
                    else:
                        raise ValueError(f"Could not extract expression data from any samples in {geo_id}")
                else:
                    raise ValueError(f"No sample data (gsms) found in {geo_id}")
            
            # Validate the extracted matrix
            if matrix is None or matrix.empty:
                raise ValueError(f"No expression matrix extracted from {geo_id}")
            
            logger.info(
                f"Successfully extracted matrix: {matrix.shape[0]} genes/probes x "
                f"{matrix.shape[1]} samples"
            )
            
            # Store for later reference
            self.geo_id = geo_id
            self.expression_matrix = matrix
            
            return matrix
        
        except ImportError as e:
            logger.error(f"Import error: {e}")
            raise
        except Exception as e:
            logger.error(f"Failed to download/extract matrix for {geo_id}: {e}")
            raise ValueError(f"Cannot retrieve expression matrix for {geo_id}: {str(e)}")
    
    def get_phenotype_data(self, geo_id: str) -> pd.DataFrame:
        """
        Fetch sample phenotype/metadata from GEO.
        
        Retrieves the sample characteristics, treatment information, and other
        phenotypic annotations stored in the GEO database.
        
        Args:
            geo_id (str): GEO accession ID (e.g., 'GSE45243').
        
        Returns:
            pd.DataFrame: Phenotype data with samples as rows and characteristics
                         as columns. Typically includes:
                         - Sample ID
                         - Treatment/Condition
                         - Strain/Organism
                         - Platform/Protocol info
                         - Other experimental metadata
        
        Raises:
            ImportError: If GEOparse is not installed.
            ValueError: If phenotype data cannot be retrieved.
        
        Examples:
            >>> fetcher = GEOFetcher()
            >>> metadata = fetcher.get_phenotype_data("GSE45243")
            >>> print(metadata.columns)
            >>> print(metadata.head())
        """
        logger.info(f"Fetching phenotype data for {geo_id}")
        
        try:
            import GEOparse
        except ImportError:
            logger.error("GEOparse is required but not installed")
            raise ImportError("Install GEOparse with: pip install GEOparse")
        
        try:
            # Fetch the GEO object
            logger.debug(f"Retrieving GEO object for metadata: {geo_id}")
            gse = GEOparse.get_GEO(geo=geo_id, destdir=self.dest_dir, silent=False)
            
            # Extract phenotype data
            phenotype_df = gse.phenotype_data
            
            if phenotype_df is None or phenotype_df.empty:
                logger.warning(f"No phenotype data found for {geo_id}")
                return pd.DataFrame()
            
            logger.info(
                f"Retrieved phenotype data: {phenotype_df.shape[0]} samples x "
                f"{phenotype_df.shape[1]} characteristics"
            )
            
            # Store for later reference
            self.phenotype_data = phenotype_df
            
            return phenotype_df
        
        except Exception as e:
            logger.error(f"Failed to retrieve phenotype data for {geo_id}: {e}")
            raise ValueError(f"Cannot retrieve phenotype data for {geo_id}: {str(e)}")
    
    def Youtube(self, geo_id: str) -> pd.DataFrame:
        """
        Alias method for get_phenotype_data (convenience wrapper).
        
        Fetches the GEO object and returns phenotype_data for sample metadata.
        This allows matching samples to their treatments, strains, and conditions.
        
        Args:
            geo_id (str): GEO accession ID.
        
        Returns:
            pd.DataFrame: Phenotype/metadata DataFrame for all samples in the GEO series.
        
        Note:
            This is a convenience alias. Use get_phenotype_data() for clarity.
        """
        logger.debug(f"Youtube() called as alias for get_phenotype_data({geo_id})")
        return self.get_phenotype_data(geo_id)


# ============================================================================
# TESTING BLOCK
# ============================================================================

if __name__ == "__main__":
    """
    Test the GEOFetcher module by downloading a real E. coli dataset from GEO.
    
    Dataset: GSE45243 - E. coli RNA-seq expression data
    Description: Small, publicly available dataset suitable for testing
    """
    
    # Configure logging for test output
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    print("\n" + "="*70)
    print("GEOFetcher Module Test")
    print("="*70)
    
    try:
        # Initialize fetcher
        print("\n[1] Initializing GEOFetcher...")
        fetcher = GEOFetcher(dest_dir="./data/raw")
        print("[PASS] GEOFetcher initialized successfully")
        
        # Download expression matrix
        # Using GSE123 (a classic microarray dataset with 12,488 probes and 3 samples)
        test_geo_id = "GSE123"
        print(f"\n[2] Downloading expression matrix for {test_geo_id}...")
        print("    (This may take a minute on first run)")
        matrix = fetcher.download_matrix(test_geo_id)
        
        print(f"[PASS] Successfully downloaded expression matrix")
        print(f"  Shape: {matrix.shape[0]} genes x {matrix.shape[1]} samples")
        
        # Display matrix head
        print("\n[3] Expression Matrix (head):")
        print("-" * 70)
        print(matrix.head())
        
        # Display matrix info
        print("\n[4] Expression Matrix Info:")
        print("-" * 70)
        print(f"Data type: {matrix.dtypes.unique()}")
        print(f"Index (genes/probes): {matrix.index[:5].tolist()}...")
        print(f"Columns (samples): {matrix.columns.tolist()}")
        print(f"Value range: [{matrix.values.min():.2f}, {matrix.values.max():.2f}]")
        
        # Fetch phenotype data
        print(f"\n[5] Fetching sample phenotype data for {test_geo_id}...")
        phenotype = fetcher.get_phenotype_data(test_geo_id)
        
        if not phenotype.empty:
            print(f"[PASS] Retrieved phenotype data")
            print(f"  Shape: {phenotype.shape[0]} samples x {phenotype.shape[1]} characteristics")
            print("\n[6] Phenotype Data (head):")
            print("-" * 70)
            print(phenotype.head())
        else:
            print("[WARN] No phenotype data available for this dataset")
        
        # Test Youtube alias method
        print(f"\n[7] Testing Youtube() alias method for {test_geo_id}...")
        phenotype_alias = fetcher.Youtube(test_geo_id)
        print(f"[PASS] Youtube() alias works correctly (returned {phenotype_alias.shape[0]} samples)")
        
        print("\n" + "="*70)
        print("[PASS] All tests completed successfully!")
        print("="*70 + "\n")
        
    except Exception as e:
        print("\n" + "="*70)
        print(f"[FAIL] Test failed with error:")
        print(f"  {type(e).__name__}: {e}")
        print("="*70 + "\n")
        raise
