"""
ControlScan: Standalone Python microservice for Systems Biology.
Process transcriptomic networks and calculate genomic co-mutation probabilities (epistasis).
"""

__version__ = "0.1.0"
__author__ = "Senior Computational Biologist"
__description__ = "Epistasis and Co-expression Analysis Tool"

from .epistasis import MutationMatrix
from .fetcher import GEOFetcher
from .normalizer import ExpressionNormalizer
from .network import NetworkBuilder

__all__ = [
    "MutationMatrix",
    "GEOFetcher",
    "ExpressionNormalizer",
    "NetworkBuilder",
]
