"""
Configuration settings for A112509 research project.
"""

import os

# Project paths
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'output')
NOTEBOOKS_DIR = os.path.join(PROJECT_ROOT, 'notebooks')

# Computation settings
DEFAULT_BEAM_WIDTH = 1000
MAX_BRUTE_FORCE_N = 25  # Beyond this, use beam search
DEFAULT_TIMEOUT_SECONDS = 3600  # 1 hour timeout for long computations

# Visualization settings
PLOT_DPI = 150
PLOT_STYLE = 'seaborn-v0_8-darkgrid'
FIGURE_SIZE = (12, 8)

# Analysis settings
MIN_VALUES_FOR_ANALYSIS = 10  # Minimum sequence values needed for statistical analysis
DIFFERENCE_TABLE_MAX_ORDER = 10

# Beam search settings
BEAM_SEARCH_CONFIG = {
    'default_width': 1000,
    'adaptive_initial_width': 500,
    'adaptive_max_width': 5000,
    'adaptive_expansion_factor': 1.5,
    'multi_start_count': 5
}

# Bounds analysis
BOUNDS_CONFIG = {
    'use_shallit_bound': True,
    'use_information_theoretic': True,
    'use_debruijn_lower_bound': True
}

# Output settings
LATEX_OUTPUT = os.path.join(OUTPUT_DIR, 'proofs')
PLOT_OUTPUT = os.path.join(OUTPUT_DIR, 'plots')
RESULTS_OUTPUT = os.path.join(OUTPUT_DIR, 'results')

# Create output directories if they don't exist
os.makedirs(LATEX_OUTPUT, exist_ok=True)
os.makedirs(PLOT_OUTPUT, exist_ok=True)
os.makedirs(RESULTS_OUTPUT, exist_ok=True)

# Logging
LOG_LEVEL = 'INFO'
LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'


def get_output_path(filename: str, output_type: str = 'plots') -> str:
    """
    Get full path for an output file.
    
    Args:
        filename: Name of the output file
        output_type: Type of output ('plots', 'latex', or 'results')
        
    Returns:
        Full path to the output file
    """
    if output_type == 'plots':
        return os.path.join(PLOT_OUTPUT, filename)
    elif output_type == 'latex':
        return os.path.join(LATEX_OUTPUT, filename)
    elif output_type == 'results':
        return os.path.join(RESULTS_OUTPUT, filename)
    else:
        return os.path.join(OUTPUT_DIR, filename)


if __name__ == "__main__":
    print("A112509 Research Configuration")
    print("=" * 50)
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Data directory: {DATA_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"\nBeam search config: {BEAM_SEARCH_CONFIG}")
