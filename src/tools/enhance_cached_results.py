"""
Enhance cached_results.json with additional OEIS sequences.

This script adds the following fields to each entry:
- A112510: Smallest n-bit number achieving A112509(n) (decimal)
- A112511: Largest n-bit number achieving A112509(n) (decimal)  
- A156022: Maximum number of positive integers (A112509(n)-1 for n>=2, 1 for n=1)
- A156023: n*(n+1)/2 - A112509(n)
- A156024: n*(n+1)/2 - A156022(n)
"""

import json
from pathlib import Path


def enhance_cached_results():
    """Add A112510, A112511, and A156023 values to cached_results.json."""
    
    # Load existing data
    data_dir = Path(__file__).parent.parent.parent / "data"
    cache_file = data_dir / "cached_results.json"
    
    print(f"Loading {cache_file}...")
    with open(cache_file, 'r') as f:
        data = json.load(f)
    
    print(f"Found {len(data)} entries")
    print("Enhancing with A112510, A112511, A156022, A156023, and A156024 values...\n")
    
    # Process each n
    results = {}
    for n_str, entry in data.items():
        n = int(n_str)
        a_n = entry["a(n)"]
        optimal_strings = entry["optimal_strings"]
        
        # Convert binary strings to decimal integers
        decimal_values = [int(binary_str, 2) for binary_str in optimal_strings]
        
        # A112510: Smallest n-bit number
        min_decimal = min(decimal_values)
        
        # A112511: Largest n-bit number
        max_decimal = max(decimal_values)
        
        # A156022: Maximum number of positive integers
        # For n=1: A156022(1) = 1, For n>=2: A156022(n) = A112509(n) - 1
        a156022 = 1 if n == 1 else a_n - 1
        
        # A156023: n*(n+1)/2 - A112509(n)
        a156023 = n * (n + 1) // 2 - a_n
        
        # A156024: n*(n+1)/2 - A156022(n)
        a156024 = n * (n + 1) // 2 - a156022
        
        # Create enhanced entry
        results[n_str] = {
            "a(n)": a_n,
            "num_optimal": entry["num_optimal"],
            "optimal_strings": optimal_strings,
            "A112510_min_decimal": min_decimal,
            "A112511_max_decimal": max_decimal,
            "A156022_positive_only": a156022,
            "A156023_complement": a156023,
            "A156024_complement_positive": a156024
        }
        
        # Print first few and last few as examples
        if n <= 5 or n >= 89:
            min_binary = optimal_strings[decimal_values.index(min_decimal)]
            max_binary = optimal_strings[decimal_values.index(max_decimal)]
            a156022 = 1 if n == 1 else a_n - 1
            a156024 = n * (n + 1) // 2 - a156022
            print(f"n={n:2d}: a(n)={a_n:4d}, "
                  f"min={min_decimal:20d} ({min_binary}), "
                  f"max={max_decimal:20d} ({max_binary}), "
                  f"A156022={a156022:4d}, A156023={a156023:4d}, A156024={a156024:4d}")
    
    # Save enhanced data atomically (write to temp then rename, so a crash
    # mid-write never truncates the original file)
    print(f"\nSaving enhanced data to {cache_file}...")
    import tempfile
    tmp_path = cache_file.with_suffix('.tmp')
    try:
        with open(tmp_path, 'w') as f:
            json.dump(results, f, indent=2)
        tmp_path.replace(cache_file)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise
    
    print("Done!")
    
    # Print summary statistics
    print("\n" + "="*70)
    print("Summary Statistics:")
    print("="*70)
    
    # Show A112510 sequence (first 20 values)
    print("\nA112510 (Smallest n-bit number) - first 20 values:")
    a112510_values = [results[str(n)]["A112510_min_decimal"] for n in range(1, min(21, len(results)+1))]
    print(", ".join(str(v) for v in a112510_values))
    
    # Show A112511 sequence (first 20 values)
    print("\nA112511 (Largest n-bit number) - first 20 values:")
    a112511_values = [results[str(n)]["A112511_max_decimal"] for n in range(1, min(21, len(results)+1))]
    print(", ".join(str(v) for v in a112511_values))
    
    # Show A156022 sequence (first 20 values)
    print("\nA156022 (Positive integers only) - first 20 values:")
    a156022_values = [results[str(n)]["A156022_positive_only"] for n in range(1, min(21, len(results)+1))]
    print(", ".join(str(v) for v in a156022_values))
    
    # Show A156023 sequence (first 20 values)
    print("\nA156023 (Complement: n*(n+1)/2 - a(n)) - first 20 values:")
    a156023_values = [results[str(n)]["A156023_complement"] for n in range(1, min(21, len(results)+1))]
    print(", ".join(str(v) for v in a156023_values))
    
    # Show A156024 sequence (first 20 values)
    print("\nA156024 (Complement positive: n*(n+1)/2 - A156022(n)) - first 20 values:")
    a156024_values = [results[str(n)]["A156024_complement_positive"] for n in range(1, min(21, len(results)+1))]
    print(", ".join(str(v) for v in a156024_values))
    
    print("\n" + "="*70)
    print("Field descriptions:")
    print("  - a(n): Maximum number of distinct integers (A112509)")
    print("  - num_optimal: Count of n-bit strings achieving a(n)")
    print("  - optimal_strings: All binary strings achieving a(n)")
    print("  - A112510_min_decimal: Smallest decimal value achieving a(n)")
    print("  - A112511_max_decimal: Largest decimal value achieving a(n)")
    print("  - A156022_positive_only: Maximum positive integers (A112509-1 for n>=2, 1 for n=1)")
    print("  - A156023_complement: n*(n+1)/2 - a(n)")
    print("  - A156024_complement_positive: n*(n+1)/2 - A156022(n)")
    print("="*70)


if __name__ == "__main__":
    enhance_cached_results()
