"""Known values of OEIS A112509 and A156025."""

from typing import Optional

# Source: https://oeis.org/A112509
# Computed by Martin Fuller up to n=80

KNOWN_VALUES = [
    1, 3, 5, 7, 10, 13, 17, 22, 27, 33, 40, 47, 55, 64, 73, 83,
    94, 106, 118, 131, 145, 160, 176, 192, 209, 227, 246, 265,
    285, 306, 328, 351, 375, 399, 424, 450, 477, 504, 532, 561,
    591, 622, 654, 687, 720, 754, 789, 825, 862, 899, 937, 977,
    1017, 1058, 1100, 1143, 1186, 1230, 1275, 1321, 1368, 1416, 
    1465, 1514, 1564, 1615, 1667, 1720, 1774, 1828, 1884, 1941, 
    1998, 2056, 2115, 2175, 2236, 2298, 2360, 2423
]

# Source: https://oeis.org/A156025
# Number of n-bit numbers achieving the maximum A112509(n) distinct integers
A156025_VALUES = [
    2, 1, 1, 3, 2, 6, 5, 1, 4, 5, 2, 8, 10, 4, 16, 22,
    12, 2, 10, 19, 17, 7, 1, 5, 9, 7, 2, 11, 24, 28,
    20, 9, 2, 10, 18, 14, 4, 26, 68, 94, 78, 44, 18, 4,
    22, 46, 46, 22, 4, 29, 104, 4, 20, 36, 28, 8, 52, 140,
    202, 168, 80, 20, 2, 14, 40, 60, 50, 22, 4, 152, 23, 2,
    14, 40, 60, 50, 22, 4, 32, 104
]

def get_known_values(max_n: Optional[int] = None) -> list:
    """
    Get known values of A112509.
    
    Args:
        max_n: Maximum n to return (None for all known values)
        
    Returns:
        List of known values
    """
    if max_n is None:
        return KNOWN_VALUES.copy()
    return KNOWN_VALUES[:min(max_n, len(KNOWN_VALUES))]


def get_value(n: int) -> int:
    """
    Get a(n) if known, otherwise raise ValueError.
    
    Args:
        n: Index (1-indexed)
        
    Returns:
        a(n)
    """
    if n < 1 or n > len(KNOWN_VALUES):
        raise ValueError(f"Value a({n}) not in known data (known up to n={len(KNOWN_VALUES)})")
    return KNOWN_VALUES[n - 1]


def max_known_n() -> int:
    """Return the maximum n for which we have a known value."""
    return len(KNOWN_VALUES)


if __name__ == "__main__":
    print(f"Known values of A112509 (n=1 to {max_known_n()}):")
    print("=" * 50)
    
    for i in range(min(20, len(KNOWN_VALUES))):
        print(f"a({i+1:2d}) = {KNOWN_VALUES[i]:4d}")
    
    if len(KNOWN_VALUES) > 20:
        print("...")
        print(f"a({len(KNOWN_VALUES)}) = {KNOWN_VALUES[-1]}")
