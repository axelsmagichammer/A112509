"""Analyze common structure across multiple n values."""

import json
import os
from collections import defaultdict


def find_runs(s):
    """Find all runs of 1s and 0s in a binary string.
    Returns list of (char, length) tuples.
    """
    if not s:
        return []
    
    runs = []
    current_char = s[0]
    current_len = 1
    
    for i in range(1, len(s)):
        if s[i] == current_char:
            current_len += 1
        else:
            runs.append((current_char, current_len))
            current_char = s[i]
            current_len = 1
    
    runs.append((current_char, current_len))
    return runs


def extract_structure(s, max_k=15):
    """Extract block structure: K blocks of 1s separated by runs of 0s."""
    runs = find_runs(s)
    
    blocks = []
    seps = []
    
    i = 0
    while i < len(runs):
        char, length = runs[i]
        if char == '1':
            blocks.append(length)
            i += 1
            # Look for separator (0s)
            if i < len(runs) and runs[i][0] == '0':
                seps.append(runs[i][1])
                i += 1
            else:
                break
        else:
            break
    
    return blocks, seps


def main():
    cache_path = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'cached_results.json')
    with open(cache_path, 'r') as f:
        cache = json.load(f)
    
    target_n_values = [137]
    
    # Collect all structures for each n
    all_structures = {}
    for n in target_n_values:
        if str(n) not in cache:
            print(f"No data for n={n}")
            continue
        
        strings = cache[str(n)]['optimal_strings']
        structures = []
        for s in strings:
            blocks, seps = extract_structure(s)
            structures.append((blocks, seps))
        
        all_structures[n] = structures
        print(f"n={n}: {len(strings)} optimal strings")
    
    print(f"\n{'='*80}")
    print(f"100% COVERAGE ANALYSIS FOR EACH N")
    print(f"{'='*80}")
    
    # Analyze each n individually for 100% coverage
    for n in target_n_values:
        if n not in all_structures:
            continue
        
        print(f"\n{'='*80}")
        print(f"n={n}")
        print(f"{'='*80}")
        
        structures = all_structures[n]
        num_strings = len(structures)
        
        # Find common separator prefix across ALL strings for this n
        all_seps = [seps for blocks, seps in structures]
        min_sep_len = min(len(seps) for seps in all_seps)
        
        common_sep_prefix = []
        for i in range(min_sep_len):
            values_at_i = set(seps[i] for seps in all_seps)
            if len(values_at_i) == 1:
                common_sep_prefix.append(list(values_at_i)[0])
            else:
                break
        
        K_common = len(common_sep_prefix) + 1  # +1 because seps are between blocks
        
        print(f"Common separator prefix (100% coverage): {common_sep_prefix}")
        print(f"K_common = {K_common} blocks")
        
        # Find block ranges for the common prefix blocks
        block_ranges = []
        for i in range(K_common):
            sizes = [blocks[i] for blocks, seps in structures if len(blocks) > i]
            if sizes:
                min_size = min(sizes)
                max_size = max(sizes)
                block_ranges.append((min_size, max_size))
                if min_size == max_size:
                    print(f"  Block {i+1}: {min_size} (fixed)")
                else:
                    print(f"  Block {i+1}: {min_size}–{max_size}")
        
        # Count total 1s in entire strings
        total_1s_list = [sum(blocks) for blocks, seps in structures]
        min_1s = min(total_1s_list)
        max_1s = max(total_1s_list)
        density_min = min_1s / n
        density_max = max_1s / n
        
        # Also show 1s in common prefix for reference
        prefix_1s_list = [sum(blocks[:K_common]) for blocks, seps in structures]
        min_prefix_1s = min(prefix_1s_list)
        max_prefix_1s = max(prefix_1s_list)
        
        print(f"\nTotal 1s in entire string: {min_1s}–{max_1s}")
        print(f"Total 1s in first {K_common} blocks: {min_prefix_1s}–{max_prefix_1s}")
        print(f"Density: {density_min:.3f}–{density_max:.3f}")
        
        # Calculate combinations
        combos = 1
        for lo, hi in block_ranges:
            combos *= (hi - lo + 1)
        print(f"Template combinations: {combos}")
        
        # Count 1-blocks in entire strings
        one_blocks_list = [len(blocks) for blocks, seps in structures]
        min_one_blk = min(one_blocks_list)
        max_one_blk = max(one_blocks_list)
        zero_blocks_list = [len(seps) for blocks, seps in structures]
        min_zero_blk = min(zero_blocks_list)
        max_zero_blk = max(zero_blocks_list)
        print(f"\n1-block counts: {min_one_blk}\u2013{max_one_blk}")
        print(f"0-block counts: {min_zero_blk}\u2013{max_zero_blk}")

        # Output for learned_bounds.json
        print(f"\nFor learned_bounds.json:")
        print(f'  "target_n": {n},')
        print(f'  "K_common": {K_common},')
        print(f'  "common_seps": {common_sep_prefix},')
        print(f'  "block_ranges": {block_ranges},')
        print(f'  "min_total_1s": {min_1s},')
        print(f'  "max_total_1s": {max_1s},')
        print(f'  "min_one_blocks": {min_one_blk},')
        print(f'  "max_one_blocks": {max_one_blk}')
    
    print(f"\n{'='*80}")
    print(f"CROSS-N COMMON SEPARATOR PREFIX")
    print(f"{'='*80}")
    
    # Find common separator prefix across ALL n values and ALL strings
    all_seps_all_n = []
    for n in target_n_values:
        if n not in all_structures:
            continue
        for blocks, seps in all_structures[n]:
            all_seps_all_n.append(seps)
    
    # Find the minimum length
    min_sep_len_all = min(len(seps) for seps in all_seps_all_n)
    
    # Find common prefix across ALL strings from ALL n
    cross_n_common_sep_prefix = []
    for i in range(min_sep_len_all):
        values_at_i = set(seps[i] for seps in all_seps_all_n)
        if len(values_at_i) == 1:
            cross_n_common_sep_prefix.append(list(values_at_i)[0])
        else:
            break
    
    K_cross_n = len(cross_n_common_sep_prefix) + 1
    
    if cross_n_common_sep_prefix:
        print(f"\nCommon separator prefix across ALL n (100% of all strings):")
        print(f"  {cross_n_common_sep_prefix}")
        print(f"  Length: {len(cross_n_common_sep_prefix)}")
        print(f"  K_common = {K_cross_n} blocks")
        
        # Find block ranges for these K blocks across all n
        print(f"\nBlock ranges for first {K_cross_n} blocks (union across all n):")
        cross_n_block_ranges = []
        for block_idx in range(K_cross_n):
            all_values = []
            for n in target_n_values:
                if n not in all_structures:
                    continue
                for blocks, seps in all_structures[n]:
                    # Check if this structure has the common separator prefix
                    if len(seps) >= len(cross_n_common_sep_prefix):
                        if list(seps[:len(cross_n_common_sep_prefix)]) == cross_n_common_sep_prefix:
                            if len(blocks) > block_idx:
                                all_values.append(blocks[block_idx])
            
            if all_values:
                min_val = min(all_values)
                max_val = max(all_values)
                cross_n_block_ranges.append((min_val, max_val))
                if min_val == max_val:
                    print(f"  Block {block_idx+1}: {min_val} (fixed)")
                else:
                    print(f"  Block {block_idx+1}: {min_val}–{max_val}")
        
        # Calculate combinations
        combos = 1
        for lo, hi in cross_n_block_ranges:
            combos *= (hi - lo + 1)
        print(f"\nTemplate combinations: {combos}")
        
        # Show which blocks are fixed
        fixed_blocks = []
        variable_blocks = []
        for i, (lo, hi) in enumerate(cross_n_block_ranges):
            if lo == hi:
                fixed_blocks.append(i+1)
            else:
                variable_blocks.append(i+1)
        
        print(f"Fixed blocks: {fixed_blocks}")
        print(f"Variable blocks: {variable_blocks}")
    else:
        print(f"\nNo common separator prefix found across all n values")


if __name__ == '__main__':
    main()
