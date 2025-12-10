"""
This script is based on the original implementation by Ondřej Huvar:
https://github.com/sybila/biodivine-algo-smt-inference/blob/main/scripts/check_fixed_points.py

The original work is licensed under the MIT License:
https://github.com/sybila/biodivine-algo-smt-inference/blob/main/LICENSE

Modifications: minor adjustments made by Timotej Šváby, 2025.
"""


from biodivine_aeon import (AsynchronousGraph, BooleanNetwork,
    ColoredVertexSet, FixedPoints)

from typing import Any


def collect_fp_combinations(fixed_points: ColoredVertexSet, print_all: bool):
    """Iterate over all colors and collect unique fixed-point combinations.
    This is a very naive way to do this, but good enough to explore toy models.

    If flag `print_all` is set, print fixed-points for all interpretations.
    """
    fixed_point_colors  = fixed_points.colors()
    unique_combinations = set()
    if print_all:
        print("Fixed point combinations per model color:")

    count = 1
    while not fixed_point_colors.is_empty():
        # Pick one color and restrict to its fixed points
        color_singleton = fixed_point_colors.pick_singleton()
        fps_single_color = fixed_points.intersect_colors(color_singleton)

        # Extract vertex values as binary strings
        fps_values = [tuple(v.values()) for v in fps_single_color.vertices()]
        binary_vectors = [''.join('1' if x else '0' for x in fp) for fp in fps_values]
        # Sets cant deal with lists, so we just make it a tuple
        unique_combinations.add(tuple(binary_vectors))

        if print_all:
            print(count, binary_vectors)
            print("\t->", next(iter(color_singleton)))  # print the only color in the set

        fixed_point_colors = fixed_point_colors.minus(color_singleton)
        count += 1

    if print_all:
        print("------")
    return unique_combinations

def print_fixed_points_on_new_stg(bn: Any) -> None:
    bn = bn.infer_valid_graph()
    stg = AsynchronousGraph(bn)

    print(f"Total colors: {stg.mk_unit_colors().cardinality()}")
    print("------")

    # Compute fixed-points across all interpretations (colors)
    fixed_points = FixedPoints.symbolic(stg)

    print(f"Total colored fixed points: {fixed_points.cardinality()}")
    print(f"Total fixed point states: {fixed_points.vertices().cardinality()}")
    print(f"Total fixed point colors: {fixed_points.colors().cardinality()}")
    print("------")

    print("Raw fixed point vertices projection (across all colors):")
    for fp in fixed_points.vertices():
        print(fp)
    print("------")

    # Iterate over all interpretations and collect the unique combinations
    # Print all colors and their fixed-points while iterating
    unique_combinations = collect_fp_combinations(fixed_points, print_all=True)

    print("Unique fixed point combinations:")
    for fp_combination in unique_combinations:
        print(list(fp_combination))