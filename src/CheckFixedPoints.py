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

    print("Fixed point vertices turned on projection (across all colors):")
    for fp in fixed_points.vertices():
        fp_vertices_on = set()
        fp_renamed = {bn.get_variable_name(variable): binary_state for variable, binary_state in fp.items()}

        for vertex, binary_state in fp_renamed.items():
            if binary_state == 1:
                fp_vertices_on.add(vertex)

        print(fp_vertices_on)
    print("------")
