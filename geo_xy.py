import pdb
import logging
from typing import TypeAlias
import operator

Point: TypeAlias = tuple[float, float]
Line: TypeAlias = tuple[Point, Point]
Polygon: TypeAlias = list[Point]

def cut_polygon(poly:Polygon, position: float, dimension: int, direction='le') -> Polygon:
    # Cut a Polygon at position, dimension indicates X or Y (0th or 1st index).
    # If direction is '<=', we keep elements that are less than or equal to the position of the given dimension.

    outside_boundary = False
    idx_boundary_pairs = list()
    cur_boundary_pair = list()
    for n_idx in range(len(poly)):
        node = poly[n_idx]
        if getattr(operator, direction)(node[dimension], position):
            if outside_boundary:  # We are re-entering
                assert len(cur_boundary_pair) == 1
                cur_boundary_pair.append(n_idx)
                idx_boundary_pairs.append(cur_boundary_pair)
                cur_boundary_pair = list()
            outside_boundary = False
        else:
            if not outside_boundary:  # We are exiting the boundary!
                assert len(cur_boundary_pair) == 0
                cur_boundary_pair.append(n_idx)
            outside_boundary = True

    # Now we have a list of "boundary pairs," indicating the first
    # node out of the box and the first node back in, respectively.
    # We iterate over our boundary pairs, create new nodes at the
    # intersection of the boundary and crossing line, and delete
    # the nodes outside of the boundary!
    for boundary_pair in idx_boundary_pairs:
        logging.debug(f"Found boundary pair: {boundary_pair}")
    return poly

def is_point_in_box(p:Point, max_x_min_y:Point, min_x_max_y:Point) -> bool:
    x_range = sorted([abs(max_x_min_y[0]), abs(min_x_max_y[0])])
    y_range = sorted([abs(max_x_min_y[1]), abs(min_x_max_y[1])])
    return (abs(p[0]) >= x_range[0] and abs(p[0]) <= x_range[1] and
            abs(p[1]) >= y_range[0] and abs(p[1]) <= y_range[1])

    
def is_point_on_line(l: Line, p: Point) -> bool:
    return (p[0] <= max(l[0][0], l[1][0]) and
            # p.x <= max(l1.p1.x, l1.p2.x)
            p[0] >= min(l[0][0], l[1][0]) and
            # && p.x <= min(l1.p1.x, l1.p2.x) BUG in source
            p[1] <= max(l[0][1], l[1][1]) and
            # && (p.y <= max(l1.p1.y, l1.p2.y)
            p[1] >= min(l[0][1], l[1][1]))
            # && (p.y <= min(l1.p1.y, l1.p2.y) BUG in source


def clock_direction(a: Point, b: Point, c: Point) -> int:
    """Clockwise direction of three sequential points.
    Returns 0 for colinear, 2 for counter-clockwise, and 1 for clockwise.
    """
    val = ((b[1] - a[1]) * (c[0] - b[0]) -
           (b[0] - a[0]) * (c[1] - b[1]))
    if val == 0:
        # Colinear
        return 0
    elif val < 0:
        # counter-clockwise direction
        return 2
    #clockwise direction
    return 1

def do_lines_intersect(l1: Line, l2: Line) -> bool:
    # Direction of line one turning to starting point of line 2
    dir_l1_l2p1 = clock_direction(l1[0], l1[1], l2[0])
    # Direction of line one turning to ending point of line 2
    dir_l1_l2p2 = clock_direction(l1[0], l1[1], l2[1])
    # Direction of line two turning to starting point of line 1
    dir_l2_l1p1 = clock_direction(l2[0], l2[1], l1[0])
    # Direction of line two turning to ending point of line 1
    dir_l2_l1p2 = clock_direction(l2[0], l2[1], l1[1])    

    # When the two lines intersect
    if (dir_l1_l2p1 != dir_l1_l2p2 and
        dir_l2_l1p1 != dir_l2_l1p2):
        return True

    # When start of line 2 is on the line l1
    if (dir_l1_l2p1 == 0
        and is_point_on_line(l1, l2[0])):
        return True

    # When end of line 2 is on line l1
    if (dir_l1_l2p2 == 0
        and is_point_on_line(l1, l2[1])):
        return True

    # When start of line 1 is on line l2
    if (dir_l2_l1p1 == 0
        and is_point_on_line(l2, l1[0])):
        return True

    # When end of line 1 is on line l2
    if (dir_l2_l1p2 == 0
        and is_point_on_line(l2, l1[1])):
        return True
    return False

def check_inside(polygon: list[Point], probe: Point) -> bool:
    if len(polygon) < 3:
        return False
    probe_line = Line(probe, (9999, probe[1]))
    intersections = 0
    for i in range(1,len(polygon)):
        cur_side = Line((polygon[i-1], polygon[i]))
        if do_lines_intersect(cur_side, probe_line):
            if clock_direction(cur_side[0], probe, cur_side[1]) == 0:
                # If you picked a point that is on a side...
                return is_point_on_line(cur_side, probe)
            intersections += 1
    # If we intersected an odd number of sides, we were on the inside.
    return (intersections & 1) 
