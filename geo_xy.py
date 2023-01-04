import pdb
import logging
from typing import TypeAlias
import operator
from math import radians, cos, sin, asin, sqrt

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
    iter_poly_idx = 0  # We iterate over the source polygon to build the amended polygon.
    amended_polygon = list()
    for boundary_pair in idx_boundary_pairs:
        # Add a vertex on the box where we exit.
        post_exit_idx = boundary_pair[0]
        if post_exit_idx > 0:
            pre_exit_idx = boundary_pair[0] - 1
        else:
            pre_exit_idx = len(poly) - 2  # We end on the initial point, so we jump 2

        post_reentry_idx = boundary_pair[1]
        if post_reentry_idx > 0:
            pre_reentry_idx = boundary_pair[1] - 1
        else:
            pre_reentry_idx = len(poly) - 2

        new_exit_coord = line_at(
            (poly[pre_exit_idx], poly[post_exit_idx]),
            position, dimension)
        new_reentry_coord = line_at(
            (poly[pre_reentry_idx], poly[post_reentry_idx]),
            position, dimension)

        # Populate nodes up to our most recent exit point into the amended polygon.
        while iter_poly_idx <= pre_exit_idx:
            amended_polygon.append(poly[iter_poly_idx])
            iter_poly_idx += 1
        amended_polygon.append(new_exit_coord)
        amended_polygon.append(new_reentry_coord)
        iter_poly_idx = post_reentry_idx
    while iter_poly_idx < len(poly):
        amended_polygon.append(poly[iter_poly_idx])
        iter_poly_idx += 1
    return amended_polygon

def line_at(l1:Line, position: float, dimension: int) -> Point:
    # Compute the slope
    if (l1[1][0] - l1[0][0]) == 0:
        # The line is vertical!
        if dimension == 0:
            # We are enforcing the X value, but the line is fixed as X. 
            pdb.set_trace()
        else:
            return (l1[1][0], position)
    if (l1[1][1] - l1[0][1]) == 0:
        # The line is horizontal!
        if dimension == 1:
            # We are enforcing the Y value, but the line is fixed as Y.
            pdb.set_trace()
        else:
            return (position, l1[1][1])
    m = (l1[1][1] - l1[0][1]) / (l1[1][0] - l1[0][0])
    b = l1[0][1] - m * l1[0][0]
    output_point = (0,0)
    if dimension == 0:  # X is fixed
        target_variable = m * position + b
        output_point = (position, target_variable)
    else:  # Y is fixed
        target_variable = (position - b) / m
        output_point = (target_variable, position)
    return output_point

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
    probe_line = Line((probe, (9999, probe[1])))
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

def haversine(line: Line):
    """
    Calculate the great circle distance in kilometers between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [line[0][0], line[0][1], line[1][0], line[0][1]])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r
