"""Great-circle geodesy on a spherical Earth.

Distances and bearings between latitude/longitude points using the haversine formula and
the standard spherical trig relations: the great-circle distance, the initial (forward)
bearing along that path, the destination reached from a point given a bearing and
distance, and the cross-track distance of a point from a great-circle path. All angles in
degrees, distances in kilometres on a mean-radius sphere (good to ~0.3% versus the
ellipsoid). Pure standard library.
"""

import math

EARTH_RADIUS_KM = 6371.0088          # IUGG mean radius


def haversine_distance(lat1, lon1, lat2, lon2, radius=EARTH_RADIUS_KM):
    """Great-circle distance between two lat/lon points (haversine), in ``radius`` units.

    Numerically stable for small distances (unlike the spherical law of cosines). Inputs
    in degrees; returns kilometres by default.
    """
    p1 = math.radians(lat1)
    p2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlam / 2) ** 2
    return 2.0 * radius * math.asin(min(1.0, math.sqrt(a)))


def initial_bearing(lat1, lon1, lat2, lon2):
    """Initial (forward) bearing from point 1 to point 2 along the great circle, degrees.

    The compass bearing to steer at the start of the path, in ``[0, 360)`` clockwise from
    north. It changes along a great circle, so this is the *initial* heading.
    """
    p1 = math.radians(lat1)
    p2 = math.radians(lat2)
    dlam = math.radians(lon2 - lon1)
    y = math.sin(dlam) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dlam)
    return (math.degrees(math.atan2(y, x)) + 360.0) % 360.0


def destination_point(lat, lon, bearing, distance, radius=EARTH_RADIUS_KM):
    """Destination ``(lat, lon)`` reached from a start point on a given bearing/distance.

    Follows the great circle from ``(lat, lon)`` heading ``bearing`` degrees for
    ``distance`` (same units as ``radius``). Returns degrees; the inverse of
    :func:`haversine_distance` / :func:`initial_bearing`.
    """
    d = distance / radius
    theta = math.radians(bearing)
    p1 = math.radians(lat)
    l1 = math.radians(lon)
    p2 = math.asin(math.sin(p1) * math.cos(d) + math.cos(p1) * math.sin(d) * math.cos(theta))
    l2 = l1 + math.atan2(math.sin(theta) * math.sin(d) * math.cos(p1),
                         math.cos(d) - math.sin(p1) * math.sin(p2))
    lon2 = (math.degrees(l2) + 540.0) % 360.0 - 180.0     # normalize to [-180, 180)
    return math.degrees(p2), lon2


def cross_track_distance(lat, lon, lat1, lon1, lat2, lon2, radius=EARTH_RADIUS_KM):
    """Signed distance of point ``(lat, lon)`` from the great circle through 1 and 2.

    Positive when the point lies to the *right* of the path direction (1 -> 2), negative
    to the left; magnitude is the perpendicular great-circle distance. Useful for "how far
    off the route am I".
    """
    d13 = haversine_distance(lat1, lon1, lat, lon, radius) / radius
    theta13 = math.radians(initial_bearing(lat1, lon1, lat, lon))
    theta12 = math.radians(initial_bearing(lat1, lon1, lat2, lon2))
    return math.asin(math.sin(d13) * math.sin(theta13 - theta12)) * radius
