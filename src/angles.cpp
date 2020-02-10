#include "angles.h"

// convert an angle from degree to radian
double degree_to_radian(double angle) { return angle * (PI / 180); }

// convert an angle from radian to degree, convert negative to positive degrees
double radian_to_degree(double angle) {
    angle *= (180 / PI);
    if (angle < 0) return angle + 360;
    return angle;
}
