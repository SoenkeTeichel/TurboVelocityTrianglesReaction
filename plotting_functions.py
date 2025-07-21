# -*- coding: utf-8 -*-
"""
Created on Sat Jun 28 08:17:54 2025

@author: Soenke
"""
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Arc, Circle
import os
os.environ["LIBGL_ALWAYS_SOFTWARE"] = "1"  # Force software rendering
import CoolProp.CoolProp as CP


def arrow_data(EndPoint, dx, dy):
    x_start = EndPoint[0] - dx
    y_start = EndPoint[1] + dy
    x_end = EndPoint[0]
    y_end = EndPoint[1]
    
    # Position arrowhead at 95% of length
    arrow_x = x_start + 0.95 * dx
    arrow_y = y_start - 0.95 * dy

    angle = 180 + np.degrees(np.arctan2(dx, dy))
    
    return x_start,y_start,x_end,y_end,arrow_x,arrow_y,angle



def draw_arrow(ax, EndPoint, dx, dy, label, color, linestyle,**kwargs): 
    """Draws a vector line with a rotated arrowhead at 95% and a label"""

    x_start,y_start,x_end,y_end,arrow_x,arrow_y,angle = arrow_data(EndPoint, dx, dy)
    
    # Line from tail to tip
    line = ax.plot([x_start, x_end], [y_start, y_end], color=color, linestyle=linestyle, linewidth=2)
    arrowhead = ax.plot(arrow_x, arrow_y, marker=(3, 0, angle), color=color, markersize=10)

    # Label position
    if kwargs.get('horizontalalignment',False):
        horizontalalignment = kwargs['horizontalalignment']
    else:
        horizontalalignment='left'
        
    text = ax.text((x_start + 0.5*dx) , y_start - 0.5*dy, label, color=color, size='large',horizontalalignment=horizontalalignment)

    return [line[0], arrowhead[0], text]


    

def redraw_arrow(obj, EndPoint, dx, dy):
    """Redraws line, moves arrowhead to 95% of the line, updates label"""
    line, arrowhead, text = obj

    x_start,y_start,x_end,y_end,arrow_x,arrow_y,angle = arrow_data(EndPoint, dx, dy)

    # Redraw main line
    line.set_xdata([x_start, x_end])
    line.set_ydata([y_start, y_end])

    arrowhead.set_data([arrow_x], [arrow_y])
    arrowhead.set_marker((3, 0, angle))

    # Move label
    text.set_x(x_start + 0.5 * dx)
    text.set_y(y_start - 0.5 * dy)
    
    return


def angle_text(u,c1,c2,w1,w2,cx,cu1,cu2,wu1,wu2,alpha1,beta1,alpha2,beta2):
    text = (
        f"$\\alpha_1$ = {alpha1:.1f}°;  $\\beta_1$ = {beta1:.1f}°\n"
        f"$\\alpha_2$ = {alpha2:.1f}°;  $\\beta_2$ = {beta2:.1f}°\n"
        f"Schaufel Umlenkung:  $\\alpha_2 - \\alpha_1$ = {alpha2-alpha1:.2f}\n"
        f"Schaufel Umlenkung:  $\\beta_2 - \\beta_1$ = {beta2-beta1:.2f}\n\n"
        f"Euler Arbeit:  $u (c_{{u2}} - c_{{u1}})$ = {u*(cu2-cu1)/1000:.2f}kJ\n"
        f"$\\frac{{w_2}}{{w_1}}$ = {w2/w1:.2f}\n"
    )
    return text




def compute_arc_geometry(P0, beta1, beta2, h, eps=1e-6):
    """
    Given:
      P0 = (x0,y0)
      beta1, beta2 in degrees (0° = down, +→ right)
      vertical rise h
    Returns:
      P1      = end-point
      C       = center (or None if straight‐line)
      R       = signed radius (np.inf if straight‐line)
      theta1, theta2 = start/end angles for Arc (None if straight‐line)
      is_line = True if we must draw a straight segment
    """
    x0, y0 = P0
    b1, b2 = np.radians(beta1), np.radians(beta2)

    # unit tangents
    T1 = np.array([ np.sin(b1), -np.cos(b1) ])
    T2 = np.array([ np.sin(b2), -np.cos(b2) ])
    # normals (90° CCW of tangent)
    N1 = np.array([ T1[1], -T1[0] ])
    N2 = np.array([ T2[1], -T2[0] ])

    # detect straight‐line cases:
    denom = N1[1] - N2[1]
    # either circle‐equation blows up, or betas identical
    if abs(denom) < eps or abs(beta2 - beta1) < eps:
        if abs(np.cos(b1)) < eps:
            # vertical tangent: just move vertically
            dx = 0
        else:
            dx = -h * np.tan(b1)
        P1 = (x0 + dx, y0 + h)
        return P1, None, np.inf, None, None, True

    # else true circular arc
    R  = h / denom
    P1 = (x0 + R*(N1[0] - N2[0]), y0 + h)
    C  = np.array([x0, y0]) + R * N1

    phi0 = np.degrees(np.arctan2(y0 - C[1], x0 - C[0])) % 360
    phi1 = np.degrees(np.arctan2(P1[1] - C[1], P1[0] - C[0])) % 360
    if R < 0:
        phi0, phi1 = phi1, phi0
    dphi = (phi1 - phi0) % 360
    if dphi > 180:
        theta1, theta2 = phi1, phi0
    else:
        theta1, theta2 = phi0, phi1

    return P1, C, R, theta1, theta2, False




def draw_camber_circles(fig,ax,is_line,arc1_patch,arc2_patch,line_patch,P0,P1,C,R,t1,t2,h,circle_patches):
    """Show either the arc or the straight‐line + update markers."""
    
    if is_line:
        # straight‐line along tangent beta1
        arc1_patch.set_visible(False)
        line_patch.set_data([P0[0], P1[0]], [P0[1], P1[1]])
        line_patch.set_visible(True)
    else:
        # circular arc
        arc1_patch.center = tuple(C)
        arc1_patch.width = 2 * abs(R)
        arc1_patch.height = 2 * abs(R)
        arc1_patch.theta1 = t1 % 360  # Ensure angles are within [0, 360]
        arc1_patch.theta2 = t2 % 360  # Ensure angles are within [0, 360]
        arc1_patch.set_visible(True)
        line_patch.set_visible(False)

        shift_x = h  # Offset for the second arc
        arc2_patch.center = (C[0] + shift_x, C[1])
        arc2_patch.width = arc1_patch.width
        arc2_patch.height = arc1_patch.height
        arc2_patch.theta1 = t1 % 360  # Ensure angles are within [0, 360]
        arc2_patch.theta2 = t2 % 360  # Ensure angles are within [0, 360]
        arc2_patch.set_visible(True)
        line_patch.set_visible(False)

    # Remove existing circles
    for c in circle_patches:
        c.remove()
    circle_patches.clear()

    if not is_line:
        C1 = np.array(arc1_patch.center)
        R1 = abs(R)

        C2 = np.array(arc2_patch.center)
        R2 = abs(R)

        # Calculate arc length and distribute points equidistantly along arc1
        theta_start_rad = np.radians(arc1_patch.theta1)
        theta_end_rad = np.radians(arc1_patch.theta2)

        # Handle cases where theta_end < theta_start (crossing 360° boundary)
        if theta_end_rad < theta_start_rad:
            theta_end_rad += 2 * np.pi

        arc_length = R1 * (theta_end_rad - theta_start_rad)  # Arc length formula
        num_points = 10  # Number of equidistant points
        relative_steps = np.linspace(0, arc_length, num_points)  # Relative length steps
        theta_range = theta_start_rad + relative_steps / R1  # Convert length steps to angles

        arc1_points = np.array([C1[0] + R1 * np.cos(theta_range),
                                C1[1] + R1 * np.sin(theta_range)]).T

        for point_on_arc1 in arc1_points:
            # Compute the closest point on arc2
            direction = (point_on_arc1 - C2) / np.linalg.norm(point_on_arc1 - C2)
            closest_point_on_arc2 = C2 + R2 * direction

            # Compute the center and diameter of the circle
            line_center = (point_on_arc1 + closest_point_on_arc2) / 2
            line_length = np.linalg.norm(point_on_arc1 - closest_point_on_arc2)
        
            circle = Circle(line_center, line_length / 2, edgecolor='green', facecolor='none', linestyle='--', linewidth=1)
            ax.add_patch(circle)
            circle_patches.append(circle)
   

    fig.canvas.draw_idle()
    return circle_patches