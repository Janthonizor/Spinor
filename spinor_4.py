import numpy as np
import math as m
import pyvista as pv
import os
import imageio
from PIL import Image
import re

def sphere2cart(r, theta, phi):
    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)
    return x, y, z

def points(theta, phi, r_inner, r_outer, r_static, r_res):
    r_values = np.linspace(r_inner, r_outer + r_static, r_res)
    total_combinations = r_res * len(theta)
    output_points = np.zeros((total_combinations, 3))
    index = 0
    for r in r_values:
        for t, p in zip(theta, phi):
            output_points[index] = [r, t, p]
            index += 1
    output_points = output_points.reshape(total_combinations, 3)
    r_combined = output_points[:, 0]
    theta_combined = output_points[:, 1]
    phi_combined = output_points[:, 2]
    x, y, z = sphere2cart(r_combined, theta_combined, phi_combined)

    cart_points = np.vstack([x, y, z]).T

    return cart_points, r_combined.ravel()

def rotatez(points, delta):
    R_z = np.array([[np.cos(delta), -np.sin(delta), 0],
                [np.sin(delta),  np.cos(delta), 0],
                [0,              0,            1]])
    rot_points = np.dot(points, R_z.T)
    return rot_points

def twist(points, r, r_inner, r_outer):
    mask = r <= r_outer
    phi = np.zeros_like(r)
    phi[mask] = m.pi/2 * (1 - np.cos((r_outer - r[mask]) / (r_outer - r_inner) * m.pi))
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    xp = np.copy(x)  
    zp = np.copy(z)
    xp[mask] = x[mask] * np.cos(phi[mask]) + z[mask] * np.sin(phi[mask])
    zp[mask] = -x[mask] * np.sin(phi[mask]) + z[mask] * np.cos(phi[mask])
    npoints = np.vstack((xp, y, zp)).T
    return npoints

   
resolution = 1000
r_inner = 1
r_outer = 3
r_flat = 1
theta = m.pi * np.concatenate((np.linspace(0,0,12), np.linspace(1, 1, 12)))
phi = np.concatenate((np.linspace(0,m.pi, 12), np.linspace(m.pi, 0, 12)))
cart_points, r = points(theta, phi, 1, 3, 1, 1000)

num_frames = 800
rotation_proportion = 1/2
num_rotations = 3
plotter = pv.Plotter(off_screen=True)
for frame in range(num_frames):
    frame_points = cart_points
    print(frame)
    delta = 2 * m.pi * frame / num_frames * num_rotations
    frame_points = rotatez(cart_points, delta)
    frame_points = twist(frame_points, r, 1,3)
    frame_points = rotatez(frame_points, -(1+rotation_proportion) * delta)
    plotter.show_axes()
    mesh = pv.PolyData(frame_points)
    norm_r = np.zeros_like(r)
    mask = r <=r_outer
    norm_r[mask] = (np.max(r[mask])-r[mask]) / (np.max(r[mask]) - np.min(r[mask]))
    mesh['radius'] = norm_r
    plotter.add_mesh(mesh,scalars = 'radius',  cmap = 'viridis' , show_scalar_bar=False, point_size = 6, render_points_as_spheres=True)
    plotter.camera.elevation = -10
    
    plotter.screenshot(f"frame_{frame}")
    plotter.clear()

image_dir = ''

image_files = [f for f in os.listdir(image_dir) if f.endswith('.png')]
image_files.sort(key=lambda x: int(re.search(r'\d+', x).group()))
images = []

for filename in image_files:
    image_path = os.path.join(image_dir, filename)
    img = Image.open(image_path)
    images.append(img.copy()) 
    img.close()

output_video = ''

fps = 30
#frame_duration = 1/fps
#imageio.mimsave(output_gif, images, fps = fps, duration = frame_duration)
images.extend(images)
imageio.mimsave(output_video, images, format='FFMPEG', fps=fps)
print(f'video saved with {fps} fps.')
