import numpy as np
import math as m
import pyvista as pv
import os
import imageio
from PIL import Image
import re

def sphere2cart(r, theta, phi):
    # too easy
    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)
    return x, y, z

def points(theta, phi, r_inner, r_outer, r_static, r_res):
    #theta phi 1d numpy
    #r_static is the length of the still fibres outside of r_outer
    r_values = np.linspace(r_inner, r_outer + r_static, r_res)
    # create combinations of each r value with each theta[i], phi[i]
    total_combinations = r_res * len(theta)
    output_points = np.zeros((total_combinations, 3))
    index = 0
    # sort the values by r as they are added to the output matrix
    for r in r_values:
        for t, p in zip(theta, phi):
            output_points[index] = [r, t, p]
            index += 1
    output_points = output_points.reshape(total_combinations, 3)
    r_combined = output_points[:, 0]
    theta_combined = output_points[:, 1]
    phi_combined = output_points[:, 2]
    # convert to cartesian for rotating and plotting
    x, y, z = sphere2cart(r_combined, theta_combined, phi_combined)

    cart_points = np.vstack([x, y, z]).T

    #return the cartesian fibre points and the corresponding r values so they only have to be calculated once
    return cart_points, r_combined
def grid(r_outer, r_inner, rad_res, line_res):
    # theta and phi arrays
    theta = m.pi * np.concatenate((np.linspace(0, 0, int(line_res / 2)), np.linspace(1, 1, int(line_res / 2))))
    phi = np.concatenate((np.linspace(0, m.pi, int(line_res / 2)), np.linspace(m.pi, 0, int(line_res / 2))))

    # Initialize lists to hold the vstacked xyz coordinates and corresponding radii
    coords_list = []
    radii_list = []

    for i in range(rad_res):
        # Calculate the radius for this shell
        radius = np.ones_like(theta) * (r_outer - (i * (r_outer - r_inner) / (rad_res - 1)))
        
        # Convert spherical to Cartesian coordinates
        x, y, z = sphere2cart(radius, theta, phi)
        
        # Stack x, y, z as rows (axis=1) and append to the list
        coords_list.append(np.vstack((x, y, z)).T)  # Transpose to get columns x, y, z
        radii_list.append(radius)

    # Concatenate all coordinates along rows
    all_coords = np.vstack(coords_list)
    all_radii = np.concatenate(radii_list)  # Flatten the radii into a single array

    # Return the concatenated coordinates and corresponding radii as columns
    return all_coords, all_radii
    

def rotatex(points, delta):
    # Rotation matrix for rotation about the x-axis
    R_x = np.array([[1, 0, 0],
                    [0, np.cos(delta), -np.sin(delta)],
                    [0, np.sin(delta),  np.cos(delta)]])
    
    # Apply the rotation matrix to the points
    rot_points = np.dot(points, R_x.T)
    return rot_points

def rotatey(points, delta):
    # Rotation matrix for rotation about the y-axis
    R_y = np.array([[ np.cos(delta), 0, np.sin(delta)],
                    [ 0, 1, 0],
                    [-np.sin(delta), 0, np.cos(delta)]])
    
    # Apply the rotation matrix to the points
    rot_points = np.dot(points, R_y.T)
    return rot_points

def rotatez(points, delta):
    #rotation matrix for rotation about the z axis
    R_z = np.array([[np.cos(delta), -np.sin(delta), 0],
                [np.sin(delta),  np.cos(delta), 0],
                [0,              0,            1]])
    rot_points = np.dot(points, R_z.T)
    return rot_points

def twist(points, r, r_inner, r_outer):
    # determine which points that are within r_inner and r_outer
    mask = r <= r_outer
    phi = np.zeros_like(r)
    # cos function provides interpolation for twist angle based on r with minima at r_outer and maximum pi radians at r_inner
    # slope of graph is 0 at endpoints which profides smooth connection to the inside and outside
    phi[mask] = m.pi/2 * (1 - np.cos((r_outer - r[mask]) / (r_outer - r_inner) * m.pi))
    # extract xyz from points
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    xp = np.copy(x)  
    zp = np.copy(z)
    #twist points within r range by phi about y axis
    xp[mask] = x[mask] * np.cos(phi[mask]) + z[mask] * np.sin(phi[mask])
    zp[mask] = -x[mask] * np.sin(phi[mask]) + z[mask] * np.cos(phi[mask])
    # replace points
    npoints = np.vstack((xp, y, zp)).T
    return npoints

   
resolution = 5000
fibres_per_quadrant = 2
num_rings = 40
ang_res = 4 * fibres_per_quadrant
r_inner = 1
r_outer = 5
r_flat  = 1



# theta and phi must be 1d numpy arrays in pairs theta[i] phi[i]
theta1 = m.pi * np.concatenate((np.linspace(0,0,int(ang_res/2)+1), np.linspace(1, 1, int(ang_res/2)+1)))
phi1 = np.concatenate((np.linspace(0,m.pi, int(ang_res/2)+1), np.linspace(m.pi, 0, int(ang_res/2)+1)))
'''
theta2 = np.linspace(0,2*m.pi, ang_res+1)
phi2 = m.pi / 2 * np.ones(ang_res+1)

theta3 = m.pi * np.concatenate((np.linspace(1/2,1/2,int(ang_res/2)+1), np.linspace(3/2, 3/2, int(ang_res/2)+1)))
phi3 = np.concatenate((np.linspace(0,m.pi, int(ang_res/2)+1), np.linspace(m.pi, 0, int(ang_res/2)+1)))

theta = np.concatenate((theta1, theta2, theta3))
phi = np.concatenate((phi1, phi2, phi3))
'''
theta = theta1
phi = phi1


cart_points, r_points = points(theta, phi, r_inner, r_outer, r_flat, 1000)

grid_rings, r_rings = grid(r_outer, r_inner, num_rings, 1000)
#grid_rings_x = rotatez(grid_rings, m.pi/2)
#grid_rings_z = rotatex(grid_rings, m.pi/2)

#grid_rings = np.concatenate((grid_rings, grid_rings_x ,grid_rings_z))
#r_rings = np.concatenate((r_rings, r_rings, r_rings))


#animation control
num_frames = 800
rotation_proportion = 1/3
num_rotations = 3


for frame in range(num_frames):
    r = r_points

    rings_r = r_rings
    plotter = pv.Plotter(off_screen=True)
    sphere = pv.Sphere(r_inner)
    plotter.add_mesh(sphere, color = "black")
    #get a fresh copy of the points
    frame_points = rotatex(cart_points, (frame / num_frames) * m.pi)

    grid_points = rotatex(grid_rings, (frame / num_frames) * m.pi)
    print(frame)

    #here is the magic
    # get twist angle from frome
    delta = 2 * m.pi * frame / num_frames * num_rotations
    #twist ccw about z
    frame_points = rotatez(frame_points, delta)
    grid_points = rotatez(grid_points, delta)
    # apply spiral twist dependent on r about y axis
    frame_points = twist(frame_points, r, r_inner,r_outer)
    grid_points = twist(grid_points, rings_r, r_inner,r_outer)
    #twist back cw about z
    frame_points = rotatez(frame_points, -(1+rotation_proportion) * delta)
    grid_points = rotatez(grid_points, -(1+rotation_proportion) * delta)

    plotter.show_axes()
    mesh = pv.PolyData(frame_points)
    mesh_r = pv.PolyData(grid_points)
   
    norm_r = np.zeros_like(r)
    ring_norm = (np.max(rings_r) - rings_r) / (np.max(rings_r) - np.min(rings_r))
    #color based on rS
    mask = r <=r_outer
    norm_r[mask] = (np.max(r[mask])-r[mask]) / (np.max(r[mask]) - np.min(r[mask]))

    mesh['radius'] = norm_r
    mesh_r['radius'] = ring_norm
    plotter.add_mesh(mesh, scalars = 'radius', cmap = 'cool', show_scalar_bar=False, point_size = 11, render_points_as_spheres=True)
    plotter.add_mesh(mesh_r, scalars = 'radius',  cmap = 'hot', show_scalar_bar = False, point_size = 6, render_points_as_spheres = True)
    
    b = int((r_inner + r_outer + r_flat) * 2.2 )
    
    position = (b,b,b*3/4) 
    focal_point = (0, 0, 0) 
    view_up = (0, 0, 1)

    plotter.camera_position = [position, focal_point, view_up]
    
    #print(plotter.camera_position.to_list())
    #plotter.camera.elevation = -10
    plotter.screenshot(f"frame_{frame}")
    plotter.clear()
    plotter.close()




'''
#folder where images were savesd
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
'''
