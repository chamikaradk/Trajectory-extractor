import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCReader
import numpy as np

topology = 'F4_bulk.gro'
trajectory = 'F4_bulk.xtc'
interval_ps = 20 # time interval has to be in ps

print(f"Time interval selected : {interval_ps} ps")
def write_to_file(universe, frame, filename, supercell):
    atoms = universe.atoms
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"Frame {frame} coords box_coordinates {supercell} \n")
        for atom in atoms:
            pos = atom.position
            f.write(f"{atom.name} {pos[0]:.3f} {pos[1]:.3f} {pos[2]:.3f} \n") 

def extract_coords(topology, trajectory, interval_ps):
    universe = mda.Universe(topology,trajectory)
    dt = universe.trajectory.dt
    interval_frames = int(interval_ps/dt)
    box_coords = []
    for ts in universe.trajectory[::interval_frames]:
        box_coords.append(ts.dimensions)
        frame, time = ts.frame, ts.time
        supercell = ts.dimensions.copy()
        supercell[:3]=supercell[:3]/10
        filename = f"frame_{frame}.xyz"
        write_to_file(universe, frame, filename, supercell)
        print(f"File written: {filename}")
    box_coords = np.array(box_coords)
    return box_coords

#box coordinates are stored in this array
box_coordinates = extract_coords(topology, trajectory, interval_ps)
box_coordinates[:,:3] = box_coordinates[:,:3]/10 #obtain a,b,c in nm units

print(box_coordinates)
