import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.lines as mlines


m = 2000  # Mass assigned to each spoke 
omega = 1.414  # Angular velocity 
r = 8.1  # Radius of the wheel 
g = 9.81  # Gravitational Acceleration 
n_spokes = 20  # Number of spokes
pre_tension = 1000 # Pre-tension (N)


total_time = 100  # seconds
time_steps = 500  # number of time steps for simulation
time = np.linspace(0, total_time, time_steps)

# Yaw angles for each spoke 
theta_spokes = np.linspace(-np.pi / 2, 3 * np.pi / 2, n_spokes, endpoint=False)

# Time-dependent angular velocity 
def angular_velocity(t):
    if t < 20:
        return (omega / 20) * t
    elif t < 80:
        return omega
    elif t <= 100:
        return (omega / 20) * (100 - t)
    else:
        return 0

# Time-dependent tilt angle 
def tilt_angle(t):
    if t < 20:
        return 0
    elif t < 40:
        return ((np.pi/2) / 20) * (t - 20)  
    elif t < 60:
        return (np.pi/2)  
    elif t < 80:
        return ((np.pi/2) / 20) * (80 - t)  
    else:
        return 0

# Tension equation
def tension(theta, phi, omega, m, r, g, pre_tension):
    centripetal_force = m * omega**2 * r   
    gravitational_force = -m * g * np.sin(phi) * np.sin(theta)  
    return pre_tension + centripetal_force + gravitational_force

# Calculate tensions at each time step
tensions = []
for t in time:
    phi_t = tilt_angle(t)
    omega_t = angular_velocity(t)
    tension_values = tension(theta_spokes, phi_t, omega_t, m, r, g, pre_tension)
    tensions.append(tension_values)

# Convert tensions to array
tensions = np.array(tensions)

cmap = get_cmap('viridis')  
colors = [cmap(i / (n_spokes - 1)) for i in range(n_spokes)]

# Key spokes for legend
key_spokes = {
    "Bottom ($\\theta$ = -π/2)": -np.pi / 2,
    "Right ($\\theta$ = 0)": 0,
    "Top ($\\theta$ = π/2)": np.pi / 2,
    "Left ($\\theta$ = π)": np.pi
}

# Create the plot
plt.figure(figsize=(9, 6))
for i, color in enumerate(colors):
    if theta_spokes[i] in key_spokes.values():
        # Add key spokes to the legend
        label = [k for k, v in key_spokes.items() if np.isclose(v, theta_spokes[i])][0]
    else:
        label = None
    plt.plot(time, tensions[:, i], color=color, label=label)

# Add labels and grid
plt.xlabel("Time (s)", fontsize=16)
plt.ylabel("Tension (N)", fontsize=16)
plt.grid()
plt.legend(loc="upper left", bbox_to_anchor=(0.75, 1), fontsize="large", title="Key Spoke Angles")
plt.tight_layout()
plt.show()

# Create a 3D plot
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111, projection='3d')
theta_grid, time_grid = np.meshgrid(theta_spokes, time)
surf = ax.plot_surface(time_grid, theta_grid, tensions, cmap='viridis', edgecolor='none')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Yaw Angle $\\theta$ (radians)')
ax.set_zlabel('Tension (N)')
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label='Tension (N)')
plt.tight_layout()
plt.show()

