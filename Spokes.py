import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# Constants
m1 = 17725           # Mass for spokes including capsule weight (kg)   
m1error = 0.1        # Relative error in m1
m2 = 26525           # Mass for spokes excluding capsule weight (kg)
m2error = 0.1        # Relative error in m2
omega = 0.00349      # Angular velocity (rad/s)
g = 9.81             # Gravitational acceleration (m/s^2)
n_spokes = 80        # Number of spokes
r = 60               # Radius of the wheel (m)
rerror = 0.05        # Relative error in r
A = np.pi * 0.035**2 # Cross-sectional area of spokes (m^2)
p = 1.225            # Air density (kg/m^3)
perror = 0.05        # Relative error in p

T0 = m2 * (g - omega ** 2 * r)  # Min Pre-tension required in spokes to gaurantee spokes are always under tension, not compression (N)
print(f"T0 = {T0}")
pre_tension = T0 + 26021.02429

# Function to find error in tension based off errors in m and r
def spokes_error(theta, m, omega, r, g):
    m_error = m * 0.1
    r_error = r * rerror

    # Partial derivatives
    dT_dm = omega**2 * r - g * np.sin(theta)
    dT_dr = m * omega**2

    # Combine as per error propgation equation
    return np.sqrt((dT_dm * m_error)**2 + (dT_dr * r_error)**2)

max_error = spokes_error(np.pi / -2, m2, omega, r, g)
print(f"max_error = {max_error}")

pre_tension = T0 + max_error # Adding the max error onto the pre-tension guarantees tension even in all cases including errors.

# Function to find tension
def tension(theta, m, omega, r, g, pre_tension):
    return pre_tension + m * omega**2 * r - m * g * np.sin(theta)

# Create a set of 80 angles to plot individual spokes with
theta_spokes = np.linspace(-np.pi/2, 3 * np.pi/2, n_spokes, endpoint=False)

# Create a large set of angles to plot 'continuous' curves with
theta_continuous = np.linspace(-np.pi / 2, 3 * np.pi / 2, 2000)

# Creating arrays to store tensions and errors for each spoke
tensions_spokes = np.zeros(n_spokes)
errors_spokes = np.zeros(n_spokes)
colors = []

# Loop through each spoke and calculate tension and error with alternating masses
for i in range(n_spokes):
    if i % 2 == 0:  # Even spokes use m1
        tensions_spokes[i] = tension(theta_spokes[i], m1, omega, r, g, pre_tension)
        errors_spokes[i] = spokes_error(theta_spokes[i], m1, omega, r, g)
        colors.append('blue')
    else:  # Odd spokes use m2
        tensions_spokes[i] = tension(theta_spokes[i], m2, omega, r, g, pre_tension)
        errors_spokes[i] = spokes_error(theta_spokes[i], m2, omega, r, g)
        colors.append('orange')

print(errors_spokes)

# Creating plot
plt.figure(figsize=(10, 6))

# Plotting separate curves for m1 and m2
plt.plot(theta_continuous, tension(theta_continuous, m1, omega, r, g, pre_tension), label='Spokes bearing capsule weight', color='blue')
plt.plot(theta_continuous, tension(theta_continuous, m2, omega, r, g, pre_tension), label='Spokes not bearing capsule weight', color='orange')

# Scatter plot with error bars
for i in range(n_spokes):
    plt.errorbar(theta_spokes[i], tensions_spokes[i], yerr=errors_spokes[i], fmt='o', color=colors[i], zorder=5, label=None)

# Axes labels and title
plt.xlabel('Yaw Angle $\\Theta$ (radians)')
plt.ylabel('Tension (N)')

# Pre-tension line
plt.axhline(pre_tension, color='red', linestyle='--', linewidth=0.8, label=f'Minimum Pre-Tension ({pre_tension} N)')

# Custom legend
blue_line = mlines.Line2D([], [], marker='o', color='blue', linestyle='None', markersize=10, label='Spokes bearing capsule weight')
orange_line = mlines.Line2D([], [], marker='o', color='orange', linestyle='None', markersize=10, label='Spokes not bearing capsule weight')
pre_tensed_line = mlines.Line2D([], [], color='red', linestyle='-', linewidth=0.8, label='Minimum required pre-tension (N)')
plt.legend(handles=[blue_line, orange_line, pre_tensed_line], loc='upper center')

# Grid and show
plt.grid(True)
plt.show()
