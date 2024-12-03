#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

# Parameters
R = 8.1  # Radius in meters
g = 9.81  # Acceleration due to gravity in m/s^2
omega_theta_max = 1.414  # Max yaw angular velocity in rad/s
omega_phi = np.pi / 40  # Pitch angular velocity in rad/s
T_total = 100  # Total simulation time in seconds
t = np.linspace(0, T_total, 8000)  # Time array
l2 = 0.7

def alpha(t):
    return np.piecewise(
        t,
        [t <= 20, 
         (t > 20) & (t <= 40), 
         (t > 40) & (t <= 60), 
         (t > 60) & (t <= 80),
         t > 80],
        [
            lambda t: 0, 
            lambda t: (np.pi / 2) * ((t - 20) / (40 - 20))**2,  
            lambda t: np.pi / 2,
            lambda t: (np.pi / 2) * ((80 - t) / (80 - 60))**2,  
            lambda t: 0
        ]
    )


def omega_theta_t(t):
    return np.piecewise(
        t,
        [t <= 20,
         (t > 20) & (t <= 80),
         (t > 80) & (t <= 100)],
        [lambda t: omega_theta_max * (t / 20),
         lambda t: omega_theta_max,
         lambda t: omega_theta_max * (1 - (t - 80) / 20)]
    )

def omega_phi_t(t):
    return np.piecewise(
        t,
        [(t > 20) & (t <= 40),
         (t > 60) & (t <= 80)],
        [omega_phi, -omega_phi],
        0
    )

def phi(t):
    return np.piecewise(
        t,
        [t <= 20,
         (t > 20) & (t <= 40),
         (t > 40) & (t <= 60),
         (t > 60) & (t <= 80)],
        [0,
         lambda t: omega_phi * (t - 20),
         np.pi / 2,
         lambda t: np.pi / 2 - omega_phi * (t - 60)]
    )

def theta_array(t, omega_theta_values):
    return cumtrapz(omega_theta_values, t, initial=0) % (2 * np.pi)

def calculate_Reff(R, alpha_values):
    return np.sqrt(R**2 + l2**2 + 2 * R * l2 * np.sin(alpha_values))

def Gx_formula(phi_values, omega_phi_values, omega_theta_values, alpha_values, theta_values, R):
    Reff = calculate_Reff(R, alpha_values)
    term1 = 0  
    gx_val = g * np.sin(phi_values) * np.cos(theta_values)
    return (term1 - gx_val) / g


def Gy_formula(phi_values, omega_phi_values, omega_theta_values, alpha_values, theta_values, R):
    Reff = calculate_Reff(R, alpha_values)
    term1 = -Reff * omega_theta_values**2 * np.cos(alpha_values) + Reff * omega_phi_values**2 * np.cos(alpha_values)
    gy_val = g * np.sin(theta_values) * np.cos(alpha_values) * np.sin(phi_values) + g * np.cos(phi_values) * np.sin(alpha_values)
    return (term1 - gy_val) / g

def Gz_formula(phi_values, omega_phi_values, omega_theta_values, alpha_values, theta_values, R):
    Reff = calculate_Reff(R, alpha_values)
    term1 = Reff * omega_theta_values**2 * np.sin(alpha_values) + Reff * omega_phi_values**2 * np.sin(alpha_values)
    gz_val = g * np.sin(theta_values) * np.sin(alpha_values) * np.sin(phi_values) - g * np.cos(phi_values) * np.cos(alpha_values)
    return (term1 + gz_val) / g




omega_theta_values = omega_theta_t(t)
theta_values = theta_array(t, omega_theta_values)
alpha_values = alpha(t)
phi_values = phi(t)
omega_phi_values = omega_phi_t(t)

Gx_values = Gx_formula(phi_values, omega_phi_values, omega_theta_values, alpha_values, theta_values, R)
Gy_values = Gy_formula(phi_values, omega_phi_values, omega_theta_values, alpha_values, theta_values, R)
Gz_values = Gz_formula(phi_values, omega_phi_values, omega_theta_values, alpha_values, theta_values, R)

# Print Results
print(f"Gx max: {Gx_values.max():.3f}, Gx min: {Gx_values.min():.3f}")
print(f"Gy max: {Gy_values.max():.3f}, Gy min: {Gy_values.min():.3f}")
print(f"Gz max: {Gz_values.max():.3f}, Gz min: {Gz_values.min():.3f}")

# Plot results
plt.figure(figsize=(10, 3))
plt.plot(t, Gx_values, label="Gx", color="#1f77b4")
plt.xlabel("Time (s)")
plt.ylabel("Gx(g)")
plt.grid()
plt.legend()
plt.show()

plt.figure(figsize=(10, 3))
plt.plot(t, Gy_values, label="Gy", color="green")
plt.xlabel("Time (s)")
plt.ylabel("Gy(g)")
plt.grid()
plt.legend()
plt.show()

plt.figure(figsize=(10, 3))
plt.plot(t, Gz_values, label="Gz", color="red")
plt.xlabel("Time (s)")
plt.ylabel("Gz(g)")
plt.grid()
plt.legend()
plt.show()


fig, ax = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
ax[0].plot(t, Gx_values, label="Gx", color="#1f77b4")
ax[0].set_ylabel("Gx (g)")
ax[0].grid()
ax[0].legend()
ax[1].plot(t, Gy_values, label="Gy", color="green")
ax[1].set_ylabel("Gy (g)")
ax[1].grid()
ax[1].legend()
ax[2].plot(t, Gz_values, label="Gz", color="red")
ax[2].set_xlabel("Time (s)")
ax[2].set_ylabel("Gz (g)")
ax[2].grid()
ax[2].legend()

plt.tight_layout()
plt.show()
