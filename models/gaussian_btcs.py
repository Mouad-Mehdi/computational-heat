import numpy as np
import utils.utilities as ut

# Physical parameters :
k = 429          # Thermal conductivity in W/(m*K)
p = 10490        # Density in kg/m^3
c = 235          # Specific heat capacity in J/(kg*K)
alpha = k/(p*c)  # Thermal diffusivity in m^2/s

# Spatial domain setup :
L = 1     # Length of the rod in m
I = 400   # Number of spatial points
dx = L/I  # Spatial step

# Temporal domain setup :
H = 110                   # Simulation time in s
target_r = 0.4            # Target r
dt = target_r*dx**2/alpha # Time step
N = round(H / dt)         # Number of temporal points

# Final calculated r
r = alpha * dt / (dx**2)

# Gaussian initial temperature distribution:
mean = 0.5
std_dev = 0.05

def f(x):
    return 1000*(1/(std_dev*np.sqrt(2*np.pi)))*np.exp(-((x-mean)**2)/(2*std_dev**2))


# Setting up the space time grid:
x = np.linspace(0,L,I+1)
t = np.linspace(0,H,N+1)


# Initialising the state vector:
T_initial = f(x)

# Numerically solving the problem:
T_numerical = ut.heat_btcs(I, N, r, T_initial)

# Plotting a heatmap of the numerical solution vs the analytical one:
ut.heatmap(T_numerical,L ,H, "Gaussian initial temperature distribution")