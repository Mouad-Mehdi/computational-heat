import numpy as np
import utils.utilities as ut

# Physical parameters :
k = 429          # Thermal conductivity in W/(m*K)
p = 10490        # Density in kg/m^3
c = 235          # Specific heat capacity in J/(kg*K)
alpha = k/(p*c)  # Thermal diffusivity in m^2/s

# Spatial domain setup :
L = 1     # Length of the rod in m
J = 100   # Number of spatial points
dx = L/J  # Spatial step

# Temporal domain setup :
H = 850                   # Simulation time in s
target_r = 0.4            # Target r
dt = target_r*dx**2/alpha # Time step
N = round(H / dt)         # Number of temporal points

# Final calculated r
r = alpha * dt / (dx**2)

#The analytical solution of the problem:
def Analytical(x,t):
    return 1000*np.sin(np.pi * x / L)*np.exp(-((np.pi / L)**2 )* alpha * t)

# Setting up the space time grid:
x = np.linspace(0,L,J+1)
t = np.linspace(0,H,N+1)

T_analytical = Analytical(x[:, None], t)

# Initialising the state vector:
T_initial = Analytical(x,0)

# Numerically solving the problem:
T_numerical = ut.heat_ftcs(J, N, r, T_initial)

# Plotting a heatmap of the numerical solution vs the analytical one:
ut.heatmap_comp(T_analytical, T_numerical,L, H, 1000)