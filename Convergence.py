import numpy as np
import utils.utilities as ut
import matplotlib.pyplot as plt

# Solver type: 'FTCS', 'BTCS' or 'CN' 
solver = 'BTCS' 

# Study type: 'time'or 'space' 
study_type = 'time'

# Physical parameters :
k = 429          # Thermal conductivity in W/(m*K)
p = 10490        # Density in kg/m^3
c = 235          # Specific heat capacity in J/(kg*K)
alpha = k/(p*c)  # Thermal diffusivity in m^2/s

# Spatial domain setup :
L = 1     # Length of the rod in m

# Temporal domain setup :
H = 100   # Simulation time in s

# The analytical solution of the problem:
def Analytical(x, t):
    return (
        1000 * np.sin(np.pi * x / L)
        * np.exp(-(np.pi / L)**2 * alpha * t)
        + 200 * np.sin(3 * np.pi * x / L)
        * np.exp(-9 * (np.pi / L)**2 * alpha * t)
    )

# Wrapper to route the convergence study to the correct solver:
def run_solver(J, N, r, T_initial, solver_type):
    if solver_type == 'CN':
        return ut.heat_cn(J, N, r, T_initial)
    if solver_type == 'FTCS':
        return ut.heat_ftcs(J, N, r, T_initial)
    if solver_type == 'BTCS':
        return ut.heat_btcs(J, N, r, T_initial)
    else:
        raise ValueError("Unknown solver")
    
def temporal_convergence(solver_name):
    print(f"Running the temporal convergence study for {solver_name}")

    # Spatial resolution, kept high to isolate temporal error:
    J = 100 if solver_name == "FTCS" else 500
    dx = L/J
    x = np.linspace(0, L, J+1)
    T_initial = Analytical(x, 0)

    # dt values:
    dt_values = [] if solver_name =="FTCS" else [0.1, 0.05, 0.025, 0.0125]
    errors = []

    for dt in dt_values:
        N = round(H / dt)
        actual_dt = H / N
        r = alpha * actual_dt / (dx**2)

        if solver_name == 'FTCS' and r > 0.5:
            raise ValueError(f"Warning, the FTCS stability condition has been violated, current r = {r:.2f}")
        
        t = np.linspace(0, H, N+1)
        T_analytical = Analytical(x[:, None], t)

        T_numerical = run_solver(J, N, r, T_initial, solver_name)

        error = np.linalg.norm(T_numerical[:, -1] - T_analytical[:, -1], ord=np.inf)

        errors.append(error)

    D = np.array(dt_values)
    errors = np.array(errors)

    fit = np.polyfit(np.log(D), np.log(errors), 1)
    order = fit[0]
    print(f"Observed temporal order : {order:.3f}")

    # CN is second order accurate in time, FTCS and BTCS are first order:
    expected_slope = 2 if solver_name == 'CN' else 1
    reference = (errors[-1] / D[-1]**expected_slope) * D**expected_slope

    # Plotting the convergence study:

    plt.figure(figsize=(8, 6))
    plt.loglog(D, errors, 'o-', label=f"{solver_name} Empirical Error")
    plt.loglog(D, reference, '--', label=f"Theoretical O(dt^{expected_slope})")
    plt.title(f"Temporal Convergence ({solver_name})")
    plt.xlabel("dt (s)")
    plt.ylabel("Max Absolute Error")
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.savefig("images/btcs_temporal_convergence.png", dpi =400)
    plt.show()

def spatial_convergence(solver_name):
    print(f"Running the spatial convergence study for {solver_name}")

    J_values = [10, 20, 40, 80, 160]
    dx_values = []
    errors = []

    for J in J_values:
        dx = L/J
        dx_values.append(dx)

        # Maintain a constant r to isolate the spatial error
        r_target = 0.4
        dt_target = r_target * (dx**2) / alpha

        N = round(H / dt_target)
        actual_dt = H / N
        r = actual_dt * alpha /(dx**2)

        x = np.linspace(0, L, J+1)
        t = np.linspace(0, H, N+1)
        
        T_initial = Analytical(x, 0)
        T_analytical = Analytical(x[:, None], t)

        T_numerical = run_solver(J, N, r, T_initial, solver_name)
            
        error = np.linalg.norm(T_numerical[:, -1] - T_analytical[:, -1], ord=np.inf)
        errors.append(error)

    D = np.array(dx_values)
    errors = np.array(errors)

    # Calculate observed order
    fit = np.polyfit(np.log(D), np.log(errors), 1)
    order = fit[0]
    print(f"Observed Spatial Order: {order:.3f}")
    
    # All three solvers are second order accurate in space
    expected_slope = 2 
    reference = (errors[-1] / D[-1]**expected_slope) * D**expected_slope
    
    plt.figure(figsize=(8, 6))
    plt.loglog(D, errors, 'o-', label=f"{solver_name} Empirical Error")
    plt.loglog(D, reference, '--', label=f"Theoretical O(dx^{expected_slope})")
    plt.title(f"Spatial Convergence ({solver_name})")
    plt.xlabel("dx (Meters)")
    plt.ylabel("Max Absolute Error")
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.show()

# Running the convergence study:
if __name__ == '__main__':
    if study_type == 'time':
        temporal_convergence(solver)
    if study_type == 'space':
        spatial_convergence(solver)