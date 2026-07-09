import numpy as np
import matplotlib.pyplot as plt

# An implementation of the FTCS approximation:
def heat_ftcs(J, N, r, T_initial):
    # Checking the stability conditions:
    if r > 0.5 :
        raise ValueError(f"Stability condition violated! r must be <= 0.5. Current r = {r:.4f}")
       
    # Initializing the space time grid matrix:
    T = np.zeros((J+1,N+1))

    # Applying the initial conditions:
    T[:,0] = T_initial

    # Forcing the boundary conditions:
    T[0,:] = 0
    T[J,:] = 0

    # FTCS loop:
    for n in range(0,N):
        T[1:J, n+1] = (1 - 2*r) * T[1:J, n] + r * (T[2:J+1, n] + T[0:J-1, n])

    return T

# An implementation of Thomas' algorithm:
def thomas(A,d):
    # Setting up the parameters:
	upper = np.diagonal(A, offset = 1).copy()
	main = np.diagonal(A).copy()
	lower = np.diagonal(A, offset = -1).copy()
	d = d.copy()
	n = len(main)
	x = np.zeros(n)
	
    # Forward elimination:
	for i in range(1,n):
		pivot = lower[i-1] / main[i-1]
		main[i] = main[i] - pivot * upper[i-1]
		d[i] = d[i] - pivot * d[i-1]
	
    # Back substitution
	x[n-1] = d[n-1] / main[n-1]
	for i in range(n-2,-1,-1):
		x[i] = (d[i] - upper[i] * x[i+1]) / main[i]
		
	return x

# An implementation of the BTCS approximation:
def heat_btcs(J, N, r, T_initial):

    # Initializing the linear system matrix:
    A = np.diag(np.full(J-1, 1 + 2*r)) + np.diag(np.full(J-2, -r), k=1) + np.diag(np.full(J-2, -r), k=-1)

    # Initializing the space-time matrix:
    T = np.zeros((J+1,N+1))

    # Applying the initial conditions:
    T[:,0] = T_initial

    # Forcing the boundary conditions:
    T[0,:] = 0
    T[J,:] = 0

    # BTCS loop:
    for n in range(1,N+1):
        T[1:J,n] = thomas(A,T[1:J,n-1])

    return T

# Plotting a heatmap:
def heatmap(data, L, H, title):
    plt.figure()
    heat = plt.imshow(data, cmap='inferno', aspect="auto", origin='lower',extent=[0,H,0,L])
    cbar = plt.colorbar(heat)
    cbar.set_label("Temperature (T)", size=11)
    plt.xlabel("Time (s)", fontsize=11)
    plt.ylabel("Position along the Rod (x)", fontsize=11)
    plt.title(title)
    plt.tight_layout()
    plt.show()

# Comparing the analytical and numerical solutions via two heatmaps:
def heatmap_comp(analytical, numerical, L, H, v_max):
    # Calculate the absolute difference matrix:
    error = np.abs(numerical - analytical)

    # Create 1 row, 3 columns
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

    limits = [0, H, 0, L]

    # Plot Numerical
    im_num = axes[0].imshow(numerical, cmap='inferno', aspect='auto', 
                        origin='lower', vmin=0, vmax=v_max, extent=limits)
    axes[0].set_title("Numerical Solution")
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Position along the Rod (m)")

    # Plot Analytical
    axes[1].imshow(analytical, cmap='inferno', aspect='auto', 
                origin='lower', vmin=0, vmax=v_max, extent=limits)
    axes[1].set_title("Analytical Solution")
    axes[1].set_xlabel("Time (s)")

    # Creating the temperature color bar:
    cbar_temp = fig.colorbar(im_num, ax=[axes[0], axes[1]])
    cbar_temp.set_label("Temperature (T)", size=12)

    # Absolute Error Heatmap 
    im_err = axes[2].imshow(error, cmap='viridis', aspect='auto', origin='lower',
                            extent=limits, vmin=0, vmax=error.max())
    axes[2].set_title("Absolute Error |Num - Ana|")
    axes[2].set_xlabel("Time (s)")

    # Add a colorbar for the error heatmap
    cbar_err = fig.colorbar(im_err, ax=axes[2])
    cbar_err.set_label("Temperature difference in C", size=11)

    plt.show()

