'''
We modify the stage 3 code to solve for an anharmonic potential.
'''

import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

# Define physical constants in SI units
ℏ = 1.05457E-34 # m2kg/s
ω = 1E14 # rad/s
m = 1E-26 # kg
L = (ℏ / (m * ω))**0.5  # Characteristic length
alpha = 0.75   # Parameter for anharmonic potential

# Define boundaries
x_start = -6.0
x_end = 6.0
x_values = np.linspace(x_start, x_end, 2000)   # Creates an array of all x-values
dx = x_values[1] - x_values[0]  # Define step size


def V(x:float):
    '''
    Define the anharmonic potential in natural units.
    '''
    return 0.5 * (x**2) + alpha * x**4


def main():
    # Set value of n (state of the system)
    n = 8

    # Creating a map to dynamically print subscript n
    subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    n_sub = str(n).translate(subscript_map)

    # Construct potential energy matrix
    V_matrix = np.diag(V(x_values), k=0)

    # Construct kinetic energy matrix
    upper = np.full((1999), 1)
    middle = np.full((2000), -2)
    lower = np.full((1999), 1)
    K_matrix = -(np.diag(middle, k=0) + np.diag(upper, k=1) + np.diag(lower, k=-1)) / (2 * dx**2)   # Add 2 off-diagonal matrices to the diagonal matrix

    # Construct the Hamiltonian of the system
    Hamiltonian = K_matrix + V_matrix
    energies, wavefunctions = la.eigh(Hamiltonian)  # Use SciPy to calculate eigenvalues and eigenfunctions

    Energy_n_SI = energies[n] * ℏ * ω   # Extract nth state energy and convert to SI units
    print(f"The energy of the {n}{'th' if 11 <= n % 100 <= 13 else ('st' if n % 10 == 1 else ('nd' if n % 10 == 2 else ('rd' if n % 10 == 3 else 'th')))} state is: {Energy_n_SI} J")
    print()


    x_axis_SI = x_values * L    # Scale x-axis to SI
    V_values_SI = V(x_axis_SI) * m * ω**2   # Scale harmonic potential to SI
    normalized_wave_n_SI = (wavefunctions[:, n] / dx**0.5) / L**0.5     # Extract nth state wavefunction, normalize, and scale to SI

    # Plotting the anharmonic potential, the nth state wavefunction and the nth state energy against the x-values
    fig, ax1 = plt.subplots(figsize=(16, 9))
    plt.title(f"{n}{'th' if 11 <= n % 100 <= 13 else ('st' if n % 10 == 1 else ('nd' if n % 10 == 2 else ('rd' if n % 10 == 3 else 'th')))} State Wavefunction for a 1D Quantum Anharmonic Oscillator")

    # Plotting potential and nth state energy
    ax1.plot(x_axis_SI, V_values_SI, label='Anharmonic Potential V(x)', color='black', linestyle='--', linewidth=1)
    ax1.axhline(Energy_n_SI, label=f"{n}{'th' if 11 <= n % 100 <= 13 else ('st' if n % 10 == 1 else ('nd' if n % 10 == 2 else ('rd' if n % 10 == 3 else 'th')))} State Energy E{n_sub} = {Energy_n_SI:.4e} J", color='red')
    ax1.set_xlabel('Position (metres)')
    ax1.set_ylabel('Energy (Joules)')
    ax1.axhline(0, color='black', linewidth=0.75)
    ax1.axvline(0, color='black', linewidth=0.75)

    # Plotting the wavefunction
    ax2 = ax1.twinx()
    ax2.plot(x_axis_SI, normalized_wave_n_SI, label=f"{n}{'th' if 11 <= n % 100 <= 13 else ('st' if n % 10 == 1 else ('nd' if n % 10 == 2 else ('rd' if n % 10 == 3 else 'th')))} State Wavefunction \u03A8{n_sub}(x)", color='blue')
    ax2.set_ylabel('Wavefunction Amplitude (1/\u221Am)')
    ax2.axhline(0, color='black', linewidth=0.75)

    ax1.grid(True)
    fig.legend()
    plt.savefig(f'plots\\anharmonic_plots\\anharmonic_plot_n={n}.png', dpi=400, bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    main()