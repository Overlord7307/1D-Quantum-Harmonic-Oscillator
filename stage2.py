'''
In stage 2, we take the foundation formed by stage 1 and rescale all the obtained values back to SI units. We also modify the physics engine so that it is now able to find energy eigenvalues and wavefunctions for any given state of the system.
'''

import numpy as np
import matplotlib.pyplot as plt

# Define physical constants in SI units
ℏ = 1.05457E-34 # m2kg/s
ω = 1E14 # rad/s
m = 1E-26 # kg
L = (ℏ / (m * ω))**0.5  # Characteristic length


# Define boundaries
x_start = -6.0
x_end = 6.0
dx = 0.01
x_values = [i for i in np.arange(x_start, x_end, dx)]   # Creates a list of all x-values


def V(x:float):
    '''
    Define the harmonic potential in natural units.
    '''
    return 0.5 * (x**2)


def ψ_ip1(x:float, E:float, ψ_i:float, ψ_im1:float):
    '''
    Solution to the Schrodinger equation using the finite difference method for the second derivative (natural units).
    '''
    return 2 * (V(x) - E) * ψ_i * dx**2 + 2*ψ_i - ψ_im1


def integrator(E_guess:float):
    '''
    Integrates the discretized solution to the Schrodinger equation using the guessed value of energy from x_start to x_end, taking increments of size dx. Finds the the wavefunction at every x and returns its value at x = x_end.
    '''
    ψ_values = [0, 0.0001]    # Sets initial values of ψ at x = x_start and x = x_start + dx
    for i in range(2, len(x_values)):
        ψ_new = ψ_ip1(x_values[i], E_guess, ψ_values[-1], ψ_values[-2])     # Finds next value of ψ using two previous values
        ψ_values.append(ψ_new)
    return ψ_values[-1]


def wavefunction_generator(E:float):
    '''
    Integrates the discretized solution to the Schrodinger equation using the ground state energy from x_start to x_end, taking increments of size dx. Returns a list of all values of ψ.
    '''
    ψ_values = [0, 0.0001]    # Sets initial values of ψ at x = x_start and x = x_start + dx
    for i in range(2, len(x_values)):
        ψ_new = ψ_ip1(x_values[i], E, ψ_values[-1], ψ_values[-2])   # Finds next value of ψ using two previous values
        ψ_values.append(ψ_new)
    return ψ_values


def count_nodes(E_guess:float):
    '''
    Counts the number of nodes formed by the wavefunction at a given energy level by counting the number of times it crosses the x-axis, i.e. the number of pairs of consecutive values having opposite signs.
    '''
    wave = wavefunction_generator(E_guess)  # Generates the wavefunction for given energy guess
    node_count = 0  # Initializes counter to count nodes
    for i in range(1, len(wave)):
        if (wave[i] * wave[i-1] < 0):
            node_count += 1     # Increments node_count of any two consecutive values have opposite signs
        else:
            pass
    return node_count


def n_scanner(n:int):
    '''
    Finds the energy bracket (El, Er), such that the wave with energy El has n nodes and the wave with energy Er has n+1 nodes. This bracket will be used in the bisection search to ensure that it hunts for exactly the nth energy state.
    '''
    E_values = [0]  # Starts checking energy from 0
    dE = 0.1    # Step value
    while True:
        if (count_nodes(E_values[-1]) <= n):
            E_values.append(E_values[-1] + dE)  # Adds the next energy value to the list
        else:
            break   # Breaks loop as soon as wave with n+1 nodes is encountered
    return (E_values[-2], E_values[-1])


def bisection_search(func, El:float, Er:float, max_error:float, n:int):
    '''
    Implements the bisection search algorithm to hunt for an energy value which lets the wavefunction vanish at the boundaries. Returns the approximate nth state energy.
    '''
    Em_values = [Er]    # Initialize a list to store guessed energy values
    error = 1   # Initialize the approximation error with its maximum value

    while (error >= max_error):
        Em = (El + Er) / 2
        Em_values.append(Em)
        divergence = func(Em) * ((-1)**n)
        if (divergence < 0):
            Er = Em
        elif (divergence > 0):
            El = Em
        else:
            return Em
        error = abs((Em_values[-1] - Em_values[-2]) / Em_values[-2])

    return Em


def main():
    # Set value of n (state of the system)
    n = 8

    # Create reduced lists of values to plot a limited range from x = -5 to x = +5
    x_axis = [(i * L) for i in np.arange(-5, 5, dx)]    # Scaling the x-axis back to SI units
    V_values = [(V(x) * m * ω**2) for x in x_axis]

    El, Er = n_scanner(n)
    Energy_n = bisection_search(integrator, El, Er, 0.0000001, n)
    Energy_n_SI = Energy_n * ℏ * ω    # Converting energy to SI units
    wavefunction = wavefunction_generator(Energy_n)

    # Calculating the normalization constant
    area_under_ψ = 0
    for val in wavefunction:
        area_under_ψ += (val**2) * dx
    N = area_under_ψ ** 0.5

    normalized_wavefunction_SI = [(val/N) / (L**0.5) for val in wavefunction]   # Normalizing and scaling wavefunction to SI units
    chopped_wavefunction = normalized_wavefunction_SI[100:-100]    # Chopping off the ends for a better plot

    # Creating a map to dynamically print subscript n
    subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    n_sub = str(n).translate(subscript_map)


    print(f"The energy of the {n}{'th' if 11 <= n % 100 <= 13 else ('st' if n % 10 == 1 else ('nd' if n % 10 == 2 else ('rd' if n % 10 == 3 else 'th')))} state is: {Energy_n_SI} J")
    print(f'Theoretical energy value using formula: {(n + 1/2) * ℏ * ω} J')
    print()


    # Plotting the harmonic potential, the nth state wavefunction and the nth state energy against the x-values

    fig, ax1 = plt.subplots(figsize=(16, 9))
    plt.title(f"{n}{'th' if 11 <= n % 100 <= 13 else ('st' if n % 10 == 1 else ('nd' if n % 10 == 2 else ('rd' if n % 10 == 3 else 'th')))} State Wavefunction for a Quantum Harmonic Oscillator")

    # Plotting potential and nth state energy
    ax1.plot(x_axis, V_values, label='Harmonic Potential V(x)', color='black', linestyle='--', linewidth=1)
    ax1.axhline(Energy_n_SI, label=f"{n}{'th' if 11 <= n % 100 <= 13 else ('st' if n % 10 == 1 else ('nd' if n % 10 == 2 else ('rd' if n % 10 == 3 else 'th')))} State Energy E{n_sub} = {Energy_n_SI:.4e} J", color='red')
    ax1.set_xlabel('Position (metres)')
    ax1.set_ylabel('Energy (Joules)')
    ax1.axhline(0, color='black', linewidth=0.75)
    ax1.axvline(0, color='black', linewidth=0.75)

    # Plotting the wavefunction
    ax2 = ax1.twinx()
    ax2.plot(x_axis, chopped_wavefunction, label=f"{n}{'th' if 11 <= n % 100 <= 13 else ('st' if n % 10 == 1 else ('nd' if n % 10 == 2 else ('rd' if n % 10 == 3 else 'th')))} State Wavefunction \u03A8{n_sub}(x)", color='blue')
    ax2.set_ylabel('Wavefunction Amplitude (1/\u221Am)')
    ax2.axhline(0, color='black', linewidth=0.75)

    ax1.grid(True)
    fig.legend()
    plt.savefig(f'stage2_plot_n={n}.png', dpi=400, bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    main()