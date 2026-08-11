'''
In stage 2, we take the foundation formed by stage 1 and rescale all the obtained values back to SI units
'''

import numpy as np
import matplotlib.pyplot as plt

ℏ = 1.05457E-34 # m2kg/s
ω = 1E14 # rad/s
m = 1E-26 # kg
L = (ℏ / (m * ω))**0.5  # Characteristic length

x_start = -5.0
x_end = 5.0
dx = 0.01
x_values = [i for i in np.arange(x_start, x_end, dx)]   # Creates a list of all x-values

def V(x:float):
    '''
    Define the harmonic potential.
    '''
    return 0.5 * (x**2)


def ψ_ip1(x:float, E:float, ψ_i:float, ψ_im1:float):
    '''
    Solution to the Schrodinger equation using the finite difference method for the second derivative.
    '''
    return 2 * (V(x) - E) * ψ_i * dx**2 + 2*ψ_i - ψ_im1


def integrator(E_guess:float):
    '''
    Integrates the discretized solution to the Schrodinger equation using the guessed value of energy from x_start to x_end, taking increments of size dx. Returns the value of the wavefunction at x = x_end.
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


def bisection_search(func, El:float, Er:float, max_error:float):
    '''
    Implements the bisection search algorithm to hunt for an energy value which lets the wavefunction vanish at the boundaries. Returns the approximate ground state energy.
    '''
    Em_values = [Er]    # Initialize a list to store guessed energy values
    error = 1   # Initialize the approximation error with its maximum value

    while (error >= max_error):
        Em = (El + Er) / 2
        Em_values.append(Em)
        if (func(Em) < 0):
            Er = Em
        elif (func(Em) > 0):
            El = Em
        else:
            return Em
        error = abs((Em_values[-1] - Em_values[-2]) / Em_values[-2])

    return Em


# Create reduced lists of values to plot a limited range from x = -4 to x = +4
x_axis = [(i * L) for i in np.arange(-4.0, 4.0, dx)]    # Scaling the x-axis back to SI units
V_values = [(V(x) * m * ω**2) for x in x_axis]

Energy_Ground = bisection_search(integrator, 0, 1, 0.0000001)
Energy_Ground_SI = Energy_Ground * ℏ * ω    # Converting energy to SI units
wavefunction = wavefunction_generator(Energy_Ground)

# Calculating the normalization constant
area_under_ψ = 0
for val in wavefunction:
    area_under_ψ += (val**2) * dx
N = area_under_ψ ** 0.5

normalized_wavefunction_SI = [(val/N) / (L**0.5) for val in wavefunction]   # Normalizing and scaling wavefunction to SI units
chopped_wavefunction = normalized_wavefunction_SI[100:-100]    # Chopping off the ends for a better plot

print('The energy of the ground state is:', Energy_Ground_SI, 'J')
print()

# Plotting the harmonic potential, the ground state wavefunction and the ground state energy against the x-values

fig, ax1 = plt.subplots(figsize=(16, 9))
plt.title('Ground State Wavefunction for a Quantum Harmonic Oscillator')

# Plotting potential and ground state energy
ax1.plot(x_axis, V_values, label='Harmonic Potential V(x)', color='black', linestyle='--', linewidth=1)
ax1.axhline(Energy_Ground_SI, label=f'Ground State Energy E\u2080 = {Energy_Ground_SI:.4e} J', color='red')
ax1.set_xlabel('Position (metres)')
ax1.set_ylabel('Energy (Joules)')
ax1.axhline(0, color='black', linewidth=0.75)
ax1.axvline(0, color='black', linewidth=0.75)

# Plotting the wavefunction
ax2 = ax1.twinx()
ax2.plot(x_axis, chopped_wavefunction, label='Ground State Wavefunction \u03A8\u2080(x)', color='blue')
ax2.set_ylabel('Wavefunction Amplitude (1/\u221Am)')

ax1.grid(True)
fig.legend()
plt.savefig('stage2_plot.png', dpi=400, bbox_inches='tight')
plt.show()