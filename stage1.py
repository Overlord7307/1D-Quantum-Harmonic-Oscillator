'''
In stage 1 of building a 1D quantum harmonic oscillator, we define the boundaries, then integrate the Schrodinger wave function using the finite difference method for the second derivative. We use the bisection search by varying energy values to check if the wave function goes to zero at the boundaries. Note that for simplicity, we have set all constants to unity (natural units). Thus we are able to obtain both the ground state energy and its corresponding wavefunction. The mathematics is detailed in the README.md file.
'''

import numpy as np
import matplotlib.pyplot as plt

# Using natural units, ℏ = ω = m = 1

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
    return 2*(V(x) - E) * ψ_i * dx**2 + 2*ψ_i - ψ_im1


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
x_axis = [i for i in np.arange(-4.0, 4.0, dx)]
V_values = [V(x) for x in x_axis]

Energy_Ground = bisection_search(integrator, 0, 1, 0.0000001)
wavefunction = wavefunction_generator(Energy_Ground)

# Calculating the normalization constant
area_under_ψ = 0
for val in wavefunction:
    area_under_ψ += (val**2) * dx
N = area_under_ψ ** 0.5

normalized_wavefunction = [val/N for val in wavefunction]
chopped_wavefunction = normalized_wavefunction[100:-100]    # Chopping off the ends for a better plot

print('The energy of the ground state is:', Energy_Ground, '\u210F\u03C9')
print()

# Plotting the harmonic potential, the ground state wavefunction and the ground state energy against the x-values
plt.figure(figsize=(16, 9))
plt.plot(x_axis, V_values, label='Harmonic Potential V(x)', color='black', linestyle='--', linewidth=1)
plt.plot(x_axis, chopped_wavefunction, label='Ground State Wavefunction \u03A8\u2080(x)', color='blue')
plt.title('Ground State Wavefunction for a Quantum Harmonic Oscillator')
plt.axhline(Energy_Ground, label=f'Ground State Energy E\u2080 = {Energy_Ground:5f} \u210F\u03C9', color='red')
plt.axhline(0, color='black', linewidth=0.75)
plt.axvline(0, color='black', linewidth=0.75)
plt.grid(True)
plt.legend()
plt.savefig('plots\\stage1_plot.png', dpi=400, bbox_inches='tight')
plt.show()