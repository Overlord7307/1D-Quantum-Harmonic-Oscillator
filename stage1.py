import numpy as np

# Using natural units, hbar = c = omega = m = 1

x_start = -5.0
x_end = 5.0
dx = 0.01


def V(x):
    return 0.5 * (x**2)


def psi_ip1(x, E, psi_i, psi_im1):
    return (x**2 - 2*E) * psi_i * dx**2 + 2*psi_i - psi_im1


def integrator(E_guess):
    psi_values = [0, 0.0001]
    for i in range(2, len(x_values)):
        psi_new = psi_ip1(x_values[i], E_guess, psi_values[-1], psi_values[-2])
        psi_values.append(psi_new)
    return psi_values[-1]

def test_func(x):
    return x**2 - 2

def bisection_search(func, El, Er, max_error):
    Em_values = [Er]
    error = 1

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


x_values = [i for i in np.arange(x_start, x_end, dx)]
V_values = [V(x) for x in x_values]

print(bisection_search(integrator, 0, 1, 0.00001))