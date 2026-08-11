# Quantum Harmonic Oscillator

In this project, I have built a quantum harmonic oscillator from scratch in order to numerically calculate its energy eigenvalues and the corresponding wavefunctions. I have divided the project into stages, starting out with calculating the ground state using simple techniques, and increasing in complexity going further.

The Hamiltonian for a 1-dimensional QHO is

$$\hat{H} = -\frac{\hbar^2}{2m} \frac{\mathrm{d}^2}{\mathrm{d}x^2} + \hat{V}(x).$$

Plugging in the harmonic potential $\hat{V}(x) = \frac{1}{2} m\omega^2\hat{x}^2$ and $\hat{H}$ into the Schrödinger equation $\hat{H}\psi = E\psi$, we get:

$$-\frac{\hbar^2}{2m}\frac{\mathrm{d}^2\psi}{\mathrm{d}x^2} + \frac{1}{2}m\omega^2x^2\psi = E\psi.$$


## Stage 1

In the first stage, the goal is very simple: to find and plot the ground state energy eigenvalue and ground state wavefunction for the QHO using only pure Python and Matplotlib (little to no usage of NumPy and SciPy). Of course, not using any external libraries means the code will be very inefficient. However, since the calculations are not very complicated, and since the main point here is to understand the physics and the computational methods required for this solution, I'm not very bothered by this. To make things even easier (and to avoid floating point errors), I will also be using **natural units**, i.e. $\hbar = \omega = m = 1$. This means that our final energy will be in terms of units of $\hbar\omega$. Substituting these natural units in the above Schrödinger equation, we get:

$$-\frac{1}{2} \frac{\mathrm{d}^2\psi}{\mathrm{d}x^2} + \frac{1}{2}x^2\psi = E\psi,$$

and rearranging we have:

$$\frac{\mathrm{d}^2\psi}{\mathrm{d}x^2} = 2(V(x) - E)\psi.$$

Now, instead of solving this 2nd order differential equation analytically, which would defeat the whole purpose of the project, we will employ the **finite difference method** for second derivatives. This method provides the following approximation for the second derivative:

$$\frac{\mathrm{d}^2\psi}{\mathrm{d}x^2} \approx \frac{\psi_{i+1} - 2\psi_i + \psi_{i-1}}{\mathrm{d}x^2}.$$

Substituting into the Schrödinger equation and rearranging the terms, we can calculate the value of $\psi_{i+1}$ using the two previous values:

$$\psi_{i+1} = 2[(V(x) - E)\mathrm{d}x^2 + 1]\psi_i - \psi_{i-1}.$$

Finally, we have the relation we need to implement an integrator function in order to iteratively find the values of $\psi(x)$ for all values of $x$.

Our aim now is to use a bisection search to find the correct value of $E$ for which $\psi(x)$ disappears at the boundaries. Since we know that the energy eigenvalues for a QHO are given by $E_n = \left(n + \frac{1}{2}\hbar\omega \right),$ we can set limits for our bisection search at 0 and 1 and be confident that the answer we get will indeed be $E_0$. Once we obtain the correct ground state energy, we can plug that back into the integrator and find the entire wavefunction. After this, plotting $V(x)$, $\psi(x)$ and $E_0$ against $x$ is a trivial matter and can be done easily using Matplotlib.