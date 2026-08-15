# 1D Quantum Harmonic Oscillator

In this project, I have built a 1D quantum harmonic oscillator from scratch in order to numerically calculate its energy eigenvalues and the corresponding wavefunctions. I have divided the project into stages, starting out with calculating the ground state using simple techniques, and increasing in complexity going further.

The Hamiltonian for a 1-dimensional QHO is

$$\begin{equation}
    \hat{H} = -\frac{\hbar^2}{2m} \frac{\mathrm{d}^2}{\mathrm{d}x^2} + \hat{V}(x).
\end{equation}$$

Plugging in the harmonic potential $\hat{V}(x) = \frac{1}{2} m\omega^2\hat{x}^2$ and $\hat{H}$ into the Schrödinger equation $\hat{H}\psi = E\psi$, we get:

$$\begin{equation}
    -\frac{\hbar^2}{2m}\frac{\mathrm{d}^2\psi}{\mathrm{d}x^2} + \frac{1}{2}m\omega^2x^2\psi = E\psi.
\end{equation}$$


## Stage 1

In the first stage, the goal is very simple: to find and plot the ground state energy eigenvalue and ground state wavefunction for the QHO using only pure Python and Matplotlib (little to no usage of NumPy and SciPy). Of course, not using any external libraries means the code will be very inefficient. However, since the calculations are not very complicated, and since the main point here is to understand the physics and the computational methods required for this solution, I'm not very bothered by this. To make things even easier (and to avoid floating point errors), I will also be using **natural units**, i.e. $\hbar = \omega = m = 1$. This means that our final energy will be in terms of units of $\hbar\omega$. Substituting these natural units in the above Schrödinger equation, we get:

$$\begin{equation}
    -\frac{1}{2} \frac{\mathrm{d}^2\psi}{\mathrm{d}x^2} + \frac{1}{2}x^2\psi = E\psi,
\end{equation}$$

and rearranging we have:

$$\begin{equation}
    \frac{\mathrm{d}^2\psi}{\mathrm{d}x^2} = 2(V(x) - E)\psi.
\end{equation}$$

Now, instead of solving this 2nd order differential equation analytically, which would defeat the whole purpose of the project, we will employ the **finite difference method** for second derivatives. This method provides the following approximation for the second derivative:

$$\begin{equation}
    \frac{\mathrm{d}^2\psi}{\mathrm{d}x^2} \approx \frac{\psi_{i+1} - 2\psi_i + \psi_{i-1}}{\mathrm{d}x^2}.
\end{equation}$$

Substituting into (4) and rearranging the terms, we can calculate the iterative value of $\psi_{i+1}$ using the two previous values:

$$\begin{equation}
    \psi_{i+1} = 2\left[\left(V(x) - E\right)\mathrm{d}x^2 + 1\right]\psi_i - \psi_{i-1}.
\end{equation}$$

Finally we have the relation we need to implement an integrator function, in order to iteratively find the values of $\psi(x)$ for all values of $x$.

Our aim now is to use a bisection search to find the correct value of $E$ for which $\psi(x)$ disappears at the boundaries. Since we know that the energy eigenvalues for a QHO are given by $E_n = \left(n + \frac{1}{2}\hbar\omega \right),$ we can set limits for our bisection search at 0 and 1 and be confident that the answer we get will indeed be $E_0$. Once we obtain the correct ground state energy, we can plug that back into the integrator and find the entire wavefunction. After this, plotting $V(x)$, $\psi(x)$ and $E_0$ against $x$ is a trivial matter and can be done easily using Matplotlib.


## Stage 2

Our goal for the 2nd stage is firstly, to present the final answers and plots in SI units. To accomplish this, we initially calculate everything in natural units as before, and then scale the results back to SI units for the purposes of plotting by using the following transformations:
$$\begin{align*}
    x &\rightarrow Lx, \\
    V(x) &\rightarrow m\omega^2V(x), \\
    E &\rightarrow \hbar\omega E, \\
    \psi(x) &\rightarrow \frac{\psi(x)}{\sqrt{L}},
\end{align*}$$
where the characteristic length of the harmonic oscillator is $L = \sqrt{\frac{\hbar}{m\omega}}$, and the physical constants $m$, $\hbar$ and $\omega$ must be in SI units.

Second and more importantly, we want to upgrade the physics engine to support higher order eigenvalues and eigenstates. We know that the number of nodes formed by a wavefunction directly corresponds to the state of the system it describes. We can use this property to create a node-counting algorithm which counts the number of nodes for each wavefunction corresponding to every energy value starting from $E = 0$ and increasing by increments of $dE = 0.1$ until the number of nodes surpasses the target $n$ (quantum state) value. Since energy levels are quantized, the number of nodes will jump to the next integer as soon as the test energy exceeds the current state's energy. This will give us an approximate energy bracket around the required node, which we can plug into our bisection search algorithm in order to find the exact energy eigenvalue.

Thus, after implementing this stage we are able to accurately determine and plot the energy eigenvalue and the wavefunction for any given state of the system as long as the wavefunction has sufficient space to go to zero in the boundary conditions we set initially.