# computational-heat

This project explores mathematical modelling and numerical simulation through the study of the heat equation, starting from the idealized 1D rod with Dirichlet boundary conditions and working through the constraints of physical parameters. The analytical solution is derived via separation of variables and Fourier series expansion, while numerical approximations are computed using the Forward-Time Central-Space (FTCS) finite differences method.

## Features
- First-principles derivation of the 1D heat equation from Fourier's law and energy conservation
- Exact analytical solution using separation of variables
- Fourier coefficient derivation for arbitrary initial temperature profiles
- Finite difference discretization grid implementation
- Forward-Time Central-Space (FTCS) derivation and implementation
- Von Neumann Stability analysis and mathematical proof of the $r \leq 1/2$ stability criteria 
- Backward-Time Central-Space (BTCS) derivation and implementation

## Table of Contents
- [1. 1D Heat equation](#1-1d-heat-equation)
  - [1.1 Deriving the equation from first principles](#11-deriving-the-equation-from-first-principles)
  - [1.2 Analytical solution](#12-analytical-solution)
- [2. Forward-Time Central-Space (FTCS)](#2-forward-time-central-space)
  - [2.1 Finite differences](#21-finite-differences)
  - [2.2 FTCS derivation](#22-ftcs-derivation)
  - [2.3 Stability Analysis](#23-stability-analysis)
  - [2.4 Implementation](#24-implementation)
  - [2.5 Results](#25-results)
- [3. Backward-Time Central-Space (BTCS)](#3-backward-time-central-space)
  - [3.1 Derivation](#31-derivation)
  - [3.2 Stability analysis](#32-stability-analysis)
  - [3.3 Thomas' Algorithm](#33-thomas-algorithm)
  - [3.4 Implementation](#34-implementation)
- [4. Future work](#4-future-work)

## 1. 1D Heat equation 

We consider a homogeneous 1D rod of size $[0, L]$, with no internal heat generation and constant physical parameters.      

We further assume that the temperature at the extremities of the rod is null: $T(0,t) = T(L,t) = 0$      
And that the temperature field at t = 0 can be expressed as $T(x,0) = f(x)$     

### 1.1 Deriving the equation from first principles 

Fourier's law tells us:   

$$
q(x,t) = - k \frac{\partial T}{\partial x}
$$

Where $T(x,t)$ is the temperature field, and _k_ is a constant due to the assumption of homogeneity.   

And conservation of Energy tells us:   

$\frac{d}{dt}(Internal \ energy)$ = (heat flux in) - (heat flux out) + (internal generation)   
   
In our case, since a slice of the rod is $[x, x + \Delta x]$ , we can derive the Internal energy as:   

$$
E = \rho cA\Delta x T(x,t)
$$

Where:   

$$
\begin{cases}
\rho \text{ is the density} \\
c \text{ is the specific heat capacity} \\
A\Delta x \text{ is the volume of the slice} \\
T \text{ is the temperature}
\end{cases}
$$

Considering $\rho$, c and A as constants, $\frac{\partial E}{\partial t}$ can thus be expressed as:   

$$
\frac{\partial E}{\partial t} = \frac{d(\rho cA\Delta x T(x,t))}{dt} = \rho cA\Delta x \frac{\partial T}{\partial t}
$$

i.e.,

$$
\rho cA\Delta x \frac{\partial T}{\partial t} = Aq(x,t) - Aq(x + \Delta x, \ t)
$$

Dividing both sides by $A \Delta x$ yields:   

$$
\rho c\frac{\partial T}{\partial t} = \frac{q(x,t) - q(x + \Delta x, \ t)}{\Delta x}
$$

Setting $\Delta x \to 0$ yields the derivative:   

$$
\lim_{\Delta x \to 0}  \frac{q(x,t) - q(x + \Delta x, \ t)}{\Delta x} = -\frac{\partial q}{\partial x}
$$

i.e.,   

$$
\frac{\partial q}{\partial x} = -\rho c\frac{\partial T}{\partial t}
$$

Differentiating Fourier's law with respect to x lands us:   

$$
\frac{\partial q}{\partial x} = -k \frac{\partial^2 T}{\partial x^2 }
$$

Substituting the formula right above yields:   

$$
\rho c \frac{\partial T}{\partial t} = k \frac{\partial^2 T}{\partial x^2 }
$$

i.e.,   

$$
\begin{cases}
\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2 } \\
\alpha = \frac {k}{\rho c}
\end{cases}
$$

Which is exactly the 1-D heat equation.   

### 1.2 Analytical solution 

We now consider the 1D heat equation:   

$$
\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2 }
$$

With the conditions:   

$$
\begin{cases}
T(0,t) = T(L,t) = 0 \\
T(x,0) = f(x)
\end{cases}
$$

We first assume that the solution can be written as two separate functions:

$$
T(x,t) = X(x)H(t)
$$

Plugging this into the heat equation yields:

$$
X(x)H'(t) = \alpha X''(x)H(t)
$$

Dividing both sides by $X(x)H(t)$ gives:

$$
\frac{H'(t)}{H(t)} = \alpha \frac{X''(x)}{X(x)}
$$

For these two functions of x and t to be equal for all (x,t), they must both be equal to the same constant $-\lambda$

i.e.,

$$
\begin{cases}
\frac{1}{\alpha} \frac{H'(t)}{H(t)} = -\lambda \\
\frac{X''(x)}{X(x)} = -\lambda
\end{cases}
$$

The first equation is a homogeneous ODE, and its solution is:

$$
H(t) = C_t e^{-\lambda \alpha t }
$$


The second one is a homogeneous second order ODE:

$$
X''(x) + \lambda X(x) = 0
$$

Solving it requires solving the characteristic equation:

$$
r^2 + \lambda = 0
$$

If $\lambda = 0$ , the solution is:

$$
\begin{cases}
X(x) = Ax + B \\
(A,B) \in \mathbb{R} 
\end{cases}
$$

Since $X(0) = 0$ and $X(L) = 0$ :

$$
X(x) = 0
$$ 

Which is a trivial solution.   

If $\lambda < 0$ : $r = \pm \sqrt{-\lambda}$, the solution is:   

$$
\begin{cases}
X(x) = Ae^{\sqrt{{-\lambda}} x} + Be^{-\sqrt{-\lambda} x} \\
X(0) = 0 \\
X(L) = 0 
\end{cases}
$$

Solving $X(0) = 0$ yields $B = -A$

i.e., 

$$
X(x) = A(e^{\sqrt{{-\lambda}} x} - e^{-\sqrt{{-\lambda}} x})
$$

Plugging in $X(L) = 0$ yields:

$$
A(e^{ \sqrt{-\lambda} L } - e^{-\sqrt{-\lambda} L}) = 0
$$

i.e.,

$$
\begin{cases}
A = 0 \\
X(x) = 0 \\
\end{cases}
$$

Which is again trivial.

If $\lambda > 0$ : $r = \pm i\sqrt{\lambda}$, the solution is therefore:

$$
X(x) = A \cos{\sqrt{\lambda} x} + B \sin{\sqrt{\lambda} x}
$$

Plugging in the initial condition $X(0) = 0$ : 

$$
A = 0
$$

i.e.,

$$
X(x) = B \sin{\sqrt{\lambda} x}
$$

As for the second condition : 

$$
B \sin{\sqrt{\lambda} L } = 0
$$

$B = 0$ yields us again with a trivial solution. Let $B \neq 0$ :

$$
\sin{\sqrt{\lambda} L } = 0
$$

$\Rightarrow$

$$
\begin{cases}
\sqrt{\lambda} L = \pi n \\
n \in \mathbb{N}
\end{cases}
$$

$\Rightarrow$

$$
\lambda_{n} = (\frac{\pi n}{L})^2
$$

i.e.,

$$
\begin{cases}
X_n(x) = b_n \sin(\frac{\pi nx}{L}) \\
B = b_n \text{ for each } X_n
\end{cases}
$$

Finally we have that:

$$
T_n(x,t) = X_n(x)H_n(t) = C_t b_n \sin(\frac{\pi nx}{L}) e^{-(\frac{\pi n}{L})^2 \alpha t } 
$$

Let $C_t b_n = C_n$

$$
T_n(x,t) = C_n \sin(\frac{\pi nx}{L}) e^{-(\frac{\pi n}{L})^2 \alpha t } 
$$

Since the heat equation is linear, a sum of its solutions is also a solution, we can generalise by taking the sum of all the solutions we've found:

$$
T(x,t) = \sum_{n=1}^{\infty} C_n \sin(\frac{\pi nx}{L}) e^{-(\frac{\pi n}{L})^2 \alpha t } 
$$

Plugging in the initial condition $T(x,0) = f(x)$ yields:

$$
f(x) = \sum_{n=1}^{\infty} C_n \sin(\frac{\pi nx}{L})
$$

Which is a Fourier Series, we must derive its coefficient to get the general solution.  

Let $(m,n) \in \mathbb{N^2}$ : 

$$
f(x) = \sum_{n=1}^{\infty} C_n \sin(\frac{\pi nx}{L})
$$

$\Rightarrow$

$$
\sin(\frac{\pi mx}{L}) f(x) = \sin(\frac{\pi mx}{L}) \sum_{n=1}^{\infty} C_n \sin(\frac{\pi nx}{L})
$$

$\Rightarrow$

$$
\int_{0}^{L} \sin(\frac{\pi mx}{L}) f(x) = \int_{0}^{L} \sin(\frac{\pi mx}{L}) \sum_{n=1}^{\infty} C_n \sin(\frac{\pi nx}{L})
$$

$\Rightarrow$

$$
\int_{0}^{L} \sin(\frac{\pi mx}{L}) f(x) dx = \sum_{n=1}^{\infty} C_n \int_{0}^{L} \sin(\frac{\pi mx}{L}) \sin(\frac{\pi nx}{L}) dx
$$

i.e.,

$$
\int_{0}^{L} \sin(\frac{\pi mx}{L}) f(x) dx = C_m \frac{L}{2}
$$

This is because : 

$$
\int_{0}^{L} \sin\left(\frac{m\pi x}{L}\right) \sin\left(\frac{n\pi x}{L}\right)\,dx =
\begin{cases}
0 & m \ne n \\
\frac{L}{2} & m = n
\end{cases}
$$

The Fourier coefficient of f can therefore be expressed as:   

$$
C_n = \frac{2}{L} \int_{0}^{L} \sin(\frac{\pi nx}{L}) f(x) dx 
$$

The general solution to the heat equation under our boundary conditions is therefore:

$$
T(x,t) = \sum_{n=1}^{\infty} \left ( \frac{2}{L}  \int_{0}^{L}  \sin(\frac{\pi nz}{L}) f(z) dz \right)\ \sin(\frac{\pi nx}{L}) e^{-(\frac{\pi n}{L})^2 \alpha t } 
$$

Now that we have the analytical solution, it will serve as a benchmark for validating our numerical solutions.

## 2. Forward-Time Central-Space

As stated previously, the one-dimensional heat equation is:

$$
\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2 }
$$

Our initial and boundary conditions are:

$$
\begin{cases}
T(0,t) = T(L,t) = 0 \\
T(x,0) = f(x)
\end{cases}
$$

The solution $T : [0,L] \times \mathbb{R^+} \to \mathbb{R}$ is a continuous function on $[0,L] \times \mathbb{R^+}$ .

Since it is impossible to calculate the solution at infinitely many points, we first discretize both time and space as such:

$$
\begin{cases}
x_j = j \Delta x \text{ with } j \in \{0, \ldots, N\} \\
\Delta x = \frac{L}{N} \\
\end{cases}
\text{ Where N is the number of spatial intervals}
$$

Similarly, the time grid is defined as such:

$$
\begin{cases}
t^n = n \Delta t \text{ with } n \in \mathbb{N}\\
\Delta t \text{ is the time step}
\end{cases}
$$

The numerical approximation will be denoted as $T_{j}^{n} =T(x_j,t^n)$

### 2.1 Finite differences

The first-order Taylor expansion of a sufficiently smooth function _f_ is:

$$
f(x+h) = f(x) + hf'(x) + O(h^2) 
$$

i.e.,

$$
f'(x) = \frac{f(x + h) - f(x)}{h} + O(h) 
$$

Finite differences use the simple idea that neglecting the $O(h)$ term gives us an approximation of $f'(x)$:

$$
f'(x) \approx \frac{f(x + h) - f(x)}{h} 
$$

This approximation is called the forward difference approximation since it uses "$f(x+h) - f(x)$", it is first-order accurate due to the truncation of the $O(h)$ term.   

Similarly, we can obtain an approximation for the second derivative of f using its fourth-order Taylor expansion:

$$
\begin{cases}
f(x+h) = f(x) + hf'(x) + \frac{h^2}{2} f''(x) + \frac{h^3}{6} f'''(x) + O(h^4) \\
f(x-h) = f(x) - hf'(x) + \frac{h^2}{2} f''(x) - \frac{h^3}{6} f'''(x) + O(h^4)
\end{cases}
$$


Adding these two equations yields:

$$
f(x+h) + f(x-h) =  2 f(x) + h^2 f''(x) + O(h^4)
$$

Solving for $f''(x)$:

$$
f''(x) = \frac{f(x+h) - 2 f(x) + f(x-h)}{h^2} + O(h^2)
$$

Discarding the $O(h^2)$ term yields:

$$
f''(x) \approx \frac{f(x+h) - 2 f(x) + f(x-h)}{h^2}
$$

This is the central difference approximation; it is second-order accurate due to the truncation of the $O(h^2)$ term.

### 2.2 FTCS derivation

Using the forward approximation on the time derivative yields:

$$
\frac{\partial T}{\partial t} \approx \frac{T_{j}^{n+1} - T_{j}^{n}}{\Delta t}
$$

Likewise, using the central difference approximation to compute the second-order space derivative yields:

$$
\frac{\partial^2 T}{\partial x^2 } \approx \frac{T_{j+1}^{n} - 2 T_{j}^{n} + T_{j-1}^{n}}{(\Delta x )^2}
$$

Plugging this into the heat equation:

$$
\frac{T_{j}^{n+1} - T_{j}^{n}}{\Delta t} = \alpha \frac{T_{j+1}^{n} - 2 T_{j}^{n} + T_{j-1}^{n}}{(\Delta x )^2}
$$

Solving for $T_{j}^{n+1}$:

$$
T_{j}^{n+1} = \frac{\alpha \Delta t}{(\Delta x )^2}(T_{j+1}^{n} - 2 T_{j}^{n} + T_{j-1}^{n}) + T_{j}^{n}
$$

Let $r = \frac{\alpha \Delta t}{(\Delta x )^2}$:

$$
T_{j}^{n+1} = r(T_{j+1}^{n} - 2 T_{j}^{n} + T_{j-1}^{n}) + T_{j}^{n}
$$

Rearranging yields:

$$
T_{j}^{n+1} = (1 - 2r) T_{j}^{n} + r(T_{j+1}^{n} + T_{j-1}^{n}) 
$$

Which is the FTCS approximation (Forward time, central space).

### 2.3 Stability Analysis

It is interesting to note that this solution is only conditionally stable: Under certain conditions, the error could grow exponentially.

To properly account for this, we model the error as a Fourier mode:

$$
\epsilon_j^n = \xi^n e^{ik j\Delta x}
$$


Where $\xi$ is the amplification factor and $k$ is the wave number. Stability requires that the error does not grow over time, i.e.:

$$
|\xi| \leq 1
$$

This works because since the error vanishes at the boundaries (just like T does), it satisfies the same conditions that made the sine series a valid basis for the solution in the analytical resolution, the error can therefore be decomposed into that same Fourier basis.

Plugging the error into the FTCS scheme:

$$
\epsilon_{j}^{n+1} = (1 - 2r) \epsilon_{j}^{n} + r(\epsilon_{j+1}^{n} + \epsilon_{j-1}^{n}) 
$$

i.e.,

$$
\xi^{n+1} e^{ik j\Delta x} = (1 - 2r) \xi^n e^{ik j\Delta x} + r(\xi^n e^{ik (j+1)\Delta x} + \xi^n e^{ik (j-1)\Delta x}) 
$$

Dividing both sides by $\xi^{n} e^{ik j\Delta x}$ yields:

$$
\xi = 1 - 2r +r(e^{ik \Delta x} + e^{-ik \Delta x})
$$

i.e.,

$$
\xi = 1 - 2r + 2r \cos(k \Delta x)
$$

$\Rightarrow$

$$
\xi = 1 - 2r(1 - \cos(k \Delta x))
$$

Knowing that $1 - \cos(\theta) = 2 \sin^2 (\frac{\theta}{2})$:

$$
\xi = 1 - 4r\sin^2(\frac{k \Delta x}{2})
$$

We need $|\xi| \leq 1$ i.e.,

$$
-1 \leq 1 - 4r\sin^2(\frac{k \Delta x}{2}) \leq 1
$$

The right side is always verified since $r > 0$, the left one is:

$$
-1 \leq 1 - 4r\sin^2(\frac{k \Delta x}{2})
$$

i.e.,

$$
4r\sin^2(\frac{k \Delta x}{2}) \leq 2
$$

The inequality must be true for all possible $k$, the worst case is $\sin^2(\frac{k \Delta x}{2}) = 1$:

$$
4r \leq 2
$$

i.e.,

$$
r \leq \frac{1}{2}
$$

Which finally gives us the stability condition:

$$
\frac{\alpha \Delta t}{(\Delta x )^2} \leq \frac{1}{2}
$$

i.e.,

$$
\Delta t \leq \frac{(\Delta x )^2}{2 \alpha}
$$

Violating this condition leads to the error growing exponentially, and making the simulation useless.

Since this is true for all Fourier modes, it must be true of any error which can be represented by a Fourier series, which in our case is all of them.

### 2.4 Implementation

To correctly apply the FTCS approximation, we must first model our system as a space time matrix:

$$
\begin{pmatrix}
T_{0}^{0} & T_{0}^{1} & \cdots & \cdots & T_{0}^{N} \\
T_{1}^{0} & \vdots & \vdots & \vdots & T_{1}^{N} \\
\vdots & \vdots & \vdots & \vdots & \vdots \\
\vdots & \vdots & \vdots & \vdots & \vdots \\
T_{I}^{0} & T_{I}^{1} & \cdots & \cdots & T_{I}^{N}
\end{pmatrix}
$$

Which with our boundary conditions translates to:

$$
\begin{pmatrix}
0 & 0 & .. & .. & 0 \\
T_{1}^{0} & .. &.. & .. & T_{1}^{N} \\
.. & .. & .. & .. & .. \\
.. & .. & .. & .. & .. \\
0 & 0 & .. & .. & 0
\end{pmatrix}
$$

The current state is represented by a state vector:

$$
T^n = 
\begin{pmatrix}
T_{0}^{n} \\
T_{1}^{n} \\
\vdots \\
\vdots \\
T_{I}^{n} 
\end{pmatrix}
$$

Since the boundary temperature is fixed, there is no need to update it. The internal state can therefore be written as:

$$
T_{int}^n = 
\begin{pmatrix}
T_{1}^{n} \\
T_{2}^{n} \\
\vdots \\
\vdots \\
T_{I-1}^{n} 
\end{pmatrix}
$$

The update loop can also be written to only include internal points:

$$
T_{int}^{n+1} =  (1 - 2r) T_{int}^{n} + r(T_{left}^{n} + T_{right}^{n})
$$

With:

$$
T_{left}^n = 
\begin{pmatrix}
T_{0}^{n} \\
T_{1}^{n} \\
\vdots \\
\vdots \\
T_{I-2}^{n} 
\end{pmatrix}
\text{ and }
T_{right}^n = 
\begin{pmatrix}
T_{2}^{n} \\
T_{3}^{n} \\
\vdots \\
\vdots \\
T_{I}^{n} 
\end{pmatrix}
$$

Which is exactly the implementation I choose to write:

```python
def heat_ftcs(I, N, r, T_initial):
    # Checking the stability conditions:
    if r > 0.5 :
        raise ValueError(f"Stability condition violated! r must be <= 0.5. Current r = {r:.4f}")
       
    # Initializing the space time grid matrix:
    T = np.zeros((I+1,N+1))

    # Applying the initial conditions:
    T[:,0] = T_initial

    # Forcing the boundary conditions:
    T[0,:] = 0
    T[I,:] = 0

    # FTCS loop:
    for n in range(0,N):
        T[1:I, n+1] = (1 - 2*r) * T[1:I, n] + r * (T[2:I+1, n] + T[0:I-1, n])

    return T
```

### 2.5 Results

As a first step, we will implement the following sinusoidal initial condition:

$$
f(x) = 1000 \sin(\frac{\pi x}{L}) 
$$

Which satisfies the boundary conditions.

The corresponding analytical solution is therefore:

$$
T(x,t) = 1000 \sin(\frac{\pi x}{L}) e^{-(\frac{\pi}{L})^2 \alpha t }
$$

This figure compares the numerical FTCS solution with the analytical solution of the one-dimensional heat equation in the sinusoidal initial condition case. The third heatmap shows the absolute error between the two solutions.

![Ideal Model](images/ftcs_comp.png)

The error is approximately 0.0035 °C at its highest, for a maximum temperature of 1000 °C ;which represents a relative error of 0.00035% ;indicating that the numerical and analytical solutions match closely.

As a second step, we will model the case of a Gaussian initial temperature distribution, i.e.,

$$
f(x) = 1000 \frac{1}{\sigma \sqrt{2\pi}} e^{-\frac{(x-\mu)^2}{2 \sigma^2 } }
$$

With $\sigma$ being the standard deviation and $\mu$ being the mean.

This second figure represents a heat map of the numerical FTCS solution in the case of Gaussian initial conditions with $\mu = 0.5$ and $\sigma = 0.05$.

![Gaussian](images/gaussian.png)

The following figure portrays the spatial Convergence of  the FTCS approximation in the case of the sinusoidal temperature distribution with a constant dt in a log-log scale, we can see that the convergence matches the second order expected rate.

![Spatial Convergence](images/spatial_convergence.png)


## 3. Backward-Time Central-Space 

### 3.1 Derivation

A major drawback of the FTCS method is that it is only conditionally stable, as shown previously, violating the stability condition $\Delta t \leq \frac{(\Delta x )^2}{2 \alpha}$ will cause the error to grow exponentially.

This is why we will introduce the BTCS approximation next. This approximation is built in the same way as FTCS, but instead of evaluating the equation at time step n, we evaluate it at the time step n+1, with backward finite differences in time. i.e.,

$$
\frac{\partial T}{\partial t} |_{n+1} \approx \frac{T_{j}^{n+1} - T_{j}^{n}}{\Delta t}
$$

Likewise, the space derivative is also evaluated at the time step n+1:

$$
\frac{\partial^2 T}{\partial x^2 } |_{n+1} \approx \frac{T_{j+1}^{n+1} - 2 T_{j}^{n+1} + T_{j-1}^{n+1}}{(\Delta x )^2}
$$

Plugging this into the heat equation again yields:

$$
\frac{T_{j}^{n+1} - T_{j}^{n}}{\Delta t} = \alpha \frac{T_{j+1}^{n+1} - 2 T_{j}^{n+1} + T_{j-1}^{n+1}}{(\Delta x )^2}
$$

Defining $r = \frac{\alpha \Delta t}{(\Delta x )^2}$:

$$
T_{j}^{n+1} - T_{j}^{n} = r ( T_{j+1}^{n+1} - 2 T_{j}^{n+1} + T_{j-1}^{n+1} )
$$

Rearranging all the unknown terms on the left side:

$$
T_{j}^{n+1} - r ( T_{j+1}^{n+1} - 2 T_{j}^{n+1} + T_{j-1}^{n+1} ) = T_{j}^{n} 
$$

i.e.,

$$
-r T_{j+1}^{n+1} + (1+2r) T_{j}^{n+1} - r T_{j-1}^{n+1} = T_{j}^{n}
$$

Which is exactly the BTCS approximation.

Defining the same internal state vector as before:

$$
T^n = 
\begin{pmatrix}
T_{1}^{n} \\
T_{2}^{n} \\
\vdots \\
\vdots \\
T_{I-1}^{n} 
\end{pmatrix}
$$

This gives us a system of equations:

$$
\begin{cases}
-r T_{2}^{n+1} + (1+2r) T_{1}^{n+1} - r T_{0}^{n+1} = T_{1}^{n} \\
-r T_{3}^{n+1} + (1+2r) T_{2}^{n+1} - r T_{1}^{n+1} = T_{2}^{n} \\
\vdots \\
\vdots \\
-r T_{I}^{n+1} + (1+2r) T_{I-1}^{n+1} - r T_{I-2}^{n+1} = T_{I-1}^{n}
\end{cases}
$$

Which exactly translates to:

$$
AT^{n+1} = T^n + b
$$

With:

$$
A =
\begin{pmatrix}
(1+2r) & -r & 0 & 0 & \cdots & 0 \\
-r & (1+2r) & -r & 0 & 0 & 0 \\
0 & -r & (1+2r) & -r  & 0 & \vdots \\
0 & 0 & -r & \ddots & \ddots & 0 \\
\vdots & \vdots & 0 & \ddots & \ddots & -r \\
0 & 0 & 0 & 0 & -r & (1+2r) \\
\end{pmatrix}
\text{ and }
b =
\begin{pmatrix}
r T_{0}^{n+1} \\
0 \\
0  \\
.. \\
0 \\
r T_{I}^{n+1}  \\
\end{pmatrix}
$$

Since our boundary temperature is always 0, the system simplifies to:

$$
AT^{n+1} = T^n
$$

It is worth noting that unlike FTCS, BTCS is an implicit method: solving for the next time step requires us to solve a linear system, rather than directly computing the new temperatures from the previous time step.

### 3.2 Stability analysis

We will again be modeling the error as a Fourier mode:

$$
\epsilon_j^n = \xi^{n} e^{ik j\Delta x}
$$

Where $\xi$ is the amplification factor and $k$ is the wave number.

Plugging this into our BTCS scheme:

$$
-r \epsilon_{j+1}^{n+1} + (1+2r) \epsilon_{j}^{n+1} - r \epsilon_{j-1}^{n+1} = \epsilon_{j}^{n}
$$

i.e.,

$$
-r \xi^{n+1} e^{ik (j+1)\Delta x} + (1+2r) \xi^{n+1} e^{ik j\Delta x} - r \xi^{n+1} e^{ik (j-1)\Delta x} = \xi^{n} e^{ik j\Delta x}
$$

Dividing both sides by $\xi^{n} e^{ikj \Delta x}$ yields:

$$
-r \xi e^{ik \Delta x} + (1+2r) \xi  - r \xi e^{-ik\Delta x} = 1
$$

i.e.,

$$
\xi ( -r e^{ik \Delta x} + (1+2r) - r e^{-ik\Delta x} )= 1
$$

$\Rightarrow$

$$
\xi ( (1+2r) -2r \cos (k\Delta x) )= 1
$$

$\Rightarrow$

$$
\xi = \frac{1}{(1+2r) -2r \cos (k\Delta x)}
$$

$\Rightarrow$

$$
\xi = \frac{1}{1+2r(1- \cos (k\Delta x))} \geq 0
$$

Since $-1 \leq \cos (k\Delta x) \leq 1$ we have that:

$$
0 \leq 1 - \cos (k\Delta x) \leq 2
$$

And since $r \in \mathbb{R^+}$:

$$
2r (1- \cos (k\Delta x)) \geq 0
$$

i.e.,

$$
1+ 2r (1- \cos (k\Delta x)) \geq 1
$$

Which finally yields:

$$
| \xi | \leq 1
$$

For all positive _r_ ; the BTCS approximation is thus unconditionally stable, meaning that no restriction on the timestep is required for the simulation to remain stable.

### 3.3 Thomas' Algorithm

The BTCS approximation requires solving the linear system:

$$
AT^{n+1} = T^n
$$

Where A is a tridiagonal matrix: a matrix with zeros everywhere except on its main diagonal and the two adjancent ones.

Instead of using a general purpose linear solver, we will take advantage of the fact that A is tridiagonal and use a specific algorithm made for this exact kind of matrices.

Thomas' algorithm is an adaptation of Gaussian elimination designed specifically for tridiagonal systems, it works by only storing and operating on the three non-zero diagonals:

$$
A =
\begin{pmatrix}
b_1 & c_1 & 0 & 0 & \cdots & 0 \\
a_2 & b_2 & c_2 & 0 & \cdots & 0 \\
0 & a_3 & b_3 & c_3 & \cdots & 0 \\
\vdots & \vdots & \ddots & \ddots & \ddots & \vdots \\
0 & 0 & \cdots & a_{n-1} & b_{n-1} & c_{n-1} \\
0 & 0 & \cdots & 0 & a_n & b_n
\end{pmatrix}
$$

where _b_ is the main diagonal, and _a_ and _c_ the adjacent ones.

### 3.4 Implementation

```Python
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
```


## 4. Future Work

- Implementing Crank-Nicolson
- Solving the 2D heat equation

