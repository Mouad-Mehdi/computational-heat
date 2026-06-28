# computational-heat
Numerical resolution of the 1D and 2D heat equation using finite difference methods.

To properly model the heat equation, we must first understand where it comes from.    

## 1D Heat equation 

We consider a homogeneous 1D rod of size $[0, L]$, with no internal heat generation and constant physical parameters.      

We further assume that the temperature at the extremities of the rod is null: $T(0,t) = T(L,t) = 0$      
And that the temperature field at t = 0 can be expressed as $T(x,0) = f(x)$     

### Deriving the equation from first principle 

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

### Analytical solution 

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

We first assume that the solution can be written as two seperate functions:

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

for these two functions of x and t to be equal for all (x,t), they must both be equal to the same constant $-\lambda$

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
(A,B) \in \mathbf{R} 
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
n \in \mathbf{N}
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

Let $(m,n) \in \mathbf{N^2}$ : 

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
\int_{0}^{L} \sin(\frac{\pi mx}{L}) f(x) dx = \sum_{n=1}^{\infty} \int_{0}^{L} C_n\sin(\frac{\pi mx}{L}) \sin(\frac{\pi nx}{L}) dx
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

the Fourier coefficient of f can therefore be expressed as:   

$$
C_n = \frac{2}{L} \int_{0}^{L} \sin(\frac{\pi nx}{L}) f(x) dx 
$$

The general solution to the heat equation under our boundary conditions is therefore:

$$
T(x,t) = \sum_{n=1}^{\infty} \left ( \frac{2}{L}  \int_{0}^{L}  \sin(\frac{\pi nz}{L}) f(z) dz \right)\ \sin(\frac{\pi nx}{L}) e^{-(\frac{\pi n}{L})^2 \alpha t } 
$$
