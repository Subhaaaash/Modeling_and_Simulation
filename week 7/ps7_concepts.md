---
title: "Modeling and Simulation: Concept Sheet for Problem Sheet 7"
subtitle: "Conserved quantities, phase space, semi-implicit Euler and Runge-Kutta 2 (WS 2026/27)"
author: "Subhash"
---

This sheet adds only the **new** ideas you need for problem sheet 7. Ideas that the earlier concept sheets already explain are **not repeated**. Here is where to find them:

| Already covered | Where |
|:--|:-:|
| Converting an $n$-th order ODE to a first-order system | week 5, section 4 |
| Forward (explicit) Euler: derivation, $\mathcal O(h)$ error | week 5, section 7 |
| Programming pattern: time grid, arrays, initial condition, loop | week 5, section 9 |
| Euler methods for **systems**, matrix form $\vec x' = A\vec x$ | week 6, sections 2, 3 |
| Backward (implicit) Euler, solving the 2×2 system by hand | week 6, section 4 |
| Trapezoidal method | week 6, section 5 |
| Stability, test equation $y' = \lambda y$, amplification factor $G$, A-stability, eigenvalues | week 6, section 6 |
| Energy as a check quantity (dissipative case) | week 6, section 8 |
| `plt.subplots`, time grid with `np.arange` | week 6, section 10 |

**New in this sheet:**

| Exercise part | New concepts | Section |
|:--|:--|:-:|
| (a) | conserved quantities, proving conservation with the derivative, phase space, circles and ellipses | 1, 2 |
| (b) | the harmonic oscillator as a system, semi-implicit (Euler-Cromer) method, RK2 midpoint method | 3, 4, 5 |
| (c) | trajectory plots in the $(x, y)$ plane, markers, equal axes | 7 |
| (d) | stability for purely imaginary eigenvalues, energy growth per step, "stable" vs. "energy-conserving", symplectic methods | 6 |

Slide numbers refer to the lecture *Part 3: Continuous models: Ordinary differential equations* (sections 9.6 to 9.11).

# 1. Conserved quantities

## 1.1 What "conserved" means

A quantity $E(x, y)$ is **conserved** (a *first integral* or *invariant* of the ODE) if it keeps its initial value along every solution:
$$
E\big(x(t), y(t)\big) = E\big(x(0), y(0)\big) \quad \text{for all } t .
$$

In mechanics this is usually the total energy of a system without friction. In week 6 the energy of the RLC circuit could only **decrease** ($dE/dt = -Ri^2$). Without a resistor ($R = 0$) it would be conserved.

## 1.2 How to prove it: differentiate along the solution

You do **not** need the solution itself. Differentiate $E$ with the chain rule and replace every $\dot x$, $\dot y$ by the right-hand side of the ODE:
$$
\frac{dE}{dt} = \frac{\partial E}{\partial x}\,\dot x + \frac{\partial E}{\partial y}\,\dot y
\;\overset{\text{ODE}}{=}\; \dots \;=\; 0 .
$$

If the result is zero for **every** $(x, y)$, then $E$ is constant on every solution.

**Worked example.** A spring with angular frequency $\omega = 2$: $\dot x = y$, $\dot y = -4x$. Try $E = 4x^2 + y^2$:
$$
\frac{dE}{dt} = 8x\,\dot x + 2y\,\dot y = 8x\,y + 2y\,(-4x) = 0 .
$$

So $E$ is conserved. The factor 4 in front of $x^2$ is needed: with $x^2 + y^2$ the two terms would not cancel. In general, for $\dot x = y$, $\dot y = -\omega^2 x$, the conserved quantity is $\omega^2 x^2 + y^2$ (twice the kinetic plus potential energy per mass).

**Typical mistakes.**

- Forgetting the chain rule factor: $\frac{d}{dt}(x^2) = 2x\,\dot x$, not $2x$.
- Inserting the exact solution and checking only one time point. That shows the value at that time, not that it is constant.

# 2. Phase space and the geometry of level curves

## 2.1 The phase plane

For a system with two state variables $(x, y)$ you can plot the solution in two ways:

- **time plots:** $x(t)$ and $y(t)$ against $t$ (as in week 6);
- **phase plot (trajectory):** $y$ against $x$, with time as a parameter along the curve.

The $(x, y)$ plane is called the **phase plane** (for $y = \dot x$: position and velocity). Each initial condition gives one curve, and a numerical method gives a sequence of points $(x_n, y_n)$ on (or near) it.

## 2.2 A conserved quantity fixes the curve

If $E$ is conserved, the trajectory must stay on the **level curve** $E(x, y) = E(x_0, y_0)$:

| Conserved quantity | Level curve $= c > 0$ | Shape |
|:--|:--|:--|
| $x^2 + y^2$ | $x^2 + y^2 = c$ | **circle**, centre $(0, 0)$, radius $r = \sqrt{c}$ |
| $\omega^2 x^2 + y^2$ | $\dfrac{x^2}{c/\omega^2} + \dfrac{y^2}{c} = 1$ | **ellipse**, semi-axes $\sqrt c/\omega$ and $\sqrt c$ |

So in the phase plane, **energy error = distance from the correct curve**. For the circle, $\sqrt{E_n} = \sqrt{x_n^2 + y_n^2}$ is simply the distance of the point from the origin. You can read the energy error directly off a trajectory plot.

## 2.3 Direction and speed of motion

The ODE tells you the velocity vector $(\dot x, \dot y)$ at every point. For $\dot x = y$, $\dot y = -\omega^2 x$ at the start point $(1, 0)$: $(\dot x, \dot y) = (0, -\omega^2)$, so the point moves **downwards**, i.e. **clockwise**. One full round takes one period $T = 2\pi/\omega$.

Useful for reading plots: with step size $h$, one round contains about $T/h$ points. If the plot shows a polygon with few corners, $h$ is large compared with the period.

# 3. The harmonic oscillator as a first-order system

The equation $\ddot x = -\omega^2 x$ (spring, LC circuit, small pendulum swings) becomes, with $y = \dot x$ (week 5, section 4):
$$
\dot x = y, \qquad \dot y = -\omega^2 x
\qquad\Longleftrightarrow\qquad
\frac{d}{dt}\begin{pmatrix}x\\y\end{pmatrix} =
\begin{pmatrix}0 & 1\\ -\omega^2 & 0\end{pmatrix}
\begin{pmatrix}x\\y\end{pmatrix}.
$$

- **Exact solution** for $x(0) = x_0$, $y(0) = 0$: $x = x_0\cos\omega t$, $y = -x_0\,\omega\sin\omega t$.
- **Eigenvalues** of the matrix: $\lambda = \pm i\,\omega$. They are **purely imaginary** ($\operatorname{Re}\lambda = 0$): the exact solution neither grows nor decays. This is the borderline case of the stability analysis in week 6, section 6, and it is the reason why oscillators are a hard test for integrators (section 6 below).

# 4. Semi-implicit Euler (Euler-Cromer) method (lecture 9.6)

## 4.1 The scheme

For a pair of ODEs of the form (slide "Euler-Cromer")
$$
v' = f(t, y, v), \qquad y' = g(t, v) \quad (g \text{ does not depend on } y),
$$
use a **forward** step for $v$ and then a **backward-looking** step for $y$ that uses the **new** $v$:
$$
\boxed{\;
\begin{aligned}
v_{n+1} &= v_n + h\,f(t_n, y_n, v_n)\\
y_{n+1} &= y_n + h\,g(t_{n+1}, v_{n+1})
\end{aligned}
\;}
$$

This form appears whenever a second-order equation $y'' = f$ is written as a system with $v = y'$ (Newton's law: velocity and position).

## 4.2 Why it is "semi-implicit" but costs nothing extra

- The first line is an ordinary forward Euler step.
- The second line looks like backward Euler, because it uses a value at $t_{n+1}$. But $v_{n+1}$ is **already known** from the first line, so there is **no equation to solve**.

So you get some of the good behaviour of implicit methods at the price of forward Euler. The only change compared with forward Euler is the **order** of the two update lines and which value is used.

## 4.3 Worked example by hand

$\dot x = y$, $\dot y = -4x$ (here $x$ plays the role of the position and $y$ the velocity), $x_0 = 1$, $y_0 = 0$, $h = 0.1$:
$$
y_{n+1} = y_n - 4h\,x_n, \qquad x_{n+1} = x_n + h\,y_{n+1}.
$$

| $n$ | forward Euler $(x_n,\ y_n)$ | $E_n$ | semi-implicit $(x_n,\ y_n)$ | $E_n$ | exact $(x,\ y)$ |
|:-:|:-:|:-:|:-:|:-:|:-:|
| 0 | (1, 0) | 4 | (1, 0) | 4 | (1, 0) |
| 1 | (1.0000, −0.4000) | 4.160 | (0.9600, −0.4000) | 3.846 | (0.9801, −0.3973) |
| 2 | (0.9600, −0.8000) | 4.326 | (0.8816, −0.7840) | 3.724 | (0.9211, −0.7788) |
| 3 | (0.8800, −1.1840) | 4.500 | (0.7679, −1.1366) | 3.651 | (0.8253, −1.1293) |

with $E = 4x^2 + y^2$. Forward Euler increases $E$ by the same factor $1 + (2h)^2 = 1.04$ in every step. The semi-implicit energy first goes down, but later it comes back up: it **oscillates** around the correct value instead of drifting away (section 6.3).

## 4.4 Typical programming mistakes

```python
# WRONG: this is forward Euler, both lines use old values
y[n + 1] = y[n] - h * x[n]
x[n + 1] = x[n] + h * y[n]

# RIGHT: semi-implicit, the second line uses the NEW y
y[n + 1] = y[n] - h * x[n]
x[n + 1] = x[n] + h * y[n + 1]
```

The order of the two lines matters. You may also update $x$ first and then $y$ with the new $x$; that is the other semi-implicit variant and behaves the same way. Just use one variant consistently and say which one you chose.

## 4.5 Properties (lecture summary)

- First order: local error $\mathcal O(h^2)$, global error $\mathcal O(h)$, like both Euler methods.
- Much better than forward Euler for oscillations: the energy error stays **bounded** (section 6.3).

# 5. Second-order Runge-Kutta (midpoint) method (lecture 9.7)

## 5.1 Idea

Forward Euler uses the slope at the **start** of the step. A much better slope for the whole step is the slope at the **midpoint** $t_n + h/2$. That midpoint value is not known, so first make a half Euler step to estimate it:
$$
\boxed{\;
\begin{aligned}
k_1 &= f(t_n,\ y_n) && \text{slope at the start}\\
k_2 &= f\!\left(t_n + \tfrac h2,\ y_n + \tfrac h2\,k_1\right) && \text{slope at the estimated midpoint}\\
y_{n+1} &= y_n + h\,k_2 && \text{full step with the midpoint slope}
\end{aligned}
\;}
$$

The method is **explicit** (no equation to solve), needs **two** evaluations of $f$ per step, and is **second order**: local error $\mathcal O(h^3)$, global error $\mathcal O(h^2)$.

## 5.2 For a system: the $k$'s are vectors

For $\vec x' = \vec f(t, \vec x)$ every $k$ has one component per variable. Compute **all** components of $k_1$ first, then the midpoint for **all** variables, then **all** components of $k_2$:

```python
k1x = f_x(x[n], y[n]);             k1y = f_y(x[n], y[n])        # 1. both slopes at the start
x_mid = x[n] + 0.5 * h * k1x;      y_mid = y[n] + 0.5 * h * k1y # 2. both midpoint values
k2x = f_x(x_mid, y_mid);           k2y = f_y(x_mid, y_mid)      # 3. both slopes at the midpoint
x[n + 1] = x[n] + h * k2x;         y[n + 1] = y[n] + h * k2y    # 4. full step
```

Mixing the order (e.g. computing `x_mid` and then using it for `k1y`) gives a different, wrong scheme.

## 5.3 Worked example by hand

$y' = -2y$, $y_0 = 1$, $h = 0.1$, the same example as in week 6, sections 4 and 5:
$$
k_1 = -2y_n,\quad k_2 = -2\,(y_n + 0.05\,k_1) = -2y_n + 0.2\,y_n = -1.8\,y_n,\quad
y_{n+1} = y_n + 0.1\,k_2 = 0.82\,y_n .
$$

| $t$ | RK2 midpoint | exact $e^{-2t}$ | error |
|:-:|:-:|:-:|:-:|
| 0.1 | 0.8200 | 0.8187 | $1.3\times10^{-3}$ |
| 0.2 | 0.6724 | 0.6703 | $2.1\times10^{-3}$ |
| 0.3 | 0.5514 | 0.5488 | $2.6\times10^{-3}$ |

The error is more than 10 times smaller than for the two Euler methods in week 6 (0.03 to 0.04 at $t = 0.3$). It is about twice the error of the trapezoidal method (0.0011), which is also second order, but RK2 gets there **without solving an equation**.

**For linear systems.** With $\vec x' = A\vec x$ one step is
$$
\vec x_{n+1} = \left(I + hA + \tfrac{h^2}{2}A^2\right)\vec x_n ,
$$
i.e. the Taylor series of $e^{hA}$ up to the $h^2$ term. Forward Euler stops after $hA$. This form is handy for checking a program and for the analysis in section 6.

## 5.4 Related schemes from the lecture (for orientation)

| Method | Slopes used | Order | $f$-evaluations per step |
|:--|:--|:-:|:-:|
| Forward Euler | start | 1 | 1 |
| RK2 midpoint (9.7) | start → midpoint | 2 | 2 |
| RK2 Heun (9.8) | average of start and (predicted) end | 2 | 2 |
| RK4 (9.9) | weighted average of 4 slopes, $\tfrac h6(k_1 + 2k_2 + 2k_3 + k_4)$ | 4 | 4 |

Heun's method is the explicit version of the trapezoidal method (week 6): it predicts $y_{n+1}$ with forward Euler and uses it in the end slope.

**Why higher order pays off (lecture 9.11).** For a global error tolerance $\delta$ and a method of order $p$, the number of steps grows like $\delta^{-1/p}$. For $\delta = 10^{-8}$ on $[0, 1]$: Euler needs about $10^8$ steps, RK2 about $10^4$ steps ($2\times10^4$ evaluations), RK4 about $10^2$ steps (400 evaluations). The extra evaluations per step are far cheaper than the extra steps.

# 6. Stability and energy for oscillations

## 6.1 Purely imaginary eigenvalues: the borderline case

In week 6 the decisive quantity was the amplification factor $G(z)$ with $z = h\lambda$. For an oscillator $\lambda = \pm i\omega$, so $z = \pm i\,h\omega$ lies **on the imaginary axis**. The exact factor is $|e^{i h\omega}| = 1$: the amplitude is kept exactly. For a numerical method:

- $|G| > 1$: the amplitude grows in every step, the trajectory **spirals outwards**. The method is unstable for this problem, even if the growth is slow.
- $|G| < 1$: the amplitude shrinks, the trajectory **spirals inwards**. The method is stable, but it adds **artificial damping**.
- $|G| = 1$: the amplitude is kept. This is what we want.

## 6.2 Energy factor per step

For $E = \omega^2 x^2 + y^2$ (proportional to the squared amplitude) one step multiplies $E$ by $|G|^2$. With $z = h\omega$:

| Method | $\lvert G\rvert^2$ = energy factor per step | after many steps |
|:--|:-:|:--|
| exact | 1 | circle/ellipse |
| Forward Euler | $1 + z^2$ | grows for **every** $h$: spiral out |
| Backward Euler | $\dfrac{1}{1 + z^2}$ | decays for every $h$: spiral in |
| Trapezoidal (week 6) | 1 | stays on the curve |
| RK2 midpoint | $\left(1 - \tfrac{z^2}{2}\right)^2 + z^2 = 1 + \dfrac{z^4}{4}$ | grows, but **very slowly** for small $z$ |
| Semi-implicit Euler | not a constant factor (see 6.3) | oscillates, **bounded** for $z < 2$ |

After $N$ steps the factor is $|G|^{2N}$. A tiny excess like $1 + z^4/4$ still adds up: compare $z^2$ (Euler) with $z^4/4$ (RK2) for $z = 0.1$: $10^{-2}$ vs. $2.5\times10^{-5}$ per step. A second-order method is not only more accurate, it also drifts much more slowly. (The full stability region of RK2 in the complex plane is a task in problem sheet 8.)

## 6.3 Why semi-implicit Euler keeps the energy bounded

For $\dot x = y$, $\dot y = -\omega^2 x$ with the update $y_{n+1} = y_n - h\omega^2 x_n$, $x_{n+1} = x_n + h\,y_{n+1}$, a short calculation shows that the **modified energy**
$$
\tilde E = \omega^2 x^2 + y^2 - h\,\omega^2\,x\,y
$$
is conserved **exactly** by the numerical scheme: $\tilde E_{n+1} = \tilde E_n$. (Check it with the table in section 4.3: after step 1, $4\cdot0.9216 + 0.16 - 0.4\cdot0.96\cdot(-0.4) = 4.000$.)

Since $\tilde E$ differs from $E$ only by the small term $h\omega^2 xy$, the true energy $E$ can only **wobble** around its initial value by a relative amount of order $h\omega$, but it can **never drift** away. The numerical trajectory is a closed, slightly tilted ellipse close to the exact curve.

A second view: the step matrix of semi-implicit Euler has determinant **1**, so it preserves **areas** in the phase plane. Methods with this property are called **symplectic**. Forward Euler (determinant $1 + z^2$) blows areas up, backward Euler shrinks them. Symplectic methods are the standard choice for long simulations of frictionless systems (planets, molecules); the Verlet scheme (lecture 9.12) is another one.

## 6.4 "Stable" and "energy-conserving" are different questions

| Question | What to look at |
|:--|:--|
| Is the method **stable**? | Does the numerical solution stay **bounded** (no spiral to infinity)? |
| Does it **conserve** $E$? | Does $E_n$ stay **close to $E_0$** for all $n$, without drift? |

Backward Euler is stable but does **not** conserve energy. A method can also look fine for a short run and still drift (RK2 with small $h$). When you answer such a question, give both the observation from the plot (spiral out/in, closed curve) and the reason (energy factor per step).

# 7. Plotting trajectories in the phase plane

```python
phi = np.linspace(0, 2 * np.pi, 400)                 # reference curve: unit circle
fig, ax = plt.subplots(figsize=(7, 7))
ax.plot(np.cos(phi), np.sin(phi), "k-", label="unit circle (exact)")
ax.plot(x_fe, y_fe, "o-", color="tab:red", ms=4, label="forward Euler, h = 0.2")
ax.plot(x_si, y_si, "^-", color="tab:green", ms=4, label="semi-implicit, h = 0.2")
ax.plot(1, 0, "k*", ms=14, label="start")            # mark the initial point
ax.set_aspect("equal")                                # circles must look round
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.legend()
plt.show()
```

- **Plot $y$ against $x$**, not against $t$: `ax.plot(x, y)`.
- **`ax.set_aspect("equal")`**: without it, the axes have different scales and a circle looks like an ellipse, which hides exactly the error you want to show.
- **Format strings** `"o-"`, `"s-"`, `"^-"`, `"x-"`, `"D-"`: a marker for every computed point, joined by lines. Different colours **and** markers keep curves apart even in black-and-white print. Use the same colour/marker for the same method in all figures.
- **One `plt.subplots` + `plt.show()` per required plot** gives separate figures (1c (i), (ii), (iii)).
- **Legend placement:** if it covers data, move it below the axes with `ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2)`.
- **Jupyter submission:** *Kernel → Restart & Run All*, then save, so that all outputs and plots are stored in the `.ipynb` file.

# Self-check before you start

1. How do you show that $E = x^2 + y^2$ is conserved without solving the ODE? *(Compute $dE/dt$ with the chain rule, insert the ODE, get 0.)*
2. What curve in the $(x, y)$ plane does a trajectory with constant $x^2 + y^2$ follow? *(A circle around the origin.)*
3. In semi-implicit Euler, which value of the first-updated variable does the second line use? *(The new one, $n+1$.)*
4. How many evaluations of $f$ does one RK2 midpoint step need, and what is its order? *(Two; second order.)*
5. For $\dot x = y$, $\dot y = -x$, by which factor does forward Euler change $E$ in one step? *($1 + h^2$.)*
6. Backward Euler spirals inwards. Is it stable? Is it energy-conserving? *(Stable: yes. Energy-conserving: no.)*
