---
title: "Modeling and Simulation: Concept Sheet for Problem Sheet 6"
subtitle: "Electric circuits, systems of ODEs, implicit methods and stability (WS 2026/27)"
author: "Subhash"
---

This sheet adds only the **new** ideas you need for problem sheet 6. Ideas that are already explained in the week 5 concept sheet (`week 5/ps_5_concepts.md`) are **not repeated**. Here is where to find them:

| Already covered in week 5 | Section there |
|:--|:-:|
| What an ODE is, and how to classify one | 1, 2 |
| Converting an $n$-th order ODE to a first-order system | 4 |
| The five-step modelling and simulation workflow | 5 |
| Forward (explicit) Euler: derivation, meaning, $\mathcal O(h)$ error | 7 |
| Programming pattern: time grid, arrays, initial condition, loop | 9 |
| Plausibility checks, convergence checks, comparison with an exact solution | 10 |

**New in this sheet:**

| Exercise part | New concepts | Section |
|:--|:--|:-:|
| (a) abstraction | circuit elements, Kirchhoff's voltage law, sign conventions, state variables | 1 |
| (b) model | a circuit as a system of two first-order ODEs, matrix form | 2 |
| (c) discretization | Euler methods for **systems**, backward Euler, trapezoidal method | 3, 4, 5 |
| (e) interpretation | stability, comparing the three methods, damped oscillation, energy check, verification vs. validation | 6, 7, 8, 9 |
| (d) plotting | `plt.subplots` with two graphs | 10 |

The outline on slide 26 lists *stability*, *backward Euler* and the *trapezoidal method* (sections 9.3 to 9.5), but the slides provided stop at slide 34, before those sections. The explanations in sections 4 to 6 below are the standard textbook versions of these methods.

# 1. Electric circuits: what you need from slide 8

## 1.1 The three passive elements

| Element | Law (slide 8) | In words |
|:--|:--|:--|
| Resistor $R$ (Ω) | $v_R = R\,i$ | voltage is proportional to the current |
| Inductor $L$ (H) | $v_L = L\,\dfrac{di}{dt}$ | voltage is proportional to how fast the current **changes** |
| Capacitor $C$ (F) | $i = C\,\dfrac{dv_C}{dt}$ | current is proportional to how fast the voltage **changes** |

"Passive" means the elements only store or use up energy; none of them supplies it. Without a source, all the energy in the circuit is the energy stored at $t = 0$.

## 1.2 Kirchhoff's laws (slide 8)

- **Voltage law (KVL):** around any closed loop, the voltages add up to zero, $\sum_l v_l = 0$.
- **Current law (KCL):** at any node, the currents add up to zero, $\sum_l i_l = 0$.

In a **series** circuit there is only one loop and no branching, so the **same current** $i$ flows through every element. KCL is then automatic, and KVL gives the one equation you need.

## 1.3 Sign conventions: the most common source of errors

1. **Choose a direction for the loop current** $i$ and draw it as an arrow. The choice is free. If the real current flows the other way, $i$ simply comes out negative.
2. **Mark the polarity** ($+$ and $-$) of each voltage you name, especially the output voltage $v$.
3. **Passive sign convention:** the laws in the table hold in that form only if the current enters the element at its $+$ terminal. If you drew the current entering at the $-$ terminal, add a minus sign, e.g. $i = -C\,\dfrac{dv}{dt}$.
4. **Walk around the loop** in the direction of $i$ and add up the voltage drops. With the passive sign convention every element is a drop, and KVL becomes "sum of drops $= 0$".

**Worked example (RC circuit).** A charged capacitor $C$ discharges through a resistor $R$. Draw $i$ so that it enters the capacitor at $+$.

- KVL: $v_R + v = 0$, so $R\,i + v = 0$.
- Capacitor law: $i = C\,\dfrac{dv}{dt}$.
- Combine: $RC\,\dfrac{dv}{dt} = -v$, so $\dfrac{dv}{dt} = -\dfrac{v}{RC}$.
- Solution: $v(t) = V_0\,e^{-t/(RC)}$. The time constant is $\tau = RC$.

## 1.4 State variables and initial conditions

The **state variables** of a circuit are the quantities that store energy and therefore cannot jump:

- the **capacitor voltage**, because the stored energy is $\tfrac12 C v^2$;
- the **inductor current**, because the stored energy is $\tfrac12 L i^2$.

These are exactly the quantities that get initial values: $v(0) = V_0$ and $i(0) = I_0$. One energy-storing element gives one first-order ODE, so a circuit with $L$ and $C$ gives **two**.

# 2. From a circuit to a system of first-order ODEs

**Recipe.**

1. Write KVL for the loop.
2. Replace each voltage using the element laws.
3. Solve for the derivative of each state variable.

**Matrix form.** A linear system without a source can be written as
$$
\frac{d\vec x}{dt} = A\,\vec x, \qquad \vec x(0) = \vec x_0 ,
$$
where $\vec x$ holds the state variables and $A$ is a constant matrix. Writing it like this makes every method in sections 3 to 5 a short formula.

**Worked example (LC circuit, slide 8).** Take an inductor and a capacitor in a loop, without a resistor. With $v$ across $C$ and $i$ entering at $+$:

- KVL: $L\,\dfrac{di}{dt} + v = 0$
- Capacitor: $i = C\,\dfrac{dv}{dt}$

$$
\frac{d}{dt}\begin{pmatrix} v \\ i\end{pmatrix} =
\begin{pmatrix} 0 & 1/C \\ -1/L & 0\end{pmatrix}
\begin{pmatrix} v \\ i\end{pmatrix}.
$$

Eliminating $v$ gives the second-order form on slide 8, $\ddot i + \frac{1}{LC}\,i = 0$: an undamped oscillation with angular frequency $\omega_0 = 1/\sqrt{LC}$. The week 5 order-reduction recipe goes the other way, from this second-order ODE back to the system.

# 3. Euler methods for a system of ODEs

Every method from week 5 works for systems: apply the same formula to **each component** at the same time.

**The important rule:** within one step, a method must use the values it promises to use. Forward Euler uses only **old** values on the right-hand side, for **every** component.

**Example.** Forward Euler for the LC circuit:
$$
v_{n+1} = v_n + \frac{h}{C}\,i_n, \qquad
i_{n+1} = i_n - \frac{h}{L}\,v_n .
$$

If you accidentally use the **new** $v_{n+1}$ in the second line, you get a different method (Euler-Cromer, slide 26, section 9.6). It is a valid method too, but not forward Euler.

In matrix form forward Euler is
$$
\vec x_{n+1} = \vec x_n + h\,A\,\vec x_n = (I + hA)\,\vec x_n .
$$

# 4. The backward (implicit) Euler method

## 4.1 Idea and derivation

Forward Euler steps along the slope at the **start** of the step. Backward Euler steps along the slope at the **end** of the step. It comes from the **backward** Taylor expansion (slide 28) of $y(t_{n+1} - h)$ around $t_{n+1}$:
$$
y(t_n) = y(t_{n+1}) - h\,y'(t_{n+1}) + \mathcal O(h^2)
\quad\Longrightarrow\quad
\boxed{\;y_{n+1} = y_n + h\,f(t_{n+1},\,y_{n+1})\;}
$$

Equivalently, the backward difference (slide 29) is used for $y'(t_{n+1})$: $\dfrac{y_{n+1} - y_n}{h} = f(t_{n+1}, y_{n+1})$.

## 4.2 What "implicit" means

The unknown $y_{n+1}$ appears on **both sides** of the equation. You cannot just evaluate the right-hand side; you have to **solve an equation** in every step.

- For a **linear** ODE this is a linear equation, which you can solve exactly.
- For a non-linear ODE you need an iterative solver such as Newton's method.

For the linear system $\vec x' = A\vec x$:
$$
\vec x_{n+1} = \vec x_n + h\,A\,\vec x_{n+1}
\quad\Longrightarrow\quad
(I - hA)\,\vec x_{n+1} = \vec x_n .
$$

**Two ways to solve this in code:**

1. **By hand (best for a 2×2 system):** substitute one equation into the other to get an explicit formula for each new value. Then the loop needs only basic arithmetic.
2. **With NumPy:** `x_new = np.linalg.solve(np.eye(2) - h*A, x_old)`.

**Worked example by hand.** $y' = -2y$, $y(0)=1$, $h = 0.1$:
$$
y_{n+1} = y_n + 0.1\cdot(-2\,y_{n+1})
\;\Longrightarrow\; 1.2\,y_{n+1} = y_n
\;\Longrightarrow\; y_{n+1} = \frac{y_n}{1.2}.
$$

| $t$ | forward Euler | backward Euler | exact $e^{-2t}$ |
|:-:|:-:|:-:|:-:|
| 0.1 | 0.8000 | 0.8333 | 0.8187 |
| 0.2 | 0.6400 | 0.6944 | 0.6703 |
| 0.3 | 0.5120 | 0.5787 | 0.5488 |

Forward Euler falls too fast and backward Euler too slowly, by about the same amount. Both are first order: local error $\mathcal O(h^2)$, global error $\mathcal O(h)$.

# 5. The trapezoidal (implicit) method

## 5.1 Idea

Use the **average** of the slope at the start and the slope at the end of the step:
$$
\boxed{\;y_{n+1} = y_n + \frac{h}{2}\Big[f(t_n, y_n) + f(t_{n+1}, y_{n+1})\Big]\;}
$$

**Why the name.** Integrate the ODE over one step (slide 27):
$$
y(t_{n+1}) = y(t_n) + \int_{t_n}^{t_{n+1}} f\,dt .
$$
Approximating the integral by the area of a **trapezoid** gives exactly this formula. Forward and backward Euler use a rectangle instead, with the height taken on the left or on the right.

The method is **implicit** because $y_{n+1}$ appears on the right. The errors of forward and backward Euler have opposite signs and largely cancel, so the trapezoidal method is **second order**: local error $\mathcal O(h^3)$, global error $\mathcal O(h^2)$. This is the same cancellation as in the central difference in week 4.

For $\vec x' = A\vec x$:
$$
\left(I - \frac h2 A\right)\vec x_{n+1} = \left(I + \frac h2 A\right)\vec x_n .
$$

**Worked example by hand.** $y' = -2y$, $h = 0.1$:
$$
y_{n+1} = y_n + 0.05\,(-2y_n - 2y_{n+1})
\;\Longrightarrow\; 1.1\,y_{n+1} = 0.9\,y_n
\;\Longrightarrow\; y_{n+1} = \tfrac{9}{11}\,y_n .
$$

| $t$ | trapezoidal | exact | error |
|:-:|:-:|:-:|:-:|
| 0.1 | 0.8182 | 0.8187 | 0.0005 |
| 0.2 | 0.6694 | 0.6703 | 0.0009 |
| 0.3 | 0.5477 | 0.5488 | 0.0011 |

The error is about 30 times smaller than for the Euler methods, at the same step size.

## 5.2 Summary of the three methods

| Method | Update | Type | Global error |
|:--|:--|:-:|:-:|
| Forward Euler | $y_{n+1} = y_n + h\,f_n$ | explicit | $\mathcal O(h)$ |
| Backward Euler | $y_{n+1} = y_n + h\,f_{n+1}$ | implicit | $\mathcal O(h)$ |
| Trapezoidal | $y_{n+1} = y_n + \frac h2 (f_n + f_{n+1})$ | implicit | $\mathcal O(h^2)$ |

Here $f_n = f(t_n, y_n)$ and $f_{n+1} = f(t_{n+1}, y_{n+1})$.

# 6. Stability of ODE integrators

## 6.1 Why accuracy is not enough

A method can be accurate for small $h$ and still **blow up** for a larger $h$: the numerical solution grows without limit, although the true solution decays. This is called **numerical instability**. It is not caused by round-off; it is a property of the method and the step size.

## 6.2 The test equation and the amplification factor

To study stability, apply each method to the **test equation**
$$
y' = \lambda\,y, \qquad \text{exact solution } y(t) = y_0\,e^{\lambda t},
$$
where $\lambda$ may be complex. If $\operatorname{Re}\lambda < 0$, the exact solution decays.

Each method turns this into $y_{n+1} = G\,y_n$, so $y_n = G^n\,y_0$. The **amplification factor** $G$ decides everything:

- $|G| < 1$: the numerical solution decays, as it should. The method is **stable**.
- $|G| > 1$: it grows every step and explodes. The method is **unstable**.

With $z = h\lambda$:

| Method | Amplification factor $G(z)$ | Stable when |
|:--|:-:|:--|
| exact | $e^{z}$ | |
| Forward Euler | $1 + z$ | $\lvert 1 + h\lambda\rvert < 1$: only for small $h$ |
| Backward Euler | $\dfrac{1}{1 - z}$ | always, if $\operatorname{Re}\lambda < 0$ |
| Trapezoidal | $\dfrac{1 + z/2}{1 - z/2}$ | always, if $\operatorname{Re}\lambda < 0$ |

Methods that are stable for every $h$ whenever the true solution decays are called **A-stable**. Backward Euler and the trapezoidal method are A-stable. Forward Euler is not. **This is the main reason to accept the extra work of an implicit method.**

## 6.3 Systems: use the eigenvalues

For $\vec x' = A\vec x$ the role of $\lambda$ is played by the **eigenvalues** of $A$. The method is stable only if $|G(h\lambda)| < 1$ for **every** eigenvalue.

For a damped oscillation the eigenvalues are complex, $\lambda = -\alpha \pm i\,\omega_d$. Then
$$
|1 + h\lambda|^2 = (1 - h\alpha)^2 + (h\omega_d)^2 < 1
\quad\Longleftrightarrow\quad
h < \frac{2\alpha}{\alpha^2 + \omega_d^2} .
$$

**Consequence:** a lightly damped oscillator (small $\alpha$) forces forward Euler to use a **very** small step size. For an undamped oscillator ($\alpha = 0$), forward Euler is unstable for **every** $h$: the amplitude grows a little in each step.

**Worked example.** For $y' = -2y$ ($\lambda = -2$), forward Euler needs $|1 - 2h| < 1$, so $h < 1$.

- $h = 1.5$: $G = -2$. The solution goes $1, -2, 4, -8, \dots$ and explodes.
- Backward Euler with $h = 1.5$: $G = 1/4$. The solution goes $1, 0.25, 0.0625, \dots$ It is not accurate, but it is stable and decays as it should.

## 6.4 Too much damping, too little damping

Even when all three methods are stable, they change the **amplitude** of an oscillation differently:

- Forward Euler has $|G|$ slightly **larger** than $|e^{h\lambda}|$: the oscillation decays **too slowly**. With an unstable $h$, it grows instead.
- Backward Euler has $|G|$ **smaller** than $|e^{h\lambda}|$: it adds artificial (numerical) damping, so the oscillation decays **too fast**.
- The trapezoidal method keeps $|G|$ very close to $|e^{h\lambda}|$. For $\operatorname{Re}\lambda = 0$ it gives $|G| = 1$ exactly, so it preserves the amplitude of an undamped oscillation.

# 7. The damped oscillation: what the exact solution looks like

The series RLC circuit leads to the same equation as a mass on a spring with friction (slide 7):
$$
\ddot y + 2\alpha\,\dot y + \omega_0^2\,y = 0 .
$$

Here $\alpha$ is the **damping rate** and $\omega_0$ is the **natural angular frequency**. Comparing coefficients for the series RLC circuit gives $\alpha = \dfrac{R}{2L}$ and $\omega_0 = \dfrac{1}{\sqrt{LC}}$.

| Case | Condition | Behaviour |
|:--|:--|:--|
| underdamped | $\alpha < \omega_0$ | oscillation inside a decaying envelope $e^{-\alpha t}$, with $\omega_d = \sqrt{\omega_0^2 - \alpha^2}$ |
| critically damped | $\alpha = \omega_0$ | fastest return to zero without oscillating |
| overdamped | $\alpha > \omega_0$ | slow decay without oscillating |

**Underdamped solution.** For $y(0) = y_0$ and $\dot y(0) = \dot y_0$:
$$
y(t) = e^{-\alpha t}\left[\,y_0\cos(\omega_d t) + \frac{\alpha y_0 + \dot y_0}{\omega_d}\sin(\omega_d t)\right].
$$

**Useful numbers to read off a plot:**

- period $T_d = 2\pi/\omega_d$, i.e. the time between two maxima;
- after one period the amplitude has shrunk by the factor $e^{-\alpha T_d}$;
- after $t = 1/\alpha$ the envelope has fallen to $e^{-1} \approx 37\,\%$.

The eigenvalues of the system matrix are $\lambda = -\alpha \pm i\,\omega_d$, which links this section to the stability condition in section 6.3.

# 8. Energy as a check quantity

In a passive circuit without a source, the stored energy
$$
E(t) = \tfrac12 C\,v^2 + \tfrac12 L\,i^2
$$
can only **decrease**. Using the circuit equations:
$$
\frac{dE}{dt} = C v\,\dot v + L i\,\dot i = v\,i + i\,(-v - Ri) = -R\,i^2 \le 0 .
$$
The resistor turns the stored energy into heat.

This gives a check that needs no exact solution:

- If $E$ ever **increases** in a simulation, the method or the step size is wrong.
- Comparing $E(t_{end})/E(0)$ between methods shows which one damps too much or too little.

The same idea appears in mechanics: with friction, the total mechanical energy can only decrease.

# 9. Verification, validation and plausibility

These three words mean different things in simulation:

| Term | Question | How you check it |
|:--|:--|:--|
| **Verification** | Are the equations solved **correctly**? (the code and method) | compare with an exact solution; check that the error shrinks with $h$ at the expected order; compare methods with each other |
| **Validation** | Are these the **right** equations? (the model) | compare with measurements or known physics; check the assumptions (ideal elements, no wire resistance, etc.) |
| **Plausibility** | Does the result look physically sensible? | signs, initial values, limits ($R = 0$ gives an undamped oscillation), energy never increases, period $\approx 2\pi\sqrt{LC}$ |

A verified program can still be useless if the model is not valid. A valid model can still give wrong numbers if the program is not verified. You need both.

# 10. Plotting two graphs in one figure

Use `plt.subplots` to put voltage and current in one figure. The rest of the programming pattern is in week 5, section 9.

```python
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)   # 2 rows, 1 column

ax1.plot(t, v_fe, label="forward Euler")
ax1.plot(t, v_be, label="backward Euler")
ax1.set_ylabel("v in V")
ax1.legend()

ax2.plot(t, i_fe, label="forward Euler")
ax2.plot(t, i_be, label="backward Euler")
ax2.set_xlabel("t in s")
ax2.set_ylabel("i in A")
ax2.legend()

plt.tight_layout()
plt.show()
```

`sharex=True` gives both graphs the same time axis. Use the **same colour for the same method** in both graphs, so the reader can connect them.

**Time grid that includes the end point.** `np.arange(0, 1, h)` stops **before** 1. Use `np.arange(0, 1 + h/2, h)` to include $t = 1$ without risking an extra point from rounding (floating-point errors from week 4).

# Self-check before you start

1. In a series circuit, how many different currents are there? *(One.)*
2. Which quantities get initial conditions in an RLC circuit, and why? *(Capacitor voltage and inductor current: they store energy and cannot jump.)*
3. Why must backward Euler solve an equation in each step? *(The new value appears on both sides.)*
4. For $y' = -10y$, what is the largest stable step for forward Euler? *($h < 0.2$.)*
5. Which of the three methods is second order? *(Trapezoidal.)*
6. A simulation of a passive circuit shows the energy increasing. What does that tell you? *(The method is unstable for this step size, or there is a sign error.)*
