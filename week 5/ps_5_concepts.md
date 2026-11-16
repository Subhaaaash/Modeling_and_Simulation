---
title: "Modeling and Simulation: Concept Sheet for Problem Sheet 5"
subtitle: "Ordinary differential equations, forward Euler, and modelling a projectile in air (WS 2026/27)"
author: "Subhash"
---

This sheet collects every idea you need for problem sheet 5, with a short explanation and a small worked example for each. The examples are deliberately **different** from the exercises, so you can practise on them first. Slide numbers refer to the lecture *Part 3: Continuous models: Ordinary differential equations*.

| Exercise | Concepts needed | Section |
|:--|:--|:-:|
| 1 (a) to (e): classify 7 ODEs | order, scalar/vector, linear, autonomous, homogeneous, source term | 1, 2 |
| 1: solve equation (3) | separation of variables, initial conditions | 3 |
| 2: 4th order to 1st order system | order reduction | 4 |
| 3 (a), (b): abstraction, model | modelling workflow, Newton's law in 1D, signs of forces | 5, 6 |
| 3 (c): discretization | forward Euler, forward finite difference | 7, 8 |
| 3 (d): programming | arrays, loop, stopping condition, plot with two y-axes | 9 |
| 3 (e): interpretation | plausibility checks, terminal velocity, analytical check | 10 |

# 1. What an ODE is

A **differential equation** is any equation that contains derivatives (slide 16).

- In an **ordinary** differential equation (ODE) the unknown function depends on **one** independent variable, e.g. $y(x)$ or $x(t)$.
- In a **partial** differential equation (PDE) it depends on several, e.g. the temperature $u(\vec r,t)$ in the heat equation.

**Watch the variable names.** The letters change between equations, so always ask *"which function is unknown, and what is it a function of?"*

- In $\frac{d}{dt}u(t)$ the unknown is $u$ and the independent variable is $t$.
- In $\frac{d^4}{dy^4}x(y)$ the unknown is $x$ and the independent variable is $y$. That is unusual, but allowed.

The general form of an ODE of order $n$ (slide 17) is
$$
F\!\left(x,\,y,\,\frac{dy}{dx},\,\frac{d^2y}{dx^2},\,\dots,\,\frac{d^ny}{dx^n}\right)=0 .
$$

# 2. The five classification questions

For every ODE, answer these five questions in order.

## 2.1 Order

The **order** is the highest derivative that appears (slide 17).

- $y'' + 3y = 0$ has order 2.
- $(y')^5 + y = 0$ has order **1**. The power 5 does not matter, only the number of dashes.
- $y = \sin y$ has order 0: it is an algebraic equation, not an ODE.

## 2.2 Scalar or vectorial

The unknown can be one number (**scalar**) or a vector of several numbers (**vectorial**) (slide 18).

- $y'(t) = -y(t)$ is scalar.
- $\frac{d}{dt}\vec r(t) = \vec v$ is vectorial. An arrow over the unknown, as in $\vec x(t)$, is the signal.

A vectorial ODE is really a set of coupled scalar ODEs, one per component (slide 19).

## 2.3 Linear or non-linear

An ODE of order $n$ is **linear** if and only if it can be written as (slide 24)
$$
\sum_{j=0}^{n} a_j(x)\,\frac{d^j y}{dx^j}+r(x)=0 .
$$

In words: every term is either **(a coefficient depending only on $x$) × (the unknown or one of its derivatives)**, or **a function of $x$ alone** (that is $r(x)$).

**Quick test.** The ODE is non-linear if the unknown $y$ or any derivative of it

- appears with a power other than 1, such as $y^2$, $y^3$ or $(y')^2$;
- appears inside a function, such as $\sin y$, $e^y$ or $\sqrt{y}$;
- is multiplied by itself or by another derivative, such as $y\,y'$.

The coefficients $a_j(x)$ may be **anything** in the independent variable, e.g. $x^2$, $1/x$ or $e^{-x}$. They do **not** make the ODE non-linear.

| ODE | Linear? | Reason |
|:--|:-:|:--|
| $x^2y'' + e^{x}y = 0$ | yes | coefficients depend only on $x$ |
| $y' = y^2$ | no | $y^2$ |
| $y'' + \sin(y) = 0$ | no | $y$ inside $\sin$ |
| $y'' + \sin(x)\,y = 0$ | yes | $\sin$ acts on the independent variable |
| $y\,y' = 1$ | no | product of $y$ and $y'$ |

## 2.4 Autonomous or non-autonomous

An ODE is **autonomous** if the independent variable does **not** appear explicitly (slide 23). It may only enter through the unknown function itself.

- $y'' + ky = 0$ is autonomous.
- $y'' + \sin(x) = 0$ and $x\,y' + c = 0$ are non-autonomous, because $x$ appears on its own.

**Physical meaning.** An autonomous system follows the same law at every moment. If you start the experiment one hour later, you get the same curve shifted by one hour.

Linear and autonomous are **independent** questions. $y' = t\,y$ is linear but non-autonomous. $y' = y^2$ is non-linear but autonomous.

## 2.5 Homogeneous or inhomogeneous, and the source term

Write the ODE in the form of slide 24 with **everything on one side**. The **source term** $r(x)$ is the part that contains **no** unknown at all.

- $r(x) = 0$: **homogeneous**.
- $r(x) \neq 0$: **inhomogeneous**, and $r(x)$ is the source term (also called the inhomogeneity).

**Example.** Take $y'' - k y + \sin(x) = 0$. Every term except $\sin(x)$ contains $y$, so $r(x)=\sin(x)$ and the ODE is inhomogeneous.

**Two common traps.**

1. **Look at the whole term, not single factors.** In $y'' = x\,y$ the factor $x$ multiplies $y$, so it is a coefficient, not a source. The ODE is homogeneous.
2. **Mind the sign convention.** Slide 24 puts $r(x)$ on the **same** side as the derivatives. For $y'' + y = \cos(x)$, bring everything to the left: $y'' + y - \cos(x) = 0$, so $r(x) = -\cos(x)$. Many books call the right-hand side $\cos(x)$ the source instead. State which convention you use, and the sign is never wrong.

Strictly, slide 24 defines homogeneous only for **linear** ODEs. For a non-linear ODE the same idea is used in practice: "is there a term that does not contain the unknown at all?" When you answer for a non-linear ODE, say that you are using this extended meaning.

## 2.6 Worked classification example

Classify $\;t\,\ddot y + 4\,\dot y^{\,2} = e^{t}$.

| Question | Answer | Why |
|:--|:--|:--|
| order | 2 | highest derivative is $\ddot y$ |
| scalar/vectorial | scalar | $y$ has no arrow |
| linear | non-linear | $\dot y^{\,2}$ |
| autonomous | non-autonomous | $t$ appears in $t\,\ddot y$ and in $e^{t}$ |
| homogeneous | inhomogeneous | $t\ddot y + 4\dot y^2 - e^t = 0$ gives $r(t) = -e^{t}$ |

# 3. Solving a first-order linear homogeneous ODE

The general first-order linear homogeneous ODE is
$$
y'(t) = a(t)\,y(t).
$$

**Method: separation of variables.** Put everything with $y$ on one side and everything with $t$ on the other, then integrate both sides.

1. Separate: $\dfrac{dy}{y} = a(t)\,dt$
2. Integrate: $\ln|y| = \displaystyle\int a(t)\,dt + C$
3. Solve for $y$: $y(t) = K\,e^{\int a(t)\,dt}$, where $K$ is any constant.

An ODE alone has **infinitely many** solutions, one for each $K$. The **initial condition** picks one of them (slide 22).

**Worked example.** Solve $y' = -2t\,y$ with $y(0)=3$.

1. Separate: $\dfrac{dy}{y} = -2t\,dt$
2. Integrate: $\ln|y| = -t^2 + C$
3. General solution: $y(t) = K\,e^{-t^2}$
4. Initial condition: $y(0) = K\cdot e^{0} = K = 3$
5. Particular solution: $y(t) = 3\,e^{-t^2}$

**Always check by substituting.** $y' = 3\,e^{-t^2}\cdot(-2t) = -2t\,y$ $\checkmark$ and $y(0)=3$ $\checkmark$.

# 4. Converting an $n$-th order ODE to a first-order system

Any ODE of order $n$ can be rewritten as $n$ first-order ODEs (slides 20, 21). This matters because numerical methods such as Euler work on first-order equations only.

**Recipe.**

1. Solve the ODE for the highest derivative.
2. Introduce new unknowns: the function and all derivatives **below** the highest one,
$$
y_1 = y,\quad y_2 = y',\quad \dots,\quad y_n = y^{(n-1)} .
$$
3. Write one equation per new unknown. All except the last just say "the derivative of one is the next one". The last one is the original ODE.

**Worked example.** Convert $y''' + 2y'' - y = \cos(t)$.

1. Highest derivative: $y''' = -2y'' + y + \cos(t)$
2. New unknowns: $y_1 = y$, $y_2 = y'$, $y_3 = y''$
3. System:
$$
\begin{aligned}
y_1' &= y_2\\
y_2' &= y_3\\
y_3' &= -2y_3 + y_1 + \cos(t)
\end{aligned}
\qquad\Longleftrightarrow\qquad
\frac{d}{dt}\begin{pmatrix}y_1\\y_2\\y_3\end{pmatrix}=
\begin{pmatrix}y_2\\y_3\\-2y_3+y_1+\cos t\end{pmatrix}.
$$

A third-order ODE gives **three** first-order equations. To solve it you need **three** initial conditions: $y(0)$, $y'(0)$ and $y''(0)$.

# 5. The modelling and simulation workflow

The first lecture introduced a fixed workflow. Exercise 3 asks you to follow it step by step.

| Step | Question it answers | What you write down |
|:--|:--|:--|
| 1. Abstraction | What matters, and what can I ignore? | sketch with forces, list of quantities and parameters, assumptions, goal |
| 2. Mathematical model | Which equations describe it? | ODEs plus initial conditions |
| 3. Discretization | How does a computer step through time? | update formulas such as $v_{i+1} = \dots$ |
| 4. Programming | How do I compute it? | code |
| 5. Interpretation | Is the answer correct and plausible? | answers, checks, conclusions |

**Typical assumptions in the abstraction step:** the projectile is a point mass, $g$ is constant, air density is constant, there is no wind, and the motion is purely vertical.

# 6. Newton's law in 1D and the signs of forces

Newton's second law in one dimension, with the $z$-axis pointing **up**, is
$$
m\,\frac{dv}{dt} = \sum F_z ,
$$
where $F_z$ is **positive** if the force points up and **negative** if it points down.

**How to get the signs right.**

1. Draw the axis first.
2. Draw each force as an arrow.
3. Give each force the sign of its arrow relative to the axis.

| Force | Direction | Sign in the ODE |
|:--|:--|:--|
| Gravity $mg$ | always down | always $-mg$ |
| Air drag $kv^2$ | always **opposite to the velocity** | $-kv^2$ when moving up, $+kv^2$ when moving down |

Drag is a friction force: it always fights the motion. $v^2$ alone is always positive, so it cannot carry the direction. You have to add the sign yourself from the direction of motion (one compact form is $-k\,v|v|$, which works for both directions).

The drag constant is $k = \tfrac12\,c_w\,\rho\,A$, with unit kg/m. Check: $k v^2$ has unit $\tfrac{\text{kg}}{\text{m}}\cdot\tfrac{\text{m}^2}{\text{s}^2} = \text{kg}\cdot\tfrac{\text{m}}{\text{s}^2} = \text{N}$ $\checkmark$.

**Velocity and position.** The velocity is the rate of change of the position:
$$
\frac{dz}{dt} = v .
$$

**Worked example: falling with linear drag.** A ball of mass $m$ is **dropped** from rest. It moves down, so $v < 0$. Use linear drag $F = -b\,v$, which already has the right sign: $v<0$ gives a force that points up.

- Model: $m\dfrac{dv}{dt} = -mg - b\,v$ with $v(0) = 0$
- At the start, $v = 0$, so $\dfrac{dv}{dt} = -g$: it falls freely at first.
- **Terminal velocity:** speed stops changing when the forces balance, $\dfrac{dv}{dt}=0$. This gives $v_\infty = -\dfrac{mg}{b}$.

For quadratic drag, the same balance $mg = kv^2$ gives the terminal speed
$$
v_t = \sqrt{\frac{mg}{k}} .
$$
This is a useful **speed scale**. When $|v| \gg v_t$, drag dominates. When $|v| \ll v_t$, gravity dominates.

# 7. The forward (explicit) Euler method

Most ODEs cannot be solved by hand (slide 25), so we approximate them on a grid of time points
$$
t_i = t_0 + i\,h,\qquad i = 0, 1, 2, \dots
$$
with step size $h$, and write $y_i \approx y(t_i)$.

**Derivation (slides 31, 32).** Take the Taylor expansion of $y(t+h)$ and drop the $h^2$ term:
$$
y(t+h) = y(t) + h\,y'(t) + \frac{h^2}{2}y''(\tau) \;\approx\; y(t) + h\,y'(t).
$$
The ODE tells us the slope, $y'(t_i) = f(t_i, y_i)$. This gives the **forward Euler scheme**:
$$
\boxed{\;y_{i+1} = y_i + h\,f(t_i, y_i),\qquad y_0 = y(t_0).\;}
$$

**Meaning (slide 33).** At each point, look at the slope the ODE gives you, and walk a step $h$ in a straight line along that slope. Then repeat from the new point.

**Errors.**

- The dropped term gives a **local** error of order $\mathcal O(h^2)$ in each step.
- After $T/h$ steps these add up to a **global** error of order $\mathcal O(h)$ (slide 30).
- So halving $h$ roughly halves the error. This is the same first-order behaviour as the forward difference in problem sheet 4.

**Worked example by hand.** Solve $y' = -2y$, $y(0)=1$ with $h = 0.1$. Here $f(t,y) = -2y$, so the scheme is $y_{i+1} = y_i + 0.1\cdot(-2y_i) = 0.8\,y_i$.

| $i$ | $t_i$ | Euler $y_i$ | exact $e^{-2t_i}$ | error |
|:-:|:-:|:-:|:-:|:-:|
| 0 | 0.0 | 1.0000 | 1.0000 | 0 |
| 1 | 0.1 | 0.8000 | 0.8187 | 0.019 |
| 2 | 0.2 | 0.6400 | 0.6703 | 0.030 |
| 3 | 0.3 | 0.5120 | 0.5488 | 0.037 |

Euler always uses the slope at the **start** of the step. Here the curve flattens as it falls, so Euler falls a bit too fast.

# 8. Forward finite difference for the position

Once the velocities $v_i$ are known, the position follows from $dz/dt = v$. Replace the derivative with the forward difference (slide 29):
$$
\frac{z_{i+1} - z_i}{h} \approx v_i
\quad\Longrightarrow\quad
\boxed{\;z_{i+1} = z_i + h\,v_i,\qquad z_0 = z(0).\;}
$$

This is just forward Euler applied to $dz/dt = v$. Geometrically, each step adds the area of a thin rectangle $h \cdot v_i$ under the velocity curve, so the position is the running area under $v(t)$.

**Worked example.** With $v_0 = 10$, $v_1 = 8$, $v_2 = 6$ m/s, $h = 0.5$ s and $z_0 = 0$:
$$
z_1 = 0 + 0.5\cdot10 = 5\text{ m},\qquad z_2 = 5 + 0.5\cdot 8 = 9\text{ m},\qquad z_3 = 9 + 0.5\cdot 6 = 12\text{ m}.
$$

# 9. Programming pattern

A forward Euler program always has the same five parts. Slide 34 shows the basic version.

```python
import numpy as np
import matplotlib.pyplot as plt

h = 0.01                       # 1. step size
t = np.arange(0, 5 + h, h)     # 2. time grid
y = np.zeros(len(t))           #    array for the solution
y[0] = 1.0                     # 3. initial condition

for i in range(len(t) - 1):    # 4. Euler loop: new value from old value
    y[i + 1] = y[i] + h * (-2 * y[i])

plt.plot(t, y, label="Euler")  # 5. plot with labels, title, legend
plt.xlabel("t"); plt.ylabel("y"); plt.title("y' = -2y"); plt.legend()
plt.show()
```

**Stopping at an event.** When the end time is not known in advance, e.g. "until the highest point", use a `while` loop that stops when the condition is reached, such as `while v[i] > 0:`. Make the arrays long enough beforehand, then cut them to the computed length.

**Two curves with different units in one figure.** Velocity (m/s) and position (m) need two y-axes, as in Figure 1 of the sheet. Use `ax2 = ax1.twinx()`, and colour each axis like its curve so the reader knows which axis belongs to which curve.

# 10. Interpretation and validation

A simulation always produces numbers, even if they are wrong. Before you trust them, run these checks.

1. **Signs and directions.** Does the velocity decrease while the object rises? Both forces point down, so it must.
2. **Limits.** Without air ($k = 0$) the model must reduce to $v = v_0 - g t$ and $z_{max} = v_0^2/(2g)$. With air, the height must be **smaller**.
3. **Scales.** Compare speeds with the terminal speed $v_t = \sqrt{mg/k}$.
4. **Convergence.** Run again with $h/2$ or $h/10$. The results should change only a little, by roughly a factor 2 or 10 less each time, for a first-order method.
5. **Analytical solution.** If one exists, compare against it.

**Analytical solution for upward motion with quadratic drag.** For $\dfrac{dv}{dt} = -g - \dfrac{k}{m}v^2$, separation of variables gives
$$
\int\frac{dv}{g + \frac{k}{m}v^2} = -\int dt
\quad\Longrightarrow\quad
v(t) = v_t \tan\!\left(\arctan\frac{v_0}{v_t} - \frac{g\,t}{v_t}\right),\qquad v_t = \sqrt{\frac{mg}{k}} .
$$
Setting $v = 0$ gives the time and height of the highest point:
$$
t_{top} = \frac{v_t}{g}\arctan\frac{v_0}{v_t},\qquad
z_{max} = \frac{m}{2k}\ln\!\left(1 + \frac{k\,v_0^2}{m\,g}\right).
$$
You don't need these to solve the exercise, but they give an exact answer to check your simulation against.

# Self-check before you start

1. Is $y' = y\cos(t)$ linear? Autonomous? *(Linear; non-autonomous.)*
2. What is the source term of $y'' + 3y' = 7$? *(With slide 24's convention, $r = -7$; the right-hand side is $7$.)*
3. How many first-order equations does a 5th-order ODE become? *(Five.)*
4. A stone moves **down** with quadratic drag. Which sign does the drag force have, with $z$ pointing up? *(Positive: it points up, against the motion.)*
5. If you halve $h$ in forward Euler, what happens to the global error? *(It roughly halves.)*
