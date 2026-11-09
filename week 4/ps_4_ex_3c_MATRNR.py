"""Problem sheet 4, exercise 3(c).

Forward and central finite differences of
    f(x) = (2 exp(ix) - 3 exp(-ix)) / (-3i)
at x = 1.0 for h = 2**(-n), n = 1...70, with relative errors, h and h**2
on one logarithmic plot.
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return (2 * np.exp(1j * x) - 3 * np.exp(-1j * x)) / (-3j)


def df_exact(x):
    # d/dx: (2i exp(ix) + 3i exp(-ix)) / (-3i) = -(2 exp(ix) + 3 exp(-ix)) / 3
    return -(2 * np.exp(1j * x) + 3 * np.exp(-1j * x)) / 3


x0 = 1.0
n = np.arange(1, 71)
h = 2.0 ** (-n)

forward = (f(x0 + h) - f(x0)) / h
central = (f(x0 + h) - f(x0 - h)) / (2 * h)

exact = df_exact(x0)
err_forward = np.abs(forward - exact) / np.abs(exact)
err_central = np.abs(central - exact) / np.abs(exact)

best_f = np.argmin(err_forward)
best_c = np.argmin(err_central)

print(f"exact f'(1.0) = {exact:.15f}")
print(f"forward: minimum relative error {err_forward[best_f]:.3e} at n = {n[best_f]}")
print(f"central: minimum relative error {err_central[best_c]:.3e} at n = {n[best_c]}")
print(f"\n{'n':>3} {'h':>11} {'forward':>11} {'central':>11}")
for k in (1, 2, 4, 8, 12, 16, 18, 20, 24, 26, 28, 30, 36, 40, 44, 48, 50,
          51, 52, 53, 54, 55, 60, 70):
    i = k - 1
    print(f"{k:>3} {h[i]:>11.3e} {err_forward[i]:>11.3e} {err_central[i]:>11.3e}")

fig, ax = plt.subplots(figsize=(9, 5.6), dpi=200)
ax.semilogy(n, err_forward, color="#2a78d6", lw=2, label="relative error, forward difference")
ax.semilogy(n, err_central, color="#eb6834", lw=2, label="relative error, central difference")
ax.semilogy(n, h, color="#1baf7a", lw=2, ls="--", label=r"$h = 2^{-n}$")
ax.semilogy(n, h ** 2, color="#eda100", lw=2, ls="--", label=r"$h^2 = 2^{-2n}$")

ax.annotate(f"forward minimum\nn = {n[best_f]}, error {err_forward[best_f]:.1e}",
            xy=(n[best_f], err_forward[best_f]), xytext=(n[best_f] + 12, 1e-18),
            fontsize=9, color="#0b0b0b",
            arrowprops=dict(arrowstyle="-", color="#52514e", lw=0.8))
ax.annotate(f"central minimum\nn = {n[best_c]}, error {err_central[best_c]:.1e}",
            xy=(n[best_c], err_central[best_c]), xytext=(n[best_c] - 15, 1e-27),
            fontsize=9, color="#0b0b0b",
            arrowprops=dict(arrowstyle="-", color="#52514e", lw=0.8))

ax.set_xlabel("n  (step size h = $2^{-n}$)")
ax.set_ylabel("relative error  /  h  /  $h^2$")
ax.set_title("Finite-difference derivative of f(x) at x = 1.0: truncation vs. round-off")
ax.set_xlim(1, 70)
ax.grid(True, which="major", color="#e4e3df", lw=0.7)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.legend(frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(Path(__file__).with_name("ps_4_ex_3c_plot.png"))
plt.show()
