import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(-4, 4, 200)
fig, ax = plt.subplots()

ax.grid(True);
ax.axis([-4, 4, -4, 4])
C=[0, 0.5, -0.5, 1, -1, 2, -2]
for c in C:
    ax.plot(x, (x-c)**3)
y = np.linspace(-4, 4, 200)
ax.plot(x, (x**2)/4+(y**2)/9 -1)
plt.show()
