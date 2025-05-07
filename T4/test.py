import numpy as np 

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

x = np.linspace(0, 1, 2)
y = np.linspace(0, 1, 2)

# Очищаем предыдущий график и строим новый
ax.clear()
scatter = ax.plot_surface(x, y, [[0, 0.5], [0.5, 2]])

# Устанавливаем подписи осей
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title(f'3D график данных при времени = {1/3:.2f}')


fig.canvas.draw_idle()

plt.tight_layout()
plt.savefig("out")