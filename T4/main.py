import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('Qt5Agg')  # Указываем интерактивный бэкенд

def get_data(filename, time):
   with open("output.idx", 'r') as f:
       for t, line in enumerate(f):
           if t == time+1:
               meta = tuple(map(float, line.split()))
                
   with open(filename, 'r') as file:
       file.seek(meta[2])
       matrix = np.fromfile(file, dtype=np.float64, count=np.prod(size[1:])).reshape(size[1:])
        
   return (meta[1], np.array(matrix))

# Функция для обновления графика при изменении времени
def update(val):
   time = int(time_slider.val)
       
   # Фильтруем данные для выбранного времени
   (time_val, data) = get_data("output.bin", time)
    
   # Обновляем данные для графика
   x = np.linspace(0, 1, size[1])
   y = np.linspace(0, 1, size[2])
   X, Y = np.meshgrid(x, y)
    
   # Очищаем предыдущий график и строим новый
   ax.clear()
   scatter = ax.plot_surface(X, Y, data, cmap='viridis', vmin=get_data("output.bin", 0)[1].min(), vmax=get_data("output.bin", 0)[1].max())
    
   # Добавляем colorbar
   if not hasattr(update, 'colorbar'):
       update.colorbar = plt.colorbar(scatter, ax=ax)
   else:
       update.colorbar.update_normal(scatter)
    
   # Устанавливаем подписи осей
   
   ax.set_xlabel('X')
   ax.set_ylabel('Y')
   ax.set_zlabel('Z')
   ax.set_zlim(get_data("output.bin", 0)[1].min(), get_data("output.bin", 0)[1].max())
   ax.set_title(f'3D график данных при времени = {time_val:.2f}')    
   fig.canvas.draw_idle()

# Функция для обработки нажатий клавиш
def on_key(event):
    if event.key == 'right' or event.key == 'up':
        # Увеличиваем значение слайдера на 1
        new_val = min(time_slider.val + 1, time_slider.valmax)
        time_slider.set_val(new_val)
    elif event.key == 'left' or event.key == 'down':
        # Уменьшаем значение слайдера на 1
        new_val = max(time_slider.val - 1, time_slider.valmin)
        time_slider.set_val(new_val)

# Считываем размеры из первой строки индексного файла
with open("output.idx", 'r') as f:
   size = f.readline().strip().split()
   size = tuple(map(int, size))

# Создаем фигуру и 3D оси
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Добавляем ползунок для времени
ax_time = plt.axes([0.25, 0.02, 0.65, 0.03])


time_slider = Slider(
   ax=ax_time,
   label='Время',
   valmin=0,
   valmax=size[0]-1,
   valinit=0,
   valstep=1
)

# Регистрируем функцию обновления при изменении ползунка
time_slider.on_changed(update)

# Регистрируем обработчик событий клавиатуры
fig.canvas.mpl_connect('key_press_event', on_key)

# Инициализируем график с начальным значением времени
update(0)

# Вместо tight_layout используем более гибкую настройку
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

plt.show()