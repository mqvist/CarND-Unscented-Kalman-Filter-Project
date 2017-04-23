import numpy as np

data = []
for line in open('data/obj_pose-laser-radar-synthetic-input.txt'):
    tokens = line.split()
    if tokens[0] == 'L':
        t, px, py, vx, vy, yaw, yaw_rate = tokens[3:]
    else:
        t, px, py, vx, vy, yaw, yaw_rate = tokens[4:]
    data.append((t, px, py, vx, vy, yaw, yaw_rate))

data = np.array(data, dtype='double')
data[:,0] /= 1e6
dt = data[1:,0] - data[:-1,0]
dvx = data[1:,3] - data[:-1,3]
dvy = data[1:,4] - data[:-1,4]
ax = dvx / dt
ay = dvy / dt
print(abs(ax).max())
print(abs(ay).max())
print('ax std dev:', ax.std())
print('ay std dev:', ay.std())
