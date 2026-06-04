import numpy as np
from scipy.linalg import norm, svd, lstsq
import matplotlib.pyplot as plt

# The number
G = 6.67e-11

# Asteroid Masses and Coordinates (kg) [position]
Asteroid_Mass = np.array([4.67e+11, 1.5e+10, 1.4e+10])
Asteroid_Location = np.array([ [0., 0., 0.], [0.0, 3110, 0.0], [240, 3110, 0.0] ])
Asteroid_Size = np.cbrt( (Asteroid_Mass / 3500) * 3 / (4 * np.pi) )

# Probe Coordinates [position, velocity]
# Probe_Coords = np.array([430e+3,  -1e+6,  0., 0., 4500., 0.])
Probe_Coords = np.array([430e+1,  0.,  0., 0., 0.09, 0.02])

# Probe Mass (kg)
Probe_Mass = 1

def Kinematics(t, s):
    # get position and velocity
    x = s[0:3]
    v = s[3:]

    # get R vectors
    R = Asteroid_Location - x
    R_norm = norm(R, axis = 1)

    # find the derivative of position
    x_dot = v[:]

    # find the derivative of velocity
    v_dot = [0.0, 0.0, 0.0]
    for i in range( len(Asteroid_Mass) ):
        a_temp = G * Asteroid_Mass[i] * ( R[i] / (R_norm[i])**3 )
        v_dot += a_temp

    #print(v_dot)

    # return the coordinate gradient
    return np.concatenate( (x_dot, v_dot) )

t0, tf = 0, 5000000
dt = 100 # pick your favourite number it is dynamically adjusted
# BUG IF dt TOO LARGE
t = t0
state0 = Probe_Coords[:]
t_vals = []
x_vals = []
v_vals = []

# Verlet Integration

# Relative, Absolute, and Petty Permitted Errors in Energy
r_dEnergy = 1e-6
a_dEnergy = 0
p_dEnergy = 1e-6

# Initial Energy
KE0 = 0.5 * Probe_Mass * norm(state0[3:6]) * norm(state0[3:6])
PE0 = 0
for i in range( len(Asteroid_Mass) ):
    PE0 -= G * Probe_Mass * Asteroid_Mass[i] / norm(state0[0:3] - Asteroid_Location[i])

E = KE0 + PE0
print("Total Energy:", E)

# If we shrink the step size, use this flag to alert the code to not grow it again
drop_flag = False

while t < tf:
    print(t)
    t += dt
    # Find the current acceleration
    state_temp = Kinematics(t, state0) 

    # Step position
    x_temp = state0[0:3] + state0[3:6] * dt + 0.5 * state_temp[3:6] * dt * dt

    # Find the next acceleration
    a_temp = Kinematics(t, np.concatenate( (x_temp, [0., 0., 0.]) ) )[3:6]

    # Step the velocity
    v_temp = state0[3:6] + 0.5 * (state_temp[3:6] + a_temp) * dt

    KE = 0.5 * Probe_Mass * norm(v_temp) * norm(v_temp)
    PE = 0
    for i in range( len(Asteroid_Mass) ):
        PE -= G * Probe_Mass * Asteroid_Mass[i] / norm(x_temp - Asteroid_Location[i])

    # If the energy drift is too large, dynamically adjust timestep
    if np.abs(KE + PE - E) > (a_dEnergy + r_dEnergy * np.abs(E)):
        t -= dt
        dt /= 2
        drop_flag = True
        continue

    # If the energy drift is sufficiently small, dynamically adjust timestep
    if np.abs(KE + PE - E) < (a_dEnergy + r_dEnergy * np.abs(E)) * p_dEnergy and drop_flag == False:
        t -= dt
        dt *= 2
        continue
    
    # Store new position, velocity, and energy
    state0[0:3] = x_temp
    state0[3:6] = v_temp

    x_vals.append(x_temp)
    v_vals.append(v_temp)
    t_vals.append(t)

    #E = KE + PE
    drop_flag = False

v_vectors = np.array(v_vals)
x_vectors = np.array(x_vals)

# Compute circle of best fit for the hodograph
x = v_vectors[:, 0]
y = v_vectors[:, 1]

# Fit the hodograph circle using a monte-carlo method
def fit_circle_3d(points):
    # Mean-center points
    centroid = np.mean(points, axis=0)
    centered = points - centroid

    # SVD to find the best fitting plane
    try:
        U, S, V = svd(centered)
        normal = V[2, :] # Normal is the last row of V (smallest singular value)
    # If scipy says too large of a matrix to do SVD, sample a random set of 2^16
    except:
        print("augh!")
        half_size = len(centered) // 2
        half_indeces = np.random.choice(len(centered), size=2**8, replace=False)
        U, S, V = svd(centered[half_indeces])
        normal = V[2, :]


    # Project 3D points to 2D plane
    # Create two orthogonal vectors in the plane
    v1 = V[0, :]
    v2 = V[1, :]
    points_2d = np.zeros((points.shape[0], 2))
    points_2d[:, 0] = np.dot(centered, v1)
    points_2d[:, 1] = np.dot(centered, v2)

    # 2D Circle Fitting (Magic of Least Squares)
    x = points_2d[:, 0]
    y = points_2d[:, 1]
    A = np.c_[x, y, np.ones(len(x))]
    b = x**2 + y**2
    c = lstsq(A, b, cond=None)[0]
    
    # Get Circle Parameters
    xc = c[0] / 2
    yc = c[1] / 2
    r = np.sqrt(c[2] + xc**2 + yc**2)

    # Transform Center back to 3D
    center_3d = centroid + xc * v1 + yc * v2
    
    return center_3d, r, normal
    
def plot_3d_circle(ax, center, radius, normal, n_points=100):
    # Normalize the normal vector
    normal = np.array(normal) / norm(normal)

    # Generate an orthogonal vector 'u'
    # Use a non-parallel vector to cross with normal
    if np.allclose(normal, [0, 0, 1]) or np.allclose(normal, [0, 0, -1]):
        u = np.array([1, 0, 0])
    else:
        u = np.cross([0, 0, 1], normal)
    u = u / np.linalg.norm(u)

    # Generate the second orthogonal vector 'v'
    v = np.cross(normal, u)

    # Calculate circle points using parametric formula
    theta = np.linspace(0, 2 * np.pi, n_points)
    circle_points = (np.array(center)[:, np.newaxis] + 
                     radius * (np.outer(u, np.cos(theta)) + 
                               np.outer(v, np.sin(theta))))

    # Plot the result
    ax.plot(circle_points[0, :], circle_points[1, :], circle_points[2, :], color='red')

# Fit
fitted_center, fitted_r, fitted_normal = fit_circle_3d(v_vectors)
print(fitted_normal)

print("Hodograph Circle Parameters:", fitted_center[0], fitted_center[1], fitted_r)
L = norm( np.cross(x_vectors[-1], v_vectors[-1]) )
M_estimate = fitted_r * L / G
print("Estimate of Asteroid Mass", M_estimate )
print("Error in Estimate of Asteroid Mass", (M_estimate - np.sum(Asteroid_Mass))/np.sum(Asteroid_Mass))


fig = plt.figure()
ax = plt.axes(111, projection='3d')
ax.plot3D(x_vectors[:, 0], x_vectors[:, 1], x_vectors[:, 2], 'blue', label='position vectors')
ax.set_title("Flyby Trajectory")
ax.set_xlabel("Position (x component)")
ax.set_ylabel("Position (y component)")
ax.set_zlabel("Position (z component)")
for i in range( len(Asteroid_Mass) ):
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = Asteroid_Size[i] * np.cos(u)*np.sin(v) + Asteroid_Location[i, 0]
    y = Asteroid_Size[i] * np.sin(u)*np.sin(v) + Asteroid_Location[i, 1]
    z = Asteroid_Size[i] * np.cos(v) + Asteroid_Location[i, 2]
    ax.plot_wireframe(x, y, z, color="black")
r = max(np.ptp(x_vectors, axis = 0)) / 2
ax.set_xlim(np.mean(x_vectors[:, 0]) - r, np.mean(x_vectors[:, 0]) + r)
ax.set_ylim(np.mean(x_vectors[:, 1]) - r, np.mean(x_vectors[:, 1]) + r)
ax.set_zlim(np.mean(x_vectors[:, 2]) - r, np.mean(x_vectors[:, 2]) + r)
ax.set_box_aspect((1,1,1))
plt.show()

fig = plt.figure()
ax = plt.axes(111, projection='3d')
ax.plot3D(v_vectors[:, 0], v_vectors[:, 1], v_vectors[:, 2], color='blue', label='velocity vectors')
ax.set_title("Velocity Hodograph")
ax.set_xlabel("Velocity (x component)")
ax.set_ylabel("Velocity (y component)")
ax.set_zlabel("Velocity (z component)")
plot_3d_circle(ax, fitted_center, fitted_r, fitted_normal)
r = max(np.ptp(v_vectors, axis = 0)) / 2
ax.set_xlim(np.mean(v_vectors[:, 0]) - r, np.mean(v_vectors[:, 0]) + r)
ax.set_ylim(np.mean(v_vectors[:, 1]) - r, np.mean(v_vectors[:, 1]) + r)
ax.set_zlim(np.mean(v_vectors[:, 2]) - r, np.mean(v_vectors[:, 2]) + r)
ax.set_box_aspect((1,1,1))
plt.show()


KE = 0.5 * Probe_Mass * norm(v_vectors, axis = 1) * norm(v_vectors, axis = 1)
PE = []
for j in range( len(t_vals) ):
    PE_temp = 0
    for i in range( len(Asteroid_Mass) ):
        PE_temp -= G * Probe_Mass * Asteroid_Mass[i] / norm(x_vectors[j] - Asteroid_Location[i])
    PE.append(PE_temp)

plt.plot(t_vals, KE, label='Kinetic Energy')
plt.plot(t_vals, PE, label='Potential Energy')
plt.title("Energy Graph")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.grid(True)
plt.show()

plt.plot(t_vals, KE + PE, label='Total Energy')
plt.title("Total Mechanical Energy")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.grid(True)
plt.show()

print('Error in Energy:', np.ptp(KE + PE) / np.mean(KE + PE))
