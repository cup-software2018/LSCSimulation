import numpy as np
import matplotlib.pyplot as plt

# R7081 Data from C++ code
z_edge = np.array([96.7, 40.0, 0.0, -40.0, -90.0, -142.0, -233.3])
rho_edge = np.array([0.0, 111.0, 126.5, 111.0, 42.25, 42.25, 42.25])
z_o = np.array([-40.0, 0.0, 0.0, 40.0, -142.0, -233.3])
n = len(z_o)

a = np.zeros(n)
b = np.zeros(n)

# 1. Translate SetAllParameters: Calculate swept radius (a) and cross-section radius (b)
for i in range(n):
    if rho_edge[i] == rho_edge[i+1]:  # Cylinder
        a[i] = rho_edge[i]
        b[i] = 0.0
    else:  # Torus section
        dz = z_edge[i+1] - z_edge[i]
        dr = rho_edge[i+1] - rho_edge[i]
        a[i] = (dz * (z_edge[i+1] + z_edge[i] - 2*z_o[i]) /
                dr + (rho_edge[i+1] + rho_edge[i])) / 2.0
        b[i] = np.sqrt((z_edge[i] - z_o[i])**2 + (rho_edge[i] - a[i])**2)
        # Check concavity
        if rho_edge[i] < a[i] and rho_edge[i+1] < a[i]:
            b[i] = -b[i]

# 2. Translate MakeSegment: Generate exact surface points
rr_list = []
zz_list = []
ns = 50  # Number of steps per segment

for i in range(n):
    rho0, rho1 = rho_edge[i], rho_edge[i+1]
    z0, z1 = z_edge[i], z_edge[i+1]
    zo = z_o[i]
    aa = a[i]
    bb = b[i]

    if rho0 == rho1:
        zz_seg = np.linspace(z0, z1, ns)
        rr_seg = np.full(ns, rho0)
    else:
        dz = abs(z0 - z1)
        dr = abs(rho0 - rho1)
        if dz > dr:
            zz_seg = np.linspace(z0, z1, ns)
            dtmp = 1.0 - ((zz_seg - zo) / bb)**2
            dtmp[dtmp < 0] = 0
            rr_seg = aa + bb * np.sqrt(dtmp)
        else:
            rr_seg = np.linspace(rho0, rho1, ns)
            dtmp = 1.0 - ((rr_seg - aa) / bb)**2
            dtmp[dtmp < 0] = 0
            sign = 1 if z0 > zo else -1
            zz_seg = zo + abs(bb) * np.sqrt(dtmp) * sign

    rr_list.extend(rr_seg)
    zz_list.extend(zz_seg)

rr = np.array(rr_list)
zz = np.array(zz_list)

# 3. Plotting
fig, ax = plt.subplots(figsize=(6, 10))

ax.plot(rr, zz, 'b-', linewidth=2, label='Exact GLG4TorusStack Profile')
ax.plot(-rr, zz, 'b-', linewidth=2)

ax.plot(rho_edge, z_edge, 'ro', markersize=5, label='Edge Points')
ax.plot(-rho_edge, z_edge, 'ro', markersize=5)

ax.fill_betweenx(zz, -rr, rr, color='lightblue', alpha=0.5)

ax.set_aspect('equal')
ax.set_xlabel('Radius (Rho) [mm]')
ax.set_ylabel('Height (Z) [mm]')
ax.set_title('Hamamatsu R7081 - Exact GEANT4 Geometry')
ax.grid(True, linestyle='--', alpha=0.7)
ax.axhline(0, color='red', linestyle='--', linewidth=1, label='Equator (Z=0)')
ax.axvline(0, color='black', linewidth=1)
#ax.legend()

plt.show()
