"""
Visualize 31 R7081 PMT positions inside the prototype buffer tank.
PMT volumes are drawn as cylinders (envelope solid).
"""

from matplotlib.patches import Patch
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ---------------------------------------------------------------------------
# Geometry constants (mm)
# ---------------------------------------------------------------------------
TANK_R_OUT = 1100.0
TANK_HH_OUT = 1093.0
WALL_T = 6.0
TANK_R_IN = TANK_R_OUT - WALL_T   # 1094.0
TANK_HH_IN = TANK_HH_OUT - WALL_T   # 1402.0

PMT_R = 126.5   # envelope radius
PMT_HH = 165.0   # envelope half-height
OFFSET = 0.0   # gap between PMT back and inner wall (mm)

Z_TOP = TANK_HH_IN - PMT_HH - OFFSET   # +1187.0
Z_BOT = -(TANK_HH_IN - PMT_HH - OFFSET)  # -1187.0
R_BARREL = TANK_R_IN - PMT_HH - OFFSET  # 879.0

R_RING = 500.5   # top/bottom ring radius

# ---------------------------------------------------------------------------
# PMT list: (cx, cy, cz, face_x, face_y, face_z)
# face = unit vector pointing from photocathode into the tank
# ---------------------------------------------------------------------------
# pmts: (cx, cy, cz, face_x, face_y, face_z, nring, region)
# region: 1=upper, 0=barrel, -1=bottom
pmts = []

# Top plate: 6 PMTs, face -z, region=1, nring=1
for i in range(6):
    phi = np.radians(i * 60)
    pmts.append((R_RING * np.cos(phi), R_RING *
                np.sin(phi), Z_TOP, 0, 0, -1, 1, 1))

# Barrel: 3 rows × 6 columns (z = +600, 0, -600), face inward, region=0, nring=2,3,4
for nring, dz in enumerate([600, 0, -600], start=2):
    for i in range(6):
        phi = np.radians(i * 60)
        cx = R_BARREL * np.cos(phi)
        cy = R_BARREL * np.sin(phi)
        pmts.append((cx, cy, dz, -np.cos(phi), -np.sin(phi), 0, nring, 0))

# Bottom plate: 7 PMTs, face +z, region=-1
for i in range(6):                                     # ring, nring=5
    phi = np.radians(i * 60)
    pmts.append((R_RING * np.cos(phi), R_RING *
                np.sin(phi), Z_BOT, 0, 0, 1, 5, -1))
pmts.append((0, 0, Z_BOT, 0, 0, 1, 6, -1))            # center, nring=6

assert len(pmts) == 31

# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------


def cylinder_surfaces(cx, cy, cz, r, hh, face):
    """Return (side_X, side_Y, side_Z, cap1_verts, cap2_verts) for a cylinder."""
    face = np.asarray(face, dtype=float)
    face /= np.linalg.norm(face)

    # Orthonormal basis
    ref = np.array([0, 0, 1]) if abs(face[2]) < 0.9 else np.array([1, 0, 0])
    u = np.cross(face, ref)
    u /= np.linalg.norm(u)
    v = np.cross(face, u)

    theta = np.linspace(0, 2 * np.pi, 40)
    center = np.array([cx, cy, cz])
    c1 = center + hh * face
    c2 = center - hh * face

    ring1 = c1[:, None] + r * \
        (np.cos(theta) * u[:, None] + np.sin(theta) * v[:, None])
    ring2 = c2[:, None] + r * \
        (np.cos(theta) * u[:, None] + np.sin(theta) * v[:, None])

    # Side surface mesh (2 × N)
    side = np.stack([ring1, ring2], axis=1)  # (3, 2, N)
    sX, sY, sZ = side[0], side[1], side[2]

    # Cap polygon vertices
    def cap_poly(center, ring):
        verts = np.column_stack([ring[0], ring[1], ring[2]])  # (N, 3)
        verts = np.vstack([verts, verts[0]])                  # close loop
        return [verts]

    cap1 = cap_poly(c1, ring1)
    cap2 = cap_poly(c2, ring2)

    return sX, sY, sZ, cap1, cap2


def draw_pmt(ax, cx, cy, cz, face, color, alpha=0.75):
    sX, sY, sZ, cap1, cap2 = cylinder_surfaces(cx, cy, cz, PMT_R, PMT_HH, face)
    ax.plot_surface(sX, sY, sZ, color=color, alpha=alpha, linewidth=0)
    for cap in (cap1, cap2):
        poly = Poly3DCollection(cap, color=color, alpha=alpha)
        ax.add_collection3d(poly)


def draw_tank(ax, r, hh, color='steelblue', alpha=0.08, n=80):
    theta = np.linspace(0, 2 * np.pi, n)
    z = np.array([-hh, hh])
    T, Z = np.meshgrid(theta, z)
    X = r * np.cos(T)
    Y = r * np.sin(T)
    ax.plot_surface(X, Y, Z, color=color, alpha=alpha, linewidth=0)

    # top and bottom discs
    rr = np.linspace(0, r, 10)
    TT, RR = np.meshgrid(theta, rr)
    Xd = RR * np.cos(TT)
    Yd = RR * np.sin(TT)
    for zd in [-hh, hh]:
        ax.plot_surface(Xd, Yd, np.full_like(Xd, zd),
                        color=color, alpha=alpha * 1.5, linewidth=0)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
COLORS = {'top': '#2196F3', 'bottom': '#F44336', 'barrel': '#4CAF50'}

fig = plt.figure(figsize=(10, 13))
ax = fig.add_subplot(111, projection='3d')

draw_tank(ax, TANK_R_OUT, TANK_HH_OUT)

for i, (cx, cy, cz, fx, fy, fz, nring, region) in enumerate(pmts):
    if region == 1:
        color = COLORS['top']
    elif region == 0:
        color = COLORS['barrel']
    else:
        color = COLORS['bottom']
    draw_pmt(ax, cx, cy, cz, (fx, fy, fz), color)

# Axis labels and limits
ax.set_xlabel('X (mm)', labelpad=8)
ax.set_ylabel('Y (mm)', labelpad=8)
ax.set_zlabel('Z (mm)', labelpad=8)
ax.set_title('Prototype Buffer Tank — 31 R7081 PMTs', pad=12)

lim = TANK_R_OUT * 1.15
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_zlim(-TANK_HH_OUT * 1.15, TANK_HH_OUT * 1.15)
ax.set_box_aspect([1, 1, TANK_HH_OUT / TANK_R_OUT])

# Legend
ax.legend(handles=[
    Patch(color=COLORS['top'],    label='Top plate  (6 PMTs)'),
    Patch(color=COLORS['bottom'], label='Bottom plate (7 PMTs)'),
    Patch(color=COLORS['barrel'], label='Barrel (18 PMTs)'),
], loc='upper left', fontsize=10)

plt.tight_layout()
#plt.savefig('pmt_layout.png', dpi=150, bbox_inches='tight')
#print('Saved: pmt_layout.png')

# ---------------------------------------------------------------------------
# Print PMT positions
# ---------------------------------------------------------------------------
print(f"\n{'#pmtno':>6}  {'x':>10}  {'y':>10}  {'z':>10}  {'nring':>6}  {'region':>7}")
print('#-' * 30)
for i, (cx, cy, cz, fx, fy, fz, nring, region) in enumerate(pmts):
    print(f"{i+1:>6}  {cx:>10.2f}  {cy:>10.2f}  {cz:>10.2f}  {nring:>6}  {region:>7}")

plt.show()
