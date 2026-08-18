import numpy as np

# ==============================================================================
# 1. PARÁMETROS DE CONFIGURACIÓN, DIMENSIONES Y RESOLUCIÓN
# ==============================================================================
# Nombre del archivo de salida
OUTPUT_FILE = "fields_dipole.txt"

# Dimensiones del dominio espacial (e.g. en cm, m o unidades de radio solar)
xmin, xmax = -10.0, 10.0
ymin, ymax = -10.0, 10.0
zmin, zmax = -10.0, 10.0

# RESOLUCIÓN: Número de puntos en cada eje (puedes modificar estos valores)
Nx = 32
Ny = 32
Nz = 32

# Parámetros del dipolo magnético
# Momento dipolar m = (mx, my, mz)
mx, my, mz = 0.0, 0.0, 1000.0  # Dipolo orientado en Z
mu0_over_4pi = 1e-7           # Factor de escala (1.0 para CGS/Gauss, 1e-7 para SI)
epsilon = 1e-3                 # Suavizado para evitar división por cero en r = 0

# ==============================================================================
# 2. GENERACIÓN DE LA MALLA TRIDIMENSIONAL
# ==============================================================================
x = np.linspace(xmin, xmax, Nx)
y = np.linspace(ymin, ymax, Ny)
z = np.linspace(zmin, zmax, Nz)

# Crear la grilla tridimensional
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Distancia radial desde el origen r = sqrt(x^2 + y^2 + z^2 + eps^2)
R_sq = X**2 + Y**2 + Z**2 + epsilon**2
R = np.sqrt(R_sq)

# ==============================================================================
# 3. CÁLCULO DEL CAMPO MAGNÉTICO (B) Y ELÉCTRICO (E)
# ==============================================================================
# Producto punto (m . r)
m_dot_r = mx * X + my * Y + mz * Z

# Ecuación dipolar vectorizada: B = (mu0/4pi) * [ 3*(m . r)*r / r^5 - m / r^3 ]
Bx = mu0_over_4pi * (3.0 * m_dot_r * X / (R**5) - mx / (R**3))
By = mu0_over_4pi * (3.0 * m_dot_r * Y / (R**5) - my / (R**3))
Bz = mu0_over_4pi * (3.0 * m_dot_r * Z / (R**5) - mz / (R**3))

# Campo eléctrico (por defecto 0.0 si es un campo magnético estático)
Ex = np.zeros_like(Bx)
Ey = np.zeros_like(By)
Ez = np.zeros_like(Bz)

# ==============================================================================
# 4. EXPORTACIÓN A ARCHIVO TXT
# ==============================================================================
# Aplanar las matrices 3D a vectores de 1 columna
data = np.column_stack([
    X.ravel(), Y.ravel(), Z.ravel(),
    Bx.ravel(), By.ravel(), Bz.ravel(),
    Ex.ravel(), Ey.ravel(), Ez.ravel()
])

# Encabezado del archivo
header = "x y z Bx By Bz Ex Ey Ez"

# Formato numérico (notación científica con 6 decimales)
fmt = "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e"

# Guardar en disco
np.savetxt(OUTPUT_FILE, data, header=header, comments='', fmt=fmt)

print(f"-> Archivo '{OUTPUT_FILE}' generado con éxito.")
print(f"-> Puntos totales: {Nx * Ny * Nz} ({Nx}x{Ny}x{Nz})")
print(f"-> Rango espacial: X[{xmin}, {xmax}], Y[{ymin}, {ymax}], Z[{zmin}, {zmax}]")
