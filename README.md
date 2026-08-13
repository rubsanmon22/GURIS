# TP-**GURIS** (Test Particle Guiding-center Unified Relativistic Integrator for MHD Simulations) 

## **GURIS** is a Fortran test-particle integrator. It interpolates 3D/2.5D MHD fields via Cloud-in-Cell and computes relativistic particle trajectories using Runge-Kutta methods, designed for solar flare and magnetic reconnection research.




## Particulas de prueba con aproximacion de centro guia y lectura TXT

Esta version evoluciona particulas de prueba en la aproximacion de centro guia. No reproduce el movimiento de giro alrededor de B.

Las variables de cada particula son:

```text
X, Y, Z, v_parallel, mu
```

La energia cinetica se calcula como:

```text
Ekin = 1/2 m v_parallel^2 + mu B
```

## Formato del TXT de campos

El lector espera una fila por punto de grilla con columnas:

```text
x y z Bx By Bz Ex Ey Ez
```

El orden de filas esperado es:

```python
for z in zgrid:
    for y in ygrid:
        for x in xgrid:
            write(x, y, z, Bx, By, Bz, Ex, Ey, Ez)
```

Es decir, el indice `i` de x corre mas rapido, luego y, luego z.

## Test incluido

El script:

```bash
python3 scripts/make_uniform_txt.py
```

genera `fields_uniform.txt` con:

```text
B = (0, 0, 1)
E = (0.1, 0, 0)
```

En unidades normalizadas, la deriva de ExB esperada es:

```text
v_E = E x B / B^2 = (0, -0.1, 0)
```

## Compilar

Con OpenMP:

```bash
make
```

Sin OpenMP:

```bash
make OMP=0
```

## Correr el ejemplo con TXT

```bash
make fields
make
./gc_particles
```

Tambien puede hacerse todo junto con:

```bash
make run
```

## Graficar energias

```bash
python3 scripts/plot_energy.py
```

Esto genera:

```text
output/gc_txt_energy_distribution.png
output/gc_txt_energy_histogram_from_fortran.png
```

## Salidas

```text
output/gc_txt_particles_final.dat
output/gc_txt_energy_table.dat
output/gc_txt_energy_histogram.dat
```

La tabla de energias tiene columnas:

```text
# id E_initial E_final active
```

El histograma escrito por Fortran tiene columnas:

```text
# E_center N_initial N_final
```

## Cambiar entre campo interno y TXT

En `input.nml`, para leer TXT:

```fortran
field_source = "file"
file_type    = "txt"
field_file   = "fields_uniform.txt"
```

Para usar el campo uniforme interno:

```fortran
field_source = "analytic"
```

## Como pasar a una simulacion real

El integrador y el CIC no deben cambiar. Lo que cambia es el modulo lector.

La regla de diseno es que cualquier lector debe llenar la misma estructura `grid_t`:

```fortran
grid%x, grid%y, grid%z
grid%bx, grid%by, grid%bz
grid%ex, grid%ey, grid%ez
```

Despues se ejecuta siempre:

```fortran
call compute_derived_fields(grid)
call initialize_particles(cfg, grid, part)
call advance_guiding_center_rk2(...)
```

Si la simulacion MHD no guarda el campo electrico, en MHD ideal se puede reconstruir como:

```text
E = - u x B
```

cuidando las unidades del archivo. En cgs puede aparecer un factor `1/c`; en unidades normalizadas muchas veces no.

## Limitaciones de esta version

* El lector TXT asume grilla cartesiana uniforme.
* El CIC usa una celda regular con pesos de primer orden.
* Para grillas no uniformes habria que reemplazar la busqueda directa de indices por una busqueda binaria.
* Para AMR conviene primero interpolar la salida de la simulacion a una grilla cartesiana uniforme.
