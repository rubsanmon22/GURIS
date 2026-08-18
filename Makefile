FC = gfortran
OMP ?= 1
SRC_DIR = src
BUILD_DIR = build
EXE = gc_particles

FFLAGS_BASE = -O3 -std=f2008 -Wall -Wextra -ffree-line-length-none -J$(BUILD_DIR)
ifeq ($(OMP),1)
  FFLAGS = $(FFLAGS_BASE) -fopenmp
  LDFLAGS = -fopenmp
else
  FFLAGS = $(FFLAGS_BASE)
  LDFLAGS =
endif

SOURCES = \
  $(SRC_DIR)/mod_precision.f90 \
  $(SRC_DIR)/mod_config.f90 \
  $(SRC_DIR)/mod_grid.f90 \
  $(SRC_DIR)/mod_fields_analytic.f90 \
  $(SRC_DIR)/mod_reader_txt.f90 \
  $(SRC_DIR)/mod_reader.f90 \
  $(SRC_DIR)/mod_derived_fields.f90 \
  $(SRC_DIR)/mod_random.f90 \
  $(SRC_DIR)/mod_cic.f90 \
  $(SRC_DIR)/mod_particles.f90 \
  $(SRC_DIR)/mod_boundary.f90 \
  $(SRC_DIR)/mod_evolution_guiding_center.f90 \
  $(SRC_DIR)/mod_diagnostics.f90 \
  $(SRC_DIR)/mod_output.f90 \
  $(SRC_DIR)/main.f90

OBJECTS = $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SOURCES))

all: $(EXE)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(EXE): $(OBJECTS)
	$(FC) $(OBJECTS) $(LDFLAGS) -o $(EXE)

fields:
	python3 scripts/make_uniform_txt.py

dipole:
	python3 scripts/make_dipole_txt.py

run: all fields dipole
	./$(EXE)

plot:
	python3 scripts/plot_energy_loglog.py

clean:
	rm -rf $(BUILD_DIR) $(EXE)

distclean: clean
	rm -rf output fields_uniform.txt fields_dipole.txt
