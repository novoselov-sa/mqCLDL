#SAGE_PATH=${HOME}/apps/SageMath-9.4
SAGE_PATH=$(shell dirname $(shell readlink -f $(shell which sage)))
$(info Sage Path: $(SAGE_PATH))

SAGE=$(SAGE_PATH)/sage

ifeq (,$(wildcard $(SAGE_PATH)/local/bin/python3))
PYTHON_PATH=python3
PYTHON_SYSTEM=true
else
PYTHON_PATH=$(SAGE_PATH)/local/bin/python3
PYTHON_SYSTEM=false
endif

#PYTHON_VERSION=3.9
$(info Python path: $(PYTHON_PATH))
PYTHON_VERSION=$(shell $(PYTHON_PATH) -c 'import platform; major, minor, patch = platform.python_version_tuple(); print(major + "." + minor)')
$(info Python version: $(PYTHON_VERSION))

PYTHON_SITE_PACKAGES_INCLUDE=$(shell $(PYTHON_PATH) -c "import site; print('-I' + ' -I'.join(site.getsitepackages()))")
#PYTHON_SITE_PACKAGES_PATH=/usr/local/lib/python$(PYTHON_VERSION)/site-packages
$(info Python site packages includes: $(PYTHON_SITE_PACKAGES_INCLUDE))

#CYTHON=$(SAGE) -cython -I $(SAGE_PATH)/local/lib64/python$(PYTHON_VERSION)/site-packages
ifeq ($(PYTHON_SYSTEM), true)
ifeq (, $(shell which cython3))
CYTHON=cython $(PYTHON_SITE_PACKAGES_INCLUDE)
else
CYTHON=cython3 $(PYTHON_SITE_PACKAGES_INCLUDE)
endif
else
CYTHON=$(SAGE) -cython $(PYTHON_SITE_PACKAGES_INCLUDE)
endif
$(info Cython: $(CYTHON))

ifeq ($(PYTHON_SYSTEM), true)
PYTHON_INCLUDE=/usr/include/python$(PYTHON_VERSION)
SAGE_INCLUDES=$(shell $(PYTHON_PATH) -c "import site; print('-I' + ' -I'.join([path + '/sage/cpython' for path in site.getsitepackages()]))") \
			  $(shell $(PYTHON_PATH) -c "import site; print('-I' + ' -I'.join([path + '/sage/libs/ntl' for path in site.getsitepackages()]))")
else
PYTHON_INCLUDE=$(SAGE_PATH)/local/include/python$(PYTHON_VERSION)
SAGE_SRC=$(SAGE_PATH)/src
SAGE_INCLUDES=-I$(SAGE_SRC) -I$(SAGE_SRC)/sage/libs/ntl -I$(SAGE_SRC)/sage/cpython
endif
$(info Python include: $(PYTHON_INCLUDE))
$(info Sage include: $(SAGE_INCLUDES))


CC=gcc -shared -pthread -fPIC -fwrapv -g -O2 -Wall -fno-strict-aliasing \
	   -I$(PYTHON_INCLUDE) \
	   $(PYTHON_SITE_PACKAGES_INCLUDE) \
	   $(SAGE_INCLUDES)

all: \
goodprime.so \
goodprimecheat.so \
hadamard.so \
centermod.so \
mult.so \
norm.so \
div.so \
subsetprod.so \
powerprod.so \
sqrt.so \
char.so \
ring.so \
field.so \
units.so \
ideal.so \
relations.so \
fb.so \
polynomial_ring.so \
prime_decomp.so \
trees.so \
cryptosystem.so \
verify.so \
saturation.so \
clgp.so \
ideal_sqrt_cyc.so \
dlog_quad.so \
field_compact.so \
sunits.so \
idealprod.so

goodprimecheat.so: goodprimecheat.c
	$(CC) -o goodprimecheat.so goodprimecheat.c

goodprimecheat.c: goodprimecheat.pyx
	$(CYTHON) goodprimecheat.pyx

goodprimecheat.pyx: goodprimecheat.sage.py
	cp goodprimecheat.sage.py goodprimecheat.pyx

goodprimecheat.sage.py: goodprimecheat.sage
	$(SAGE) -preparse goodprimecheat.sage

goodprime.so: goodprime.c
	$(CC) -o goodprime.so goodprime.c

goodprime.c: goodprime.pyx
	$(CYTHON) goodprime.pyx

goodprime.pyx: goodprime.sage.py
	cp goodprime.sage.py goodprime.pyx

goodprime.sage.py: goodprime.sage
	$(SAGE) -preparse goodprime.sage

hadamard.so: hadamard.c
	$(CC) -o hadamard.so hadamard.c

hadamard.c: hadamard.pyx
	$(CYTHON) hadamard.pyx

centermod.so: centermod.c
	$(CC) -o centermod.so centermod.c

centermod.c: centermod.pyx
	$(CYTHON) centermod.pyx

mult.so: mult.c
	$(CC) -o mult.so mult.c

mult.c: mult.pyx
	$(CYTHON) mult.pyx

mult.pyx: mult.sage.py
	cp mult.sage.py mult.pyx

mult.sage.py: mult.sage
	$(SAGE) -preparse mult.sage

norm.so: norm.c
	$(CC) -o norm.so norm.c

norm.c: norm.pyx
	$(CYTHON) norm.pyx

norm.pyx: norm.sage.py
	cp norm.sage.py norm.pyx

norm.sage.py: norm.sage
	$(SAGE) -preparse norm.sage

div.so: div.c
	$(CC) -o div.so div.c

div.c: div.pyx
	$(CYTHON) div.pyx

div.pyx: div.sage.py
	cp div.sage.py div.pyx

div.sage.py: div.sage
	$(SAGE) -preparse div.sage

subsetprod.so: subsetprod.c
	$(CC) -o subsetprod.so subsetprod.c

subsetprod.c: subsetprod.pyx
	$(CYTHON) subsetprod.pyx

subsetprod.pyx: subsetprod.sage.py
	cp subsetprod.sage.py subsetprod.pyx

subsetprod.sage.py: subsetprod.sage
	$(SAGE) -preparse subsetprod.sage

powerprod.so: powerprod.c
	$(CC) -o powerprod.so powerprod.c

powerprod.c: powerprod.pyx
	$(CYTHON) powerprod.pyx

powerprod.pyx: powerprod.sage.py
	cp powerprod.sage.py powerprod.pyx

powerprod.sage.py: powerprod.sage
	$(SAGE) -preparse powerprod.sage

sqrt.so: sqrt.c
	$(CC) -o sqrt.so sqrt.c

sqrt.c: sqrt.pyx
	$(CYTHON) sqrt.pyx

sqrt.pyx: sqrt.sage.py
	cp sqrt.sage.py sqrt.pyx

sqrt.sage.py: sqrt.sage
	$(SAGE) -preparse sqrt.sage

char.so: char.c
	$(CC) -o char.so char.c

char.c: char.pyx
	$(CYTHON) char.pyx

char.pyx: char.sage.py
	cp char.sage.py char.pyx

char.sage.py: char.sage
	$(SAGE) -preparse char.sage

ring.so: ring.c
	$(CC) -o ring.so ring.c

ring.c: ring.pyx
	$(CYTHON) ring.pyx

ring.pyx: ring.sage.py
	cp ring.sage.py ring.pyx

ring.sage.py: ring.sage
	$(SAGE) -preparse ring.sage

field.so: field.c
	$(CC) -o field.so field.c

field.c: field.pyx
	$(CYTHON) -3 field.pyx

field.pyx: field.sage.py
	cp field.sage.py field.pyx

field.sage.py: field.sage
	$(SAGE) -preparse field.sage

units.so: units.c
	$(CC) -o units.so units.c

units.c: units.pyx
	$(CYTHON) units.pyx

units.pyx: units.sage.py
	cp units.sage.py units.pyx

units.sage.py: units.sage
	$(SAGE) -preparse units.sage

ideal.so: ideal.c
	$(CC) -o ideal.so ideal.c

ideal.c: ideal.pyx
	$(CYTHON) ideal.pyx

ideal.pyx: ideal.sage.py
	cp ideal.sage.py ideal.pyx

ideal.sage.py: ideal.sage
	$(SAGE) -preparse ideal.sage

cryptosystem.so: cryptosystem.c
	$(CC) -o cryptosystem.so cryptosystem.c

cryptosystem.c: cryptosystem.pyx
	$(CYTHON) cryptosystem.pyx

cryptosystem.pyx: cryptosystem.sage.py
	cp cryptosystem.sage.py cryptosystem.pyx

cryptosystem.sage.py: cryptosystem.sage
	$(SAGE) -preparse cryptosystem.sage

relations.so: relations.c
	$(CC) -o relations.so relations.c

relations.c: relations.pyx
	$(CYTHON) relations.pyx

relations.pyx: relations.sage.py
	cp relations.sage.py relations.pyx

relations.sage.py: relations.sage
	$(SAGE) -preparse relations.sage

fb.so: fb.c
	$(CC) -o fb.so fb.c

fb.c: fb.pyx
	$(CYTHON) fb.pyx

fb.pyx: fb.sage.py
	cp fb.sage.py fb.pyx

fb.sage.py: fb.sage
	$(SAGE) -preparse fb.sage

polynomial_ring.so: polynomial_ring.c
	$(CC) -o polynomial_ring.so polynomial_ring.c

polynomial_ring.c: polynomial_ring.pyx
	$(CYTHON) polynomial_ring.pyx

polynomial_ring.pyx: polynomial_ring.sage.py
	cp polynomial_ring.sage.py polynomial_ring.pyx

polynomial_ring.sage.py: polynomial_ring.sage
	$(SAGE) -preparse polynomial_ring.sage

prime_decomp.so: prime_decomp.c
	$(CC) -o prime_decomp.so prime_decomp.c

prime_decomp.c: prime_decomp.pyx
	$(CYTHON) prime_decomp.pyx

prime_decomp.pyx: prime_decomp.sage.py
	cp prime_decomp.sage.py prime_decomp.pyx

prime_decomp.sage.py: prime_decomp.sage
	$(SAGE) -preparse prime_decomp.sage

trees.so: trees.c
	$(CC) -o trees.so trees.c

trees.c: trees.pyx
	$(CYTHON) -3 trees.pyx

trees.pyx: trees.sage.py
	cp trees.sage.py trees.pyx

trees.sage.py: trees.sage
	$(SAGE) -preparse trees.sage

clgp.so: clgp.c
	$(CC) -o clgp.so clgp.c

clgp.c: clgp.pyx
	$(CYTHON) -3 clgp.pyx

clgp.pyx: clgp.sage.py
	cp clgp.sage.py clgp.pyx

clgp.sage.py: clgp.sage
	$(SAGE) -preparse clgp.sage

verify.so: verify.c
	$(CC) -o verify.so verify.c

verify.c: verify.pyx
	$(CYTHON) -3 verify.pyx

verify.pyx: verify.sage.py
	cp verify.sage.py verify.pyx

verify.sage.py: verify.sage
	$(SAGE) -preparse verify.sage

saturation.so: saturation.c
	$(CC) -o saturation.so saturation.c

saturation.c: saturation.pyx
	$(CYTHON) -3 saturation.pyx

saturation.pyx: saturation.sage.py
	cp saturation.sage.py saturation.pyx

saturation.sage.py: saturation.sage
	$(SAGE) -preparse saturation.sage

ideal_sqrt_cyc.so: ideal_sqrt_cyc.c
	$(CC) -o ideal_sqrt_cyc.so ideal_sqrt_cyc.c

ideal_sqrt_cyc.c: ideal_sqrt_cyc.pyx
	$(CYTHON) -3 ideal_sqrt_cyc.pyx

ideal_sqrt_cyc.pyx: ideal_sqrt_cyc.sage.py
	cp ideal_sqrt_cyc.sage.py ideal_sqrt_cyc.pyx

ideal_sqrt_cyc.sage.py: ideal_sqrt_cyc.sage
	$(SAGE) -preparse ideal_sqrt_cyc.sage

dlog_quad.so: dlog_quad.c
	$(CC) -o dlog_quad.so dlog_quad.c

dlog_quad.c: dlog_quad.pyx
	$(CYTHON) -3 dlog_quad.pyx

dlog_quad.pyx: dlog_quad.sage.py
	cp dlog_quad.sage.py dlog_quad.pyx

dlog_quad.sage.py: dlog_quad.sage
	$(SAGE) -preparse dlog_quad.sage

field_compact.so: field_compact.c
	$(CC) -o field_compact.so field_compact.c

field_compact.c: field_compact.pyx
	$(CYTHON) -3 field_compact.pyx

field_compact.pyx: field_compact.sage.py
	cp field_compact.sage.py field_compact.pyx

field_compact.sage.py: field_compact.sage
	$(SAGE) -preparse field_compact.sage

sunits.so: sunits.c
	$(CC) -o sunits.so sunits.c

sunits.c: sunits.pyx
	$(CYTHON) -3 sunits.pyx

sunits.pyx: sunits.sage.py
	cp sunits.sage.py sunits.pyx

sunits.sage.py: sunits.sage
	$(SAGE) -preparse sunits.sage

idealprod.so: idealprod.c
	$(CC) -o idealprod.so idealprod.c

idealprod.c: idealprod.pyx
	$(CYTHON) -3 idealprod.pyx

idealprod.pyx: idealprod.sage.py
	cp idealprod.sage.py idealprod.pyx

idealprod.sage.py: idealprod.sage
	$(SAGE) -preparse idealprod.sage

clean:
	rm -f goodprimecheat.so goodprimecheat.c goodprimecheat.pyx
	rm -f goodprime.so goodprime.c goodprime.pyx
	rm -f hadamard.so hadamard.c
	rm -f centermod.so centermod.c
	rm -f mult.so mult.c mult.pyx
	rm -f norm.so norm.c norm.pyx
	rm -f div.so div.c div.pyx
	rm -f subsetprod.so subsetprod.c subsetprod.pyx
	rm -f powerprod.so powerprod.c powerprod.pyx
	rm -f sqrt.so sqrt.c sqrt.pyx
	rm -f char.so char.c char.pyx
	rm -f ring.so ring.c ring.pyx
	rm -f field.so field.c field.pyx
	rm -f units.so units.c units.pyx
	rm -f relations.so relations.c relations.pyx
	rm -f fb.so fb.c fb.pyx
	rm -f polynomial_ring.so polynomial_ring.c polynomial_ring.pyx
	rm -f prime_decomp.so prime_decomp.c prime_decomp.pyx
	rm -f trees.so trees.c trees.pyx
	rm -f ideal.so ideal.c ideal.pyx
	rm -f clgp.so clgp.c clgp.pyx
	rm -f cryptosystem.so cryptosystem.c cryptosystem.pyx
	rm -f ideal_sqrt_cyc.so ideal_sqrt_cyc.c ideal_sqrt_cyc.pyx
	rm -f dlog_quad.so dlog_quad.c dlog_quad.pyx
	rm -f field_compact.so field_compact.c field_compact.pyx
	rm -f verify.so verify.c verify.pyx
	rm -f saturation.so saturation.c saturation.pyx
	rm -f sunits.so sunits.c sunits.pyx
	rm -f idealprod.so idealprod.c idealprod.pyx
	for f in in tests/*test_sage.py; do\
		rm -f "$$f";\
	done
	rm -f *.sage.py

# Running unittests. Unittest discover doesn't support "." in patterns, so we have to rename files first.
tests: all
	$(SAGE) -preparse tests/*_test.sage
	for f in tests/*.sage.py; do\
        mv "$$f" "$$(echo $${f} | sed s/.sage.py/_sage.py/)";\
    done
	$(SAGE) -python -m unittest discover -s tests/ -p "*test_sage.py"
