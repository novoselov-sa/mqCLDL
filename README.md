This is a source code for the paper "On the discrete logarithm problem in the ideal class group of multiquadratic fields".

It is built on top of the implementation of the Biasse-van Vredendaal algorithm for the class group computation for real multiquadratic fields:
* https://scarecryptow.org/publications/multiclass.html
and on it's adaptation to the imaginary multiquadratic fields from:
* https://github.com/novoselov-sa/multiclass-im

# Requirements
* Sage 9.5+

# How to compile
1. Run ```make```
2. Edit paths in Makefile according to your installation if code above fails.

# Folder structure
* tests/ - unit tests
* trees/ - trees describing prime ideal splitting over subfields. Factor base is loaded from these files. They are generated by ```trees_generation.sage``` script.
* relations/ - class group relations matrices for every subfield. They are generated by ```testrelations.sage``` script.

# How to use
1. Set the target field in the file ```trees_generation.sage``` (```food``` variable). Alternatively, you can pass d_1 d_2 ... d_n as command line arguments of this script.
2. Run ```sage trees_generation.sage``` to generate trees.
3. Run ```sage testrelations.sage``` to compute class group and relation matrices for subfields.
4. Run ```sage dlog_cyc_PoC.sage``` to compute discrete logarithm for a random ideal using algorithms from the paper.
5. Run ```sage dlog_smooth.sage``` to compute discrete logarithm assuming that target ideal factors over factor base. You can add primes dividing norm of the target ideal to ```PRIMES``` variable in ```trees_generation.sage``` to include them to the factor base. This requires recomputation of trees and relation matrices.

# Versions
To download all code from this repo use link:
* [mqCLDL-v2.0.1.zip](releases/mqCLDL-v2.0.1.zip)
