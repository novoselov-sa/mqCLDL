This is a mqCLDL code build on top of:
* https://scarecryptow.org/publications/multiclass.html

# How to compile
1. Run ```cp "Makefile(template)" Makefile```
2. Run ```make```
3. Edit paths in Makefile according to your installation if code above fails.

# Folder structure
* tests/ - unit tests
* trees/ - trees describing prime ideal splitting over subfields. Factor base is loaded from this files.
* relations/ - class group relation matrices for every subfield
