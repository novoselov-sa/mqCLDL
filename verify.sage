# Levels of verification checks

# Disable all checks
NONE = -1

# Light checks and checks based on class group computation for subfields with n = len(d) <= 5
LIGHT = 0

# Checks based on class group computation for all subfields
MEDIUM = 1

# Checks based on explicit computation of ideal products
HEAVY = 2

_level = NONE

def level():
    return _level

def set(level):
    global _level
    assert level in [NONE, LIGHT, MEDIUM, HEAVY]
    _level = level
