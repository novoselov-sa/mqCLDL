
food = {'a':-11, 'b':-19, 'c':-23, 'd':-31, 'e':-43}
d = list(food.values())
rho = var('rho', n=2^len(d) + 1)
