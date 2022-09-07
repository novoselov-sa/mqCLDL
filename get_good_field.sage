def get_good_field(deg, lbound_start):

	continue_flag = true
	P = Primes()
	field = {}
	count = 0
	lbound_this = lbound_start
	for i in range(deg):
		while true:
			d = P.next(lbound_this)
			lbound_this = d + 5
			if ((-d) % 4 == 1):
				field[chr(ord('a')+count)] = -d
				count +=1
				break
	return field

print(get_good_field(3, 700))

