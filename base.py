
digits: str = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+/"

def to_positional_base(i: int, b: int) -> str:
	sign: str
	ret: str
	if 1 < b <= len(digits):
		if i > 0:
			sign = ""
		elif i == 0:
			return "0"
		else:
			i = -i
			sign = "-"
		ret = ""
		while i:
			ret = digits[i%b] + ret
			i //= b
		return sign + ret
	elif b == 1:
		if i >= 0:
			return "1"*i
		else:
			raise ValueError("Cannot convert negative integer to base 1")
	elif b == -1:
		if i in [0, -1]:
			return "1"*-i
		else:
			raise ValueError("Cannot convert integer != 0 or 1 to base -1")
	elif -1 > b >= -len(digits):
		if i == 0:
			return "0"
		else:
			ret = ""
			while i:
				d = i % b
				i //= b
				if d < 0:
					d -= b
					i += 1
				ret = digits[abs(d)] + ret
			return ret
	else:
		raise ValueError("Base %d not supported" % b)

def from_positional_base(s: str, b: int) -> int:
	if b and -len(digits) <= b <= len(digits):
		if s.startswith("-"):
			return -sum(digits.index(digit) * b**e for e, digit in enumerate(s[1::-1]))
		else:
			return sum(digits.index(digit) * b**e for e, digit in enumerate(s[::-1]))
	else:
		raise ValueError("Base %d not supported" % b)
