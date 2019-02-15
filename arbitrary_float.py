from bits import bits
import struct
import math
import base as base_module

color = True

def continued_fraction(a, b, base=2): # https://possiblywrong.wordpress.com/2017/09/30/digits-of-pi-and-python-generators/
	"""Generate digits of continued fraction a(0)+b(1)/(a(1)+b(2)/(...)."""
	(p0, q0), (p1, q1) = (a(0), 1), (a(1) * a(0) + b(1), a(1))
	k = 1
	while True:
		(d0, r0), (d1, r1) = divmod(p0, q0), divmod(p1, q1)
		if d0 == d1:
			yield d1
			p0, p1 = base * r0, base * r1
		else:
			k = k + 1
			x, y = a(k), b(k)
			(p0, q0), (p1, q1) = (p1, q1), (x * p1 + y * p0, x * q1 + y * q0)

class TODO(Exception):
	pass

class ArbitraryFloatType(type):
	types = dict() # mapping from tuple(exp_len, mant_len) to class
	base = None
	def __new__(cls, *args):
		if cls.base is None: # ArbitraryFloatBase creation
			cls.base = type.__new__(cls, *args)
			return cls.base
		elif len(args) == 2 and isinstance(args[0], int) and args[0] > 0 and isinstance(args[1], int) and args[1] > 0:
			if args in cls.types:
				return cls.types[args]
			new = type.__new__(cls, "ArbitraryFloat(%d,%d)" % args, (cls.base,), {"exp_len": args[0], "mant_len": args[1]})
			cls.types[args] = new
			return new
		else: # things like long double
			try:
				if len(args) == 3:
					return type.__new__(cls, *args)
			except TypeError:
				raise
			raise ValueError(args)
	def __init__(cls, *args):
		pass
	@staticmethod
	def least_precision(*args):
		"""
		If given one argument:
			Return a type that can represent the argument exactly
			Always ArbitraryFloatType(11, 52) for floats
			For ints, returns a type that can represent every integer in [0,val]
		If given two arguments:
			Return a type that can represent the arguments exactly as a normal value
			arg[0] is an int exp
			arg[2] is a normalized bits mantissa
		"""
		if len(args) == 0:
			return ArbitraryFloatType(1, 1)
		elif len(args) == 1:
			val = args[0]
			if isinstance(val, float):
				return ArbitraryFloatType(11, 52) # IEEE 754 Double
			elif isinstance(val, int):
				if abs(val) <= 1:
					return ArbitraryFloatType(1, 1)
				return ArbitraryFloatType(1+(val.bit_length()-1).bit_length(), val.bit_length()-1)
			elif isinstance(val, ArbitraryFloatBase):
				return type(val)
			else:
				raise TypeError(val)
		elif len(args) == 2:
			exp, mant = args
			if not isinstance(exp, int) or not isinstance(mant, bits):
				raise TypeError("least_precision(exp, mant) takes (int, bits) not %r" % (tuple(type(arg) for arg in args),))
			if len(mant) == 0 or not mant[0]:
				raise ValueError("Mantissa must be normalized for least_precision(exp, mant), not %r" % mant)
			if exp < 0:
				exp = -exp + 1 # max_normal_exp = -min_normal_exp + 1
			exp_len = exp.bit_length() + 1 if exp else 2 # exp_len of 1 has only subnormal values
			return ArbitraryFloatType(exp_len, len(mant)-1)
			
	@staticmethod
	def best_precision(a, b):
		"Return a type that can represent both a and b exactly"
		if not isinstance(a, (ArbitraryFloatBase, ArbitraryFloatType)):
			a = ArbitraryFloatType.least_precision(a)
		if not isinstance(b, (ArbitraryFloatBase, ArbitraryFloatType)):
			b = ArbitraryFloatType.least_precision(b)
		return ArbitraryFloatType(max(a.exp_len, b.exp_len), max(a.mant_len, b.mant_len))
	@property
	def bias(self):
		"The (positive) exponent bias\nThe real exponent is (exp - bias)"
		return 2**(self.exp_len-1)-1
	@property
	def min_normal_exp(self):
		return 1 - self.bias
	@property
	def max_normal_exp(self):
		# This used to be wrong, so I need to check everywhere that uses it
		#return 2**self.exp_len - 1 - self.bias
		return 2**self.exp_len - 2 - self.bias
	@property
	def min_subnormal_exp(self):
		return 1 - self.bias - self.mant_len
	@property
	def max_subnormal_exp(self):
		return -self.bias
	
	@property
	def inf(self):
		return self(~bits(self.exp_len), bits(self.mant_len))
	@property
	def nan(self):
		return self(~bits(self.exp_len), ~bits(self.mant_len))
	@property
	def zero(self):
		return self(bits(self.exp_len), bits(self.mant_len))
	@property
	def maxint(self):
		if self.max_normal_exp <= self.mant_len:
			return 2**(self.max_normal_exp+1)-1
		else:
			return 2**(self.mant_len+1)
	@property
	def pi(self):
		pi_cont = continued_fraction(
			lambda k: 0 if k == 0 else 2 * k - 1,
			lambda k: 4 if k == 1 else (k - 1)**2,
			2
		)
		next(pi_cont) # 3, all following are [0,1]
		mant = bits([1,1] + [next(pi_cont) for i in range(self.mant_len-1)]) # 1+mant_len is the necesary length == 2 + mant_len-1
		if next(pi_cont): # just after the end
			mant = mant.inc() # pi is irrational, so there will always be another set bit afterwards to cause a round
		return self(False, 1, mant)
	@property
	def e(self):
		e_cont = continued_fraction(
			lambda k: 2 if k == 0 else 2*(k//3+1) if k%3 == 2 else 1,
			lambda k: 1,
			2
		)
		next(e_cont) # 2, all following are [0,1]
		mant = bits([1,0] + [next(e_cont) for i in range(self.mant_len-1)]) # 1+mant_len is the necesary length == 2 + mant_len-1
		if next(e_cont): # just after the end
			mant = mant.inc() # e is irrational, so there will always be another set bit afterwards to cause a round
		return self(False, 1, mant)
	@property
	def phi(self):
		phi_cont = continued_fraction(
			lambda k: 1,
			lambda k: 1,
			2
		)
		mant = bits([next(phi_cont) for i in range(self.mant_len+1)]) # 1+mant_len is the necesary length
		if next(phi_cont): # just after the end
			mant = mant.inc() # phi is irrational, so there will always be another set bit afterwards to cause a round
		return self(False, 0, mant)
	
	
class ArbitraryFloatBase(metaclass=ArbitraryFloatType):
	def __getattr__(self, attr):
		return type(self).__getattribute__(type(self), attr)
	@property
	def exp(self):
		if self.issubnormal:
			return 1 - self.bias
		else: #if self.isnormal:
			return self.exp_bits.decode_int() - self.bias
	@property
	def normalized(self):
		"Returns a 3-tuple of (bool sign, int exp, bits normalized_mant)"
		if self.isinf or self.isnan or self.iszero:
			raise ValueError("Cannot normalize inf, nan, or zero value %r" % self)
		elif self.isnormal:
			return (self.sign, self.exp, [1] + self.mant_bits)
		else: # subnormal
			exp = self.exp - 1
			exp -= self.mant_bits.find(1)
			return (self.sign, exp, self.mant_bits.lstrip())
	@property
	def normalized_exp(self):
		"Returns an int exp"
		if self.isinf or self.isnan or self.iszero:
			raise ValueError("Cannot normalize inf, nan, or zero value %r" % self)
		elif self.isnormal:
			return self.exp
		else: # subnormal
			exp = self.exp - 1
			exp -= self.mant_bits.find(1)
			return exp
	@property
	def normalized_mant(self):
		"Returns a bits normalized_mant"
		if self.isinf or self.isnan or self.iszero:
			raise ValueError("Cannot normalize inf, nan, or zero value %r" % self)
		elif self.isnormal:
			return [1] + self.mant_bits
		else: # subnormal
			return self.mant_bits.lstrip()
	@property
	def ispositive(self):
		return not self.sign and not self.iszero and not self.isnan
	@property
	def isnegative(self):
		return self.sign and not self.iszero and not self.isnan
	@property
	def issignless(self):
		return self.iszero or self.isnan
	@property
	def iszero(self):
		return not any(self.exp_bits) and not any(self.mant_bits)
	@property
	def issubnormal(self):
		return not any(self.exp_bits) and any(self.mant_bits)
	@property
	def isnormal(self):
		return any(self.exp_bits) and not all(self.exp_bits)
	@property
	def isfinite(self):
		return self.isnormal or self.issubnormal or self.iszero
	@property
	def isinf(self):
		return all(self.exp_bits) and not any(self.mant_bits)
	@property
	def isnan(self):
		return all(self.exp_bits) and any(self.mant_bits)
	@property
	def isint(self):
		if self.isinf or self.isnan:
			return False
		elif self.iszero:
			return True
		else:
			sign, exp, mant = self.normalized
			if exp < 0:
				return False
			elif any(mant[exp+1:]):
				return False
			else:
				return True
	def __init_helper__(self, *args):
		"Init self with normalized bool sign, int exp, and bits mant, rounding and over/underflowing accordingly"
		if len(args) == 2:
			sign, exp, mant = False, *args
		else:
			sign, exp, mant = args
		
		if isinstance(sign, bool):
			self.sign = sign
		else:
			raise TypeError("sign must be bool, not %r" % type(sign))
		
		if not isinstance(exp, int):
			raise TypeError("exp must be int, not %r" % type(exp))
		
		if not isinstance(mant, bits):
			raise TypeError("mant must be bits, not %r" % type(mant))
		elif len(mant) < 1 or not mant[0]:
			raise ValueError("mant must be normalized with a leading 1, not %r" % mant)
		
		if self.exp_len == 1: # special case: causes problems b/c all finite values are subnormal, some general assumptions dont hold 
			if exp > 0: # round up to inf
				self.exp_bits = ~bits(self.exp_len)
				self.mant_bits = bits(self.mant_len)
			else:
				mant = bits(-exp) + mant
				if len(mant) <= self.mant_len: # no rounding needed
					self.exp_bits = bits(self.exp_len)
					self.mant_bits = mant.extend(self.mant_len)
				else: # could round
					if all(mant[:self.mant_len+1]): # half-even round up to inf
						self.exp_bits = ~bits(self.exp_len)
						self.mant_bits = bits(self.mant_len)
					elif mant[self.mant_len] and \
						(mant[self.mant_len-1] or  # half-even up 
						 any(mant[self.mant_len+1:])): # round up normally
						self.exp_bits = bits(self.exp_len)
						self.mant_bits = mant[:self.mant_len].inc()
					else: # no rounding up
						self.exp_bits = bits(self.exp_len)
						self.mant_bits = mant[:self.mant_len].inc()
		elif exp > self.max_normal_exp: # round up to inf
			self.exp_bits = ~bits(self.exp_len)
			self.mant_bits = bits(self.mant_len)
		elif exp >= self.min_normal_exp:
			if len(mant)-1 > self.mant_len: # could cause rouding
				if all(mant[1:].crop(self.mant_len+1)): # half-even round up, increasing exp
					if exp+1 > self.max_normal_exp: # round up to inf
						self.exp_bits = ~bits(self.exp_len)
						self.mant_bits = bits(self.mant_len)
					else: # round up to next exp
						self.exp_bits = bits.encode_int(exp+1 + self.bias, self.exp_len)
						self.mant_bits = bits(self.mant_len)
				else:
					self.exp_bits = bits.encode_int(exp + self.bias, self.exp_len)
					if mant[1+self.mant_len] and (mant[self.mant_len] or # half-even up
														 any(mant[self.mant_len+2:])): # other up
						self.mant_bits = mant[1:self.mant_len+1].inc()
					else: # no rounding
						self.mant_bits = mant[1:self.mant_len+1]
			else: # no rounding possible
				self.mant_bits = mant[1:].extend(self.mant_len)
				self.exp_bits = bits.encode_int(exp + self.bias, self.exp_len)
		elif exp >= self.min_subnormal_exp:
			mant = [0]*(self.min_normal_exp - exp - 1) + mant.extend(self.mant_len)
			if exp == self.max_subnormal_exp and all(mant.crop(self.mant_len+1)): # round up, going to normal
				self.exp_bits = bits.encode_int(1, self.exp_len)
				self.mant_bits = bits(self.mant_len)
			else: # stay subnormal, rounding or not
				self.exp_bits = bits(self.exp_len)
				if mant[self.mant_len] and (mant[self.mant_len-1] or # half-even up
											any(mant[self.mant_len+1:])): # other up
					self.mant_bits = mant[:self.mant_len].inc()
				else:
					self.mant_bits = mant[:self.mant_len]
		elif exp == self.min_subnormal_exp-1 and any(mant[1:]): # round up (val is normal/ized)
			self.exp_bits = bits(self.exp_len)
			self.mant_bits = bits.encode_int(1, self.mant_len)
		else: # val.exp < self.min_subnorman_exp, cannot round up
			self.exp_bits = bits(self.exp_len)
			self.mant_bits = bits(self.mant_len)
		
	def __init__(self, *args):
		if len(args) == 1:
			val = args[0]
			if isinstance(val, float):
				val_bits = bits(struct.pack(">d", val))
				sign = val_bits[0]
				exp_bits = val_bits[1:12]
				mant = val_bits[12:]
				if all(exp_bits): # inf or nan double
					self.exp_bits = ~bits(self.exp_len)
					if any(mant): # nan
						self.mant_bits = ~bits(self.mant_len)
					else: # inf
						self.mant_bits = bits(self.mant_len)
					return
				elif not any(exp_bits): # zero or denormal double
					if any(mant): # denormal double
						exp = -1022 - 1 # double bias is 1023, min normal is -1022, we are normalizing
						while not mant[0]: # normalize
							exp -= 1
							mant = mant[1:]
					else: # zero double
						self.exp_bits = bits(self.exp_len)
						self.mant_bits = bits(self.mant_len)
						return
				else: # normal double
					exp = exp_bits.decode_int() - 1023 # double bias is 1023
					mant = [1] + mant
				
				self.__init_helper__(sign, exp, mant)
			elif isinstance(val, int):
				sign = val < 0
				val = abs(val)
				if val == 0:
					self.sign = False
					self.exp_bits = bits(self.exp_len)
					self.mant_bits = bits(self.mant_len)
					return
				exp = val.bit_length() - 1
				mant = bits.encode_int(val)
				self.__init_helper__(sign, exp, mant)
				return
			elif isinstance(val, ArbitraryFloatBase):
				self.sign = val.sign
				# Following this we will assume positive (i.e. "round up" means round down for negatives)
				if val.isinf:
					self.exp_bits = ~bits(self.exp_len)
					self.mant_bits = bits(self.mant_len)
				elif val.isnan:
					self.exp_bits = ~bits(self.exp_len)
					if val.mant_len < self.mant_len or any(val.mant_bits[:self.mant_len]): # copy the NAN signal bits if possible
						self.mant_bits = val.mant_bits.crop(self.mant_len)
					else:
						self.mant_bits = ~bits(self.mant_len)
				elif val.iszero:
					self.exp_bits = bits(self.exp_len)
					self.mant_bits = bits(self.mant_len)
				else: # normal or subnormal
					self.__init_helper__(*val.normalized)
			else:
				raise TypeError(args[0])
		elif 2 <= len(args) <= 3:
			try:
				self.__init_helper__(*args)
				return
			except TypeError: # exp is bits or error
				pass 
			if len(args) == 3 and not isinstance(args[0], bool):
				raise TypeError("sign must be bool object in 2/3-arg %s.__init__" % type(self).__name__)
			if not isinstance(args[-1], bits):
				raise TypeError("mantissa must be bits object in 2/3-arg %s.__init__" % type(self).__name__)
			if not isinstance(args[-2], bits):
				raise TypeError("exp must be int or bits object in 2/3-arg %s.__init__" % type(self).__name__)
			self.sign = bool(args[0]) if len(args) == 3 else False # positive if no sign given
			self.exp_bits = bits(args[-2])
			self.mant_bits = bits(args[-1])
			if len(self.exp_bits) != self.exp_len:
				raise ValueError(args[-2])
			if len(self.mant_bits) != self.mant_len:
				raise ValueError(args[-1])
		elif len(args) == 0:
			self.sign = False
			self.exp_bits = bits(self.exp_len)
			self.mant_bits = bits(self.mant_len)
		else:
			raise TypeError("%s takes 0 to 3 arguments" % type(self).__name__)
	def __int__(self):
		if self.isnan:
			raise ValueError("Cannot convert %s NaN to integer" % type(self).__name__)
		if self.isinf:
			raise OverflowError("Cannot convert %s infinity to integer" % type(self).__name__)
		if self.iszero:
			return 0
		sign, exp, mant = self.normalized
		mant = mant.crop(exp+1)
		return (-1)**self.sign * mant.decode_int()
	def __float__(self):
		if self.exp_len == 11 and self.mant_len == 52: # IEEE 754 double
			return struct.unpack(">d", bytes([self.sign] + self.exp_bits + self.mant_bits))[0]
		return float(ArbitraryFloatType(11, 52)(self))
	def __bits__(self):
		return [self.sign] + self.exp_bits + self.mant_bits
	def hex(self):
		s = '-'*self.sign
		if self.isinf:
			return s + 'inf'
		elif self.isnan:
			return 'nan'
		elif self.iszero:
			return s + '0x0.0p+0'
		elif self.isnormal:
			s += '0x1.'
		else: # self.issubnormal
			s += '0x0.'
		s += self.mant_bits.hex()
		s += 'p'
		s += '%+d' % self.exp
		return s
	@classmethod
	def fromhex(cls, hexstr):
		import re
		hex_re = re.compile(
			"""
			(^
				([+-])?
				(0x)?
				([0-9a-fA-F]*)
				(.
					([0-9a-fA-F]*)
				)?
				(p
					([+-][0-9]+)
				)?
			$)|(^
				([+-]?)(inf|nan)
			$)
			""",
			re.X
		)
		match = hex_re.match(hexstr)
		if not match:
			raise ValueError(hexstr)
		groups = match.groups()
		if groups[0] is not None: # not inf or nan
			sign = groups[1] == '-'
			baseexp = 4*len(groups[3] or '') - 1 # -1 for the leading bit not in the fractional part
			mant = bits.fromhex(groups[3] or '') + bits.fromhex(groups[5] or '')
			if not any(mant):
				return cls(sign, bits(cls.exp_len), bits(cls.mant_len))
			baseexp -= mant.find(1) # wont return -1 because it has a 1
			mant = mant.lstrip()
			exp = baseexp + int(groups[7] or 0)
			return cls(sign, exp, mant)
		else: # if groups[8] is not None: # inf or nan
			if groups[10] == 'inf':
				if groups[9] == '-':
					return -cls.inf
				else:
					return cls.inf
			else:
				if groups[9] == '-':
					return -cls.nan
				else:
					return cls.nan
	
	def totalOrder(self, other):
		"""
		Is self ordered before other or the same as other (-0 is not the same as +0 for this purpose)?
		-NaN < -inf < -finite < -0 < +0 < +finite < +inf < +NaN
		"""
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		
		if self.isnan or other.isnan:
			if self.isnan and self.sign:
				if other.isnan and other.sign:
					return self.mant_bits >= other.mant_bits # IEEE 754 5.10.d.3.iii
				else:
					return True # IEEE 754 5.10.d.1 and IEEE 754 5.10.d.3.i
			elif self.isnan:
				if other.isnan and not other.sign:
					return self.mant_bits <= other.mant_bits # IEEE 754 5.10.d.3.iii
				else:
					return False # IEEE 754 5.10.d.1 and IEEE 754 5.10.d.3.i
			elif other.isnan and other.sign:
				return False # self is not NaN, everything else > -NaN
			elif other.isnan:
				return True # self is not NaN, everything else <= +Nan
		# No NaN below here
		elif self.isinf and self.sign:
			return True # -inf <= everything other than -NaN
		elif self.isinf:
			if other.isinf and not other.sign:
				return True # +inf <= +inf
			else:
				return False # +inf > everything other than +NaN and +inf
		elif other.isinf and other.sign:
			return False # everything other than -Nan and -inf > -inf
		elif other.isinf:
			return True # everything other than +Nan and +inf <= +inf
		# No NaN or inf below here
		elif self.iszero:
			if other.iszero and (self.sign or not other.sign):
				return True # -0 <= +0 and -0 <= -0 and +0 <= +0
			elif other.iszero:
				return False # +0 > -0
			elif other.sign:
				return False # 0 > -finite
			else:
				return True # 0 <= +finite
		elif other.iszero:
			if self.sign:
				return True # -finite <= 0
			else:
				return False # +finite > 0
		# Only finite below here
		s_sign, s_exp, s_mant = self.normalized
		o_sign, o_exp, o_mant = other.normalized
		if len(s_mant) > len(o_mant):
			o_mant = o_mant.extend(len(s_mant))
		else:
			s_mant = s_mant.extend(len(o_mant))
		if s_sign and not o_sign:
			return True #-anything < +anything
		elif o_sign and not s_sign:
			return False # +anything > -anything
		# Signs are the same
		if s_sign: # swap and abs() if negative to simplify logic
			s_sign, o_sign = False, False
			s_exp, o_exp = o_exp, s_exp
			s_mant, o_mant = o_mant, s_mant
		# Both are positive
		if s_exp > o_exp:
			return False
		elif s_exp < o_exp:
			return True
		elif s_mant > o_mant:
			return False
		else:
			return True
	def totalOrderMag(self, other):
		return abs(self).totalOrder(abs(other))
	
	
	
	def inc(self):
		# TODO: count>1
		if self.isnan:
			return self
		elif self.ispositive:
			if self.isinf:
				return self
			if all(self.mant_bits):
				if self.exp == self.max_normal_exp:
					return self.inf
				return type(self)(False, self.exp_bits.inc(), bits(self.mant_len))
			return type(self)(False, self.exp_bits, self.mant_bits.inc())
		elif self.isnegative:
			if not any(self.mant_bits):
				if self.exp == self.min_normal_exp:
					return type(self)(True, bits(self.exp_len), ~bits(self.mant_len))
				return type(self)(True, self.exp_bits.dec(), ~bits(self.mant_len))
			return type(self)(True, self.exp_bits, self.mant_bits.dec())
		else: # zero
			return type(self)(False, bits(self.exp_len), bits.encode_int(1, self.mant_len))
	def dec(self):
		# TODO
		if self.isnan:
			return self
		elif self.isnegative:
			if self.isinf:
				return self
			if all(self.mant_bits):
				if self.exp == self.max_normal_exp:
					return self.inf
				return type(self)(True, self.exp_bits.inc(), bits(self.mant_len))
			return type(self)(True, self.exp_bits, self.mant_bits.inc())
		elif self.ispositive:
			if not any(self.mant_bits):
				if self.exp == self.min_normal_exp:
					return type(self)(False, bits(self.exp_len), ~bits(self.mant_len))
				return type(self)(False, self.exp_bits.dec(), ~bits(self.mant_len))
			return type(self)(False, self.exp_bits, self.mant_bits.dec())
		else: # zero
			return type(self)(True, bits(self.exp_len), bits.encode_int(1, self.mant_len))
	def __neg__(self):
		return type(self)(not self.sign, self.exp_bits[:], self.mant_bits[:])
	def __pos__(self):
		return self
	def __abs__(self):
		return type(self)(False, self.exp_bits[:], self.mant_bits[:])
		
	def __mul_helper__(self, other):
		if self.isnan or other.isnan or \
			(self.isinf and other.iszero) or \
			(self.iszero and other.isinf):
				return ArbitraryFloatType.best_precision(self, other).nan
		elif self.iszero or other.iszero:
			if self.sign ^ other.sign:
				return -ArbitraryFloatType.best_precision(self, other).zero
			else:
				return +ArbitraryFloatType.best_precision(self, other).zero
		elif self.isinf or other.isinf:
			if self.sign ^ other.sign:
				return -ArbitraryFloatType.best_precision(self, other).inf
			else:
				return +ArbitraryFloatType.best_precision(self, other).inf
		else: # both are finite
			return None
	def __mul__(self, other):
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__mul_helper__(other)
		if ret is not None: # Special case happened
			return ret
		# Both finite
		s_sign, s_exp, s_mant = self.normalized
		o_sign, o_exp, o_mant = other.normalized
		sign = s_sign ^ o_sign
		exp = s_exp + o_exp
		mant = s_mant.multiply_unsigned(o_mant)
		if mant[0]:
			exp += 1
		else:
			mant = mant[1:]
		return ArbitraryFloatType.best_precision(self, other)(sign, exp, mant)
	def mul(self, other):
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__mul_helper__(other)
		if ret is not None: # Special case happened
			return ret
		# Both finite
		s_sign, s_exp, s_mant = self.normalized
		o_sign, o_exp, o_mant = other.normalized
		sign = s_sign ^ o_sign
		exp = s_exp + o_exp
		mant = s_mant.multiply_unsigned(o_mant)
		if mant[0]:
			exp += 1
		else:
			mant = mant[1:]
		return ArbitraryFloatType.least_precision(exp, mant)(sign, exp, mant)
		
	def __add_helper__(self, other):
		if self.isnan or other.isnan or \
			(self.isinf and other.isinf and (self.sign ^ other.sign)):
			return ArbitraryFloatType.best_precision(self, other).nan
		elif self.isinf:
			return ArbitraryFloatType.best_precision(self, other)(self)
		elif other.isinf:
			return ArbitraryFloatType.best_precision(self, other)(other)
		elif self.iszero:
			return ArbitraryFloatType.best_precision(self, other)(other)
		elif other.iszero:
			return ArbitraryFloatType.best_precision(self, other)(self)
		else: # both finite
			return None
	def __add__(self, other):
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__add_helper__(other)
		if ret is not None: # Speical case happened
			return ret
		# Both finite
		if self.sign ^ other.sign:
			return self - -other
		else: # same sign
			sign, s_exp, s_mant = self.normalized
			_   , o_exp, o_mant = other.normalized
			if s_exp >= o_exp:
				exp = s_exp
				mant = s_mant.add_unsigned_ljust(bits(s_exp - o_exp) + o_mant)
			else:
				exp = o_exp
				mant = o_mant.add_unsigned_ljust(bits(o_exp - s_exp) + s_mant)

			if mant[0]:
				exp += 1
			else:
				mant = mant[1:]
			return ArbitraryFloatType.best_precision(self, other)(sign, exp, mant)
	def add(self, other):
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__add_helper__(other)
		if ret is not None: # Speical case happened
			return ret
		# Both finite
		if self.sign ^ other.sign:
			return self.sub(-other)
		else: # same sign
			sign, s_exp, s_mant = self.normalized
			_   , o_exp, o_mant = other.normalized
			if s_exp >= o_exp:
				exp = s_exp
				mant = s_mant.add_unsigned_ljust(bits(s_exp - o_exp) + o_mant)
			else:
				exp = o_exp
				mant = o_mant.add_unsigned_ljust(bits(o_exp - s_exp) + s_mant)

			if mant[0]:
				exp += 1
			else:
				mant = mant[1:]
			return ArbitraryFloatType.least_precision(exp, mant)(sign, exp, mant)
	def __sub_helper__(self, other):
		if self.isnan or other.isnan or \
			(self.isinf and other.isinf and not (self.sign ^ other.sign)):
			return ArbitraryFloatType.best_precision(self, other).nan
		elif self.isinf:
			return ArbitraryFloatType.best_precision(self, other)(self)
		elif other.isinf:
			return ArbitraryFloatType.best_precision(self, other)(-other)
		elif self.iszero:
			return ArbitraryFloatType.best_precision(self, other)(-other)
		elif other.iszero:
			return ArbitraryFloatType.best_precision(self, other)(self)
		else: # both finite
			return None
	def __sub__(self, other):
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__sub_helper__(other)
		if ret is not None: # Special case happened
			return ret
		# Both finite
		if self.sign ^ other.sign: # -5 - +6 == -5 + -6), +5 - -6 == +5 + +6
			return self + -other
		else: # same sign
			sign, s_exp, s_mant = self.normalized
			_   , o_exp, o_mant = other.normalized

			if s_exp > o_exp:
				exp = s_exp
				mant = s_mant.sub_unsigned_ljust(bits(s_exp - o_exp) + o_mant)
			elif s_exp < o_exp:
				sign = not sign
				exp = o_exp
				mant = o_mant.sub_unsigned_ljust(bits(o_exp - s_exp) + s_mant)
			elif s_mant > o_mant:
				exp = s_exp
				mant = s_mant.sub_unsigned_ljust(o_mant)
			elif s_mant < o_mant:
				sign = not sign
				exp = o_exp
				mant = o_mant.sub_unsigned_ljust(s_mant)
			else: # equal inputs, a-a == 0
				return ArbitraryFloatType.best_precision(self, other).zero

			if any(mant):
				exp -= mant.find(1)
				mant = mant.lstrip()
				return ArbitraryFloatType.best_precision(self, other)(sign, exp, mant)
			else: # zero
				return ArbitraryFloatType.best_precision(self, other).zero
	def sub(self, other):
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__sub_helper__(other)
		if ret is not None: # Special case happened
			return ret
		# Both finite
		if self.sign ^ other.sign: # -5 - +6 == -5 + (-6), +5 - -6 == +5 + +6
			return self.add(-other)
		else: # same sign
			sign, s_exp, s_mant = self.normalized
			_   , o_exp, o_mant = other.normalized

			if s_exp > o_exp:
				exp = s_exp
				mant = s_mant.sub_unsigned_ljust(bits(s_exp - o_exp) + o_mant)
			elif s_exp < o_exp:
				sign = not sign
				exp = o_exp
				mant = o_mant.sub_unsigned_ljust(bits(o_exp - s_exp) + s_mant)
			elif s_mant > o_mant:
				exp = s_exp
				mant = s_mant.sub_unsigned_ljust(o_mant)
			elif s_mant < o_mant:
				sign = not sign
				exp = o_exp
				mant = o_mant.sub_unsigned_ljust(s_mant)
			else: # equal inputs, a-a == 0
				return ArbitraryFloatType.best_precision(self, other).zero

			if any(mant):
				exp -= mant.find(1)
				mant = mant.lstrip()
				return ArbitraryFloatType.least_precision(exp, mant)(sign, exp, mant)
			else: # zero
				return ArbitraryFloatType.best_precision(self, other).zero
	
	def __truediv__(self, other):
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		if other.iszero:
			raise ZeroDivisionError("ArbitraryFloat division by zero")
		elif self.isnan or other.isnan or \
			(self.isinf and other.isinf):
			return ArbitraryFloatType.best_precision(self, other).nan
		elif other.isinf or self.iszero:
			if self.sign == other.sign:
				return ArbitraryFloatType.best_precision(self, other).zero
			else:
				return -ArbitraryFloatType.best_precision(self, other).zero
		elif self.isinf:
			if self.sign == other.sign:
				return ArbitraryFloatType.best_precision(self, other).inf
			else:
				return -ArbitraryFloatType.best_precision(self, other).inf
		else:
			def round_undecided(mant, length):
				"""Return true if an additional 1 added to mant could cause a change when rounding to length+1 bits (+1 for the leading normalized bit)
				Only used when len(mant) > length + 2"""
				if not mant[length+1]: #the first bit after used, if not set, cannot round up
					return False
				if mant[length]: # half-even, will round up
					return False
				if any(mant[length+2:]): # round up
					return False
				return True # could still change
			s_sign, s_exp, s_mant = self.normalized
			o_sign, o_exp, o_mant = other.normalized
			sign = s_sign ^ o_sign
			exp = s_exp - o_exp
			length = 1+max(self.mant_len, other.mant_len)
			mant = bits()
			s_mant = s_mant.extend(length)
			o_mant = o_mant.extend(length)
			s_mant = s_mant.decode_int()
			o_mant = o_mant.decode_int()
			while s_mant and (len(mant) <= length+1 or round_undecided(mant, length-1)):
				if s_mant >= o_mant:
					mant += [1]
					s_mant -= o_mant
				else:
					mant += [0]
				s_mant *= 2
			exp -= mant.find(1)
			mant = mant.lstrip()
			return ArbitraryFloatType.best_precision(self, other)(sign, exp, mant)
	def __rtruediv__(self, other):
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		return other.__truediv__(self)
	
	def __pow_helper__(self, other):
		raise TODO
		# shoule implement these checks in the way specified by man pow(3)
		if (self.normalized_exp == 0 and self.isint) or other.iszero: # self == 1 or other == 0
			return ArbitraryFloatType.best_precision(self, other)(1)
		elif self.isnan:
			return ArbitraryFloatType.best_precision(self, other).nan
		elif other.isnan:
			if self.normalized_exp == 0 and self.isint: # self is 1, 1**nan is 1 for some reason
				raise TODO
			return ArbitraryFloatType.best_precision(self, other).nan
		else:
			raise TODO
	def __pow__(self, other):
		if isinstance(other, int) or (isinstance(other, ArbitraryFloatBase) and other.isint) or (isinstance(other, float) and other.is_integer()):
			other = int(other)
			ret = type(self)(1)
			while other > 0:
				ret *= self
				other -= 1
			while other < 0:
				ret /= self
				other += 1
			return ret
		elif not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__pow_helper__(other)
		if ret is not None: # Special case happened
			return ret
		raise TODO
	def pow(self, other):
		if isinstance(other, int) or (isinstance(other, ArbitraryFloatBase) and other.isint) or (isinstance(other, float) and other.is_integer()):
			other = int(other)
			ret = ArbitraryFloatType(1, 1)(1)
			while other > 0:
				ret = ret.mul(self)
				other -= 1
			if other < 0:
				inv = ArbitraryFloatType(self.exp_len+1, self.mant_len**2)(1)/self
				while other < 0:
					ret = ret.mul(inv)
					other += 1
			return ret
		elif not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__pow_helper__(other)
		if ret is not None: # Special case happened
			return ret
		raise TODO
		
	def __and_helper__(self, other):
		if self.isnan or other.isnan or \
			self.isinf or other.isinf:
			return ArbitraryFloatType.best_precision(self, other).nan
		elif self.iszero:
			return abs(ArbitraryFloatType.best_precision(self, other)(self))
		elif other.iszero:
			return abs(ArbitraryFloatType.best_precision(self, other)(other))
		else: # both finite
			return None
	def __and__(self, other):
		"""
		Gives the positive result of binary ANDing the theoretical infinite mantissas of abs(self) and abs(other)
		nan & x == nan
		inf & x == nan
		0 & x == 0
		5.5 & 3.7 == 1.5
			101.1(0) & 11.10110011(0011) == 001.1(0)
		"""
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__and_helper__(other)
		if ret is not None:
			return ret
		
		_, s_exp, s_mant = self.normalized
		_, o_exp, o_mant = other.normalized

		if s_exp == o_exp: # normalized values with same exp, leading 1 will stay
			mant = s_mant.extend(max(len(s_mant), len(o_mant))) & o_mant.extend(max(len(s_mant), len(o_mant)))
			return ArbitraryFloatType.best_precision(self, other)(False, s_exp, mant)
		elif s_exp > o_exp:
			exp = s_exp
			exp_diff = s_exp - o_exp
			big_mant, small_mant = s_mant, o_mant
		else:
			exp = o_exp
			exp_diff = o_exp - s_exp
			big_mant, small_mant = o_mant, s_mant
		small_mant = bits(exp_diff) + small_mant
		mant = big_mant.extend(max(len(big_mant), len(small_mant))) & small_mant.extend(max(len(big_mant), len(small_mant)))
		if not any(mant):
			return ArbitraryFloatType.best_precision(self, other).zero
		exp -= mant.find(1)
		mant = mant.lstrip()
		return ArbitraryFloatType.best_precision(self, other)(False, exp, mant)
	def and_(self, other):
		"""
		Gives the positive result of binary ANDing the theoretical infinite mantissas of abs(self) and abs(other)
		nan & x == nan
		inf & x == nan
		0 & x == 0
		5.5 & 3.7 == 1.5
			101.1(0) & 11.10110011(0011) == 001.1(0)
		"""
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__and_helper__(other)
		if ret is not None:
			return ret
		
		_, s_exp, s_mant = self.normalized
		_, o_exp, o_mant = other.normalized

		if s_exp == o_exp: # normalized values with same exp, leading 1 will stay
			mant = s_mant.extend(max(len(s_mant), len(o_mant))) & o_mant.extend(max(len(s_mant), len(o_mant)))
			return ArbitraryFloatType.least_precision(exp, mant)(False, s_exp, mant)
		elif s_exp > o_exp:
			exp = s_exp
			exp_diff = s_exp - o_exp
			big_mant, small_mant = s_mant, o_mant
		else:
			exp = o_exp
			exp_diff = o_exp - s_exp
			big_mant, small_mant = o_mant, s_mant
		small_mant = bits(exp_diff) + small_mant
		mant = big_mant.extend(max(len(big_mant), len(small_mant))) & small_mant.extend(max(len(big_mant), len(small_mant)))
		if not any(mant):
			return ArbitraryFloatType.best_precision(self, other).zero
		exp -= mant.find(1)
		mant = mant.lstrip()
		return ArbitraryFloatType.least_precision(exp, mant)(False, exp, mant)
	
	def __or_helper__(self, other):
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		if self.isnan or other.isnan or \
			self.isinf or other.isinf:
			return self.nan
		elif self.iszero:
			return abs(ArbitraryFloatType.best_precision(self, other)(other))
		elif other.iszero:
			return abs(ArbitraryFloatType.best_precision(self, other)(self))
		else: # both finite
			return None
	def __or__(self, other):
		"""
		Gives the positive result of binary ORing the theoretical infinite mantissas of abs(self) and abs(other)
		nan | x == nan
		inf | x == nan
		0 | x == x
		5.5 & 3.7 == 7.7
			101.1(0) | 11.1011(0011) == 111.1011(0011)
		"""
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__or_helper__(other)
		if ret is not None:
			return ret
		
		_, s_exp, s_mant = self.normalized
		_, o_exp, o_mant = other.normalized

		if s_exp >= o_exp:
			exp = s_exp
			exp_diff = s_exp - o_exp
			big_mant, small_mant = s_mant, o_mant
		else:
			exp = o_exp
			exp_diff = o_exp - s_exp
			big_mant, small_mant = o_mant, s_mant
		small_mant = bits(exp_diff) + small_mant
		mant = big_mant.extend(max(len(big_mant), len(small_mant))) | small_mant.extend(max(len(big_mant), len(small_mant)))
		# don't need to re-normalize since ORing cant remove the leading 1
		return ArbitraryFloatType.best_precision(self, other)(False, exp, mant)
	def or_(self, other):
		"""
		Gives the positive result of binary ORing the theoretical infinite mantissas of abs(self) and abs(other)
		nan | x == nan
		inf | x == nan
		0 | x == x
		5.5 & 3.7 == 7.7
			101.1(0) | 11.1011(0011) == 111.1011(0011)
		"""
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__or_helper__(other)
		if ret is not None:
			return ret
		
		_, s_exp, s_mant = self.normalized
		_, o_exp, o_mant = other.normalized

		if s_exp >= o_exp:
			exp = s_exp
			exp_diff = s_exp - o_exp
			big_mant, small_mant = s_mant, o_mant
		else:
			exp = o_exp
			exp_diff = o_exp - s_exp
			big_mant, small_mant = o_mant, s_mant
		small_mant = bits(exp_diff) + small_mant
		mant = big_mant.extend(max(len(big_mant), len(small_mant))) | small_mant.extend(max(len(big_mant), len(small_mant)))
		# don't need to re-normalize since ORing cant remove the leading 1
		return ArbitraryFloatType.least_precision(exp, mant)(False, exp, mant)
	
	def __xor_helper__(self, other):
		if self.isnan or other.isnan or \
			self.isinf or other.isinf:
			return self.nan
		elif self.iszero:
			return abs(ArbitraryFloatType.best_precision(self, other)(other))
		elif other.iszero:
			return abs(ArbitraryFloatType.best_precision(self, other)(self))
		else: # both finite
			return None
	def __xor__(self, other):
		"""
		Gives the positive result of binary XORing the theoretical infinite mantissas of abs(self) and abs(other)
		nan ^ x == nan
		inf ^ x == nan
		0 ^ x == x
		5.5 ^ 3.7 == 6.2
			101.1(0) & 11.10(1100) == 110.00(1100)
		"""
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__xor_helper__(other)
		if ret is not None:
			return ret
		
		_, s_exp, s_mant = self.normalized
		_, o_exp, o_mant = other.normalized

		if s_exp == o_exp: # normalized values with same exp, leading 1 will not stay
			mant = s_mant.extend(max(len(s_mant), len(o_mant))) ^ o_mant.extend(max(len(s_mant), len(o_mant)))
			if not any(mant):
				return ArbitraryFloatType.best_precision(self, other).zero
			exp = s_exp - mant.find(1)
			mant = mant.lstrip()
			return ArbitraryFloatType.best_precision(self, other)(False, exp, mant)
		elif s_exp > o_exp: # leading 1 will stay
			exp = s_exp
			exp_diff = s_exp - o_exp
			big_mant, small_mant = s_mant, o_mant
		else: # leading 1 will stay
			exp = o_exp
			exp_diff = o_exp - s_exp
			big_mant, small_mant = o_mant, s_mant
		small_mant = bits(exp_diff) + small_mant
		mant = big_mant.extend(max(len(big_mant), len(small_mant))) ^ small_mant.extend(max(len(big_mant), len(small_mant)))
		# don't need to re-normalize since XORing cant remove the leading 1 from unaligned mantissas
		return ArbitraryFloatType.best_precision(self, other)(False, exp, mant)
	def xor(self, other):
		"""
		Gives the positive result of binary XORing the theoretical infinite mantissas of abs(self) and abs(other)
		nan ^ x == nan
		inf ^ x == nan
		0 ^ x == x
		5.5 ^ 3.7 == 6.2
			101.1(0) & 11.10(1100) == 110.00(1100)
		"""
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		ret = self.__xor_helper__(other)
		if ret is not None:
			return ret
		
		_, s_exp, s_mant = self.normalized
		_, o_exp, o_mant = other.normalized

		if s_exp == o_exp: # normalized values with same exp, leading 1 will not stay
			mant = s_mant.extend(max(len(s_mant), len(o_mant))) ^ o_mant.extend(max(len(s_mant), len(o_mant)))
			if not any(mant):
				return ArbitraryFloatType.best_precision(self, other).zero
			exp = s_exp - mant.find(1)
			mant = mant.lstrip()
			return ArbitraryFloatType.least_precision(exp, mant)(False, exp, mant)
		elif s_exp > o_exp: # leading 1 will stay
			exp = s_exp
			exp_diff = s_exp - o_exp
			big_mant, small_mant = s_mant, o_mant
		else: # leading 1 will stay
			exp = o_exp
			exp_diff = o_exp - s_exp
			big_mant, small_mant = o_mant, s_mant
		small_mant = bits(exp_diff) + small_mant
		mant = big_mant.extend(max(len(big_mant), len(small_mant))) ^ small_mant.extend(max(len(big_mant), len(small_mant)))
		# don't need to re-normalize since XORing cant remove the leading 1 from unaligned mantissas
		return ArbitraryFloatType.least_precision(exp, mant)(False, exp, mant)

	#in another function (to retur nthe shortest unique repr, compare digits to self.inc and self.dec
	
	def _evenbase_str(self, base):
		if base % 2:
			raise ValueError(base,"is not a multiple of 2")
		elif not base:
			raise ValueError("Cannot convert to base zero")
		elif base > len(base_module.digits):
			raise ValueError("Base %d is too large" % base)
			
		if self.isinf:
			return "-"*self.sign + "inf"
		elif self.isnan:
			return "nan"
		elif self.iszero:
			return "-"*self.sign + "0.0"
		
		
		
		sign, exp, mant = self.normalized
		
		num = mant.decode_int()
		effective_exp = len(mant) - exp - 1
		print(num, effective_exp, *self.normalized)
		if effective_exp >= 0:
			num *= (base//2)**effective_exp
			num_str = base_module.to_positional_base(num, base)
			print(num_str)
			num_str = num_str.rjust(effective_exp+1, "0")
			print(num_str)
			if exp >= 0:
				return "-"*sign + num_str[:-effective_exp] + "." + num_str[-effective_exp:]
			else:
				return "-"*sign + num_str[:-effective_exp] + "." + num_str[-effective_exp:]
		else:
			num *= 2**-effective_exp
			return "-"*sign + base_module.to_positional_base(num, base) + ".0"
		
		print(num)
		
		#decimal_num = mant.decode_int()
		#for i in range(len(mant)-1):
		#	decimal_num *= 5
		#exp -= len(mant)-1
		#decimal_str = str(decimal_num)
		#decimal_str = "-"*sign + decimal_str[:-exp] + "." + decimal_str[-exp:]
		#return decimal_str
	def __repr__(self):
		try:
			if not color:
				return type(self).__name__ + "(%s) # %r" % (bits(self), float(self))
			else:
				return type(self).__name__ + "(\x1b[31m%s\x1b[32m%s\x1b[33m%s\x1b[0m) # %r" % (bits([self.sign]), self.exp_bits, self.mant_bits, float(self))
		except Exception as e:
			import traceback
			if not color:
				return type(self).__name__ + "(%s) # error \n\t%s" % (bits(self), '\n\t'.join(traceback.format_exc().split('\n')))
			else:
				return type(self).__name__ + "(\x1b[31m%s\x1b[32m%s\x1b[33m%s\x1b[0m) # error \n\t%s" % (bits([self.sign]), self.exp_bits, self.mant_bits, '\n\t'.join(traceback.format_exc().split('\n')))
	
	def sin(self, terms=4):
		"Inaccurate"
		ret = self
		for i in range(1,terms):
			ret += (-1)**i * self.pow(2*i+1) / math.factorial(2*i+1)
		return ret
	
	# Commutative right-side operators
	__radd__ = __add__
	__rmul__ = __mul__
	__rand__ = __and__
	__ror__ = __or__
	__rxor__ = __xor__
	
	# Anticommutative right-side operators
	__rsub__ = lambda s, o: -s + o
	

Quarter = ArbitraryFloatType(4, 3)
Half = ArbitraryFloatType(5, 10)
Single = ArbitraryFloatType(8, 23)
Double = ArbitraryFloatType(11, 52)
Quadruple = ArbitraryFloatType(15, 112)
Octuple = ArbitraryFloatType(19, 236)

class LongDouble(ArbitraryFloatBase):
	exp_len = 15
	mant_len = 63
	def __bits__(self):
		return [self.sign] + self.exp_bits + [self.isnormal] + self.mant_bits
	def __repr__(self):
		try:
			if not color:
				return type(self).__name__ + "(%s) # %r" % (bits(self), float(self))
			else:
				return type(self).__name__ + "(\x1b[31m%s\x1b[32m%s\x1b[34m%s\x1b[33m%s\x1b[0m) # %r" % (bits([self.sign]), self.exp_bits, bits([self.isnormal]), self.mant_bits, float(self))
		except Exception as e:
			import traceback
			if not color:
				return type(self).__name__ + "(%s) # error \n\t%s" % (bits(self), '\n\t'.join(traceback.format_exc().split('\n')))
			else:
				return type(self).__name__ + "(\x1b[31m%s\x1b[32m%s\x1b[34m%s\x1b[33m%s\x1b[0m) # error \n\t%s" % (bits([self.sign]), self.exp_bits, bits([self.isnormal]), self.mant_bits, '\n\t'.join(traceback.format_exc().split('\n')))
		
