from bits import bits
import struct

color = True

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
	def least_precision(val):
		"Return a type that can represent val exactly\nAlways ArbitraryFloatType(11, 52) for floats\nFor ints, returns a type that can represent every integer in [0,val]"
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
	def isinf(self):
		return all(self.exp_bits) and not any(self.mant_bits)
	@property
	def isnan(self):
		return all(self.exp_bits) and any(self.mant_bits)
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
		
		if exp > self.max_normal_exp: # round up to inf
			self.exp_bits = ~bits(self.exp_len)
			self.mant_bits = bits(self.mant_len)
		elif exp > self.min_normal_exp:
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
	
	def __mul__(self, other):
		if not isinstance(other, ArbitraryFloatBase):
			other = ArbitraryFloatType.least_precision(other)(other)
		if self.isnan or other.isnan or \
			(self.isinf and other.iszero) or \
			(self.iszero and other.isinf):
				return self.nan
		elif self.iszero or other.iszero:
			return -self.zero if (self.sign ^ other.sign) else self.zero
		elif self.isinf or other.isinf:
			return -self.inf if (self.sign ^ other.sign) else self.inf
		else: # both are finite
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
			
	
	def _decimal_str(self):
		pass
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


def multiplication_helper(mant1: bits, mant2: bits) -> (bits, int):
	"Takes in two normalized mantissas and returnsa normalized mantissa and an exp adjustment (0 if not needed, 1 if carried)"
	out = mant1.multiply_unsigned(mant2)
	if not out[0]:
		return out[1:], 0
	else:
		return out, 1
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
		