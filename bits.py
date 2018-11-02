from bitarray import bitarray
import struct

class bits:
	"Immutable(ish) layer on the bitarray module"
	def __init__(self, val=0):
		if hasattr(val, "__bits__"):
			self._val = val.__bits__()._val
		elif isinstance(val, bytes):
			self._val = bitarray()
			self._val.frombytes(val)
		elif isinstance(val, int):
			self._val = bitarray(val)
			self._val.setall(0)
		else:
			self._val = bitarray(val)
	def __getitem__(self, index):
		val = self._val[index]
		if isinstance(val, bitarray):
			return bits(self._val[index])
		return val
	def __str__(self):
		return str(self._val).replace("bitarray('","").replace("')","")
		#return str(self._val).replace('array', 's')
	def __repr__(self):
		return repr(self._val).replace('array', 's')
	def __len__(self):
		return len(self._val)
	def __eq__(self, other):
		if isinstance(other, bits):
			return other._val == self._val
		return bitarray(other) == self._val
	def __lt__(self, other):
		if isinstance(other, bits):
			return other._val > self._val #Todo, swap these for readability
		return bitarray(other) > self._val
	def __gt__(self, other):
		if isinstance(other, bits):
			return other._val < self._val #Todo, swap these for readability
		return bitarray(other) < self._val
		
	def __add__(self, other):
		if isinstance(other, bits):
			return bits(self._val + other._val)
		return bits(self._val + other)
	def __radd__(self, other):
		if isinstance(other, bits):
			return bits(other._val + self._val)
		return bits(bitarray(other) + self._val)
	def __or__(self, other):
		if isinstance(other, bits):
			return bits(self._val | other._val)
		return bits(self._val | other)
	def __ror__(self, other):
		if isinstance(other, bits):
			return bits(other._val | self._val)
		return bits(other | self._val)
	def __xor__(self, other):
		if isinstance(other, bits):
			return bits(self._val ^ other._val)
		return bits(self._val ^ other)
	def __rxor__(self, other):
		if isinstance(other, bits):
			return bits(other._val ^ self._val)
		return bits(other ^ self._val)
	def __invert__(self):
		return bits(~self._val)
	
	def __bytes__(self):
		return self._val.tobytes()
	def cut(self, length):
		for i in range(0,len(self),length):
			yield self[i:length+i].extend(length)
	def hex(self):
		s = ''
		for nibble in self.cut(4):
			s += '0123456789abcdef'[nibble.decode_int()]
		return s
	@classmethod
	def fromhex(cls, hexstr):
		if hexstr:
			return cls.encode_int(int(hexstr, 16), 4*len(hexstr))
		return cls()
	def find(self, sub, start=None, end=None):
		if start is not None or end is not None:
			return start + self[start:end].find(sub)
		if not isinstance(sub, bits):
			sub = bits(sub) if sub not in [0,1] else bits([sub])
		for i in range(len(self)):
			if self[i:len(sub)+i] == sub:
				return i
		return -1
	def strip(self, val=0):
		return self[self.find(not val):~self[::-1].find(not val)]
	def rstrip(self, val=0):
		return self[:~self[::-1].find(not val)]
	def lstrip(self, val=0):
		return self[self.find(not val):]
		
	@classmethod
	def encode_int(cls, val, length=None, signed=False):
		if not signed and val < 0:
			raise ValueError("Cannot encode %d as unsigned" % val)
		if length is None:
			if not signed:
				length = val.bit_length()
			elif signed and val >= 0:
				length = val.bit_length() + 1
			else:
				length = (~val).bit_length() + 1
		if signed and (val < -2**(length-1) or val >= 2**(length-1)) or \
			not signed and val >= 2**length:
			raise ValueError("Cannot encode %d as a %d-bit %s integer" % (val, length, "signed" if signed else "unsigned"))
		val += 2**(length+1)
		# set the negative bit and adjust 2s-complement for signed
		# leave unsigned alone
		# (after we chop off the leading bits)
		byteval = val.to_bytes(length//8+2, 'big')
		bitval = bitarray()
		bitval.frombytes(byteval)
		del bitval[:-length] # delete leading bits
		return bits(bitval)
	def decode_int(self, signed=False):
		if len(self) % 8:
			val = (bits((-len(self))%8) + self).decode_int()
		else:
			val = int.from_bytes(self._val.tobytes(), 'big')
		if signed and self[0]:
			val = ~val
		return val
	
	def inc(self, count=1):
		length = len(self)
		val = self.decode_int() + count
		return self.encode_int(val % 2**length, length)
	def dec(self, count=1):
		length = len(self)
		val = self.decode_int() - count
		return self.encode_int(val % 2**length, length)
	
	def extend(self, length):
		"Extend to at least length bits"
		if len(self) > length:
			return self
		return self + bits(length-len(self))
	def crop(self, length):
		"Extend or truncate, depending on length"
		if length < 0:
			return bits()
		if len(self) > length:
			return self[:length]
		return self + bits(length-len(self))
	def truncate(self, length):
		if length < 0:
			return bits()
		"Truncate to at most length bits (same as self[:length])"
		return self[:length]
	
	def multiply_unsigned(self, other):
		"Multiply self and other, returning a bits object of length len(self)+len(other)"
		if not isinstance(other, bits):
			raise TypeError(other)
		return bits.encode_int(self.decode_int()*other.decode_int(), len(self)+len(other))
		
