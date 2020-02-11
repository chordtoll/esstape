def byte2bits(b):
    return [1==(b>>i)&0x01 for i in range(8)][::-1]

def chr2bits(c):
    return [1==(ord(c)>>i)&0x01 for i in range(8)][::-1]

def short2bits(s):
    return [1==(s>>i)&0x01 for i in range(16)][::-1]

def flip16(b):
    for i in range(0,len(b),16):
        b[i:i+16]=b[i:i+16:-1]
    return b

def str2bits(s):
    bits=[]
    for c in s:
        bits+=chr2bits(c)
    return bits

def bytes2bits(s):
    bits=[]
    for b in s:
        bits+=byte2bits(b)
    return bits

def poly2bits(p):
    return [True]+[1==(p>>i)&0x01 for i in range(16)][::-1]

def xorround(bits,poly,i):
    if bits[i]:
        for o in range(17):
            bits[i+o]^=poly[o]

def str2bstr(s):
    return bits2bstr(str2bits(s))

def bits2bstr(bits):
    return ''.join(['1' if i else '0' for i in bits])


def printbits(bits):
    print(''.join(['1' if i else '0' for i in bits]))

def bits2val(bits):
    v=0
    for i in range(16):
        v<<=1
        if bits[i]: v|=0x1
    return v

def takecrc(bits,poly):
    poly=[True]+short2bits(poly)
    bits=bits+[False]*16
    for i in range(len(bits)-16):
        xorround(bits,poly,i)
    return bits[-16:]

def compute_crc(bits):
    return bits2val(takecrc(bits,0x8005)[::-1])
