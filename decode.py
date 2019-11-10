#!/usr/bin/env python3
from scipy.io import wavfile
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import struct
import time
fs, data = wavfile.read(sys.argv[1])
samples=len(data)


WINDOW=10000


state='IBG'

blocks=[]
print("Locating blocks...")
i=WINDOW
sumabs=sum(abs(data[:WINDOW]))
while True:
    if i>=samples:
        break
    if i%10000==0:
        print('\033[1G    [%s] %10d/%10d (%2d%%): %6d blocks'%(state,i,samples,100*i/samples,len(blocks)),end=' ',flush=True)
    if state=='IBG':
        sumabs+=abs(data[i])
        sumabs-=abs(data[i-WINDOW])
        if sumabs>2500*WINDOW:
            start=i-WINDOW
            state='BLK'
            print('\033[1G    [%s] %10d/%10d (%2d%%): %6d blocks'%(state,i,samples,100*i/samples,len(blocks)),end=' ',flush=True)
    else:
        sumabs+=abs(data[i])
        sumabs-=abs(data[i-WINDOW])
        if sumabs<2500*WINDOW:
            end=i
            blocks.append([start,end])
            state='IBG'
            print('\033[1G    [%s] %10d/%10d (%2d%%): %6d blocks'%(state,i,samples,100*i/samples,len(blocks)),end=' ',flush=True)
    i+=1
print()           
print()           
print()           

def decode(bits,n):
    res=0
    mask=1
    for i in bits:
        if i:
            res|=mask
        mask<<=1
    return res

def tf(b):
    return 'T' if b else 'F'

for j,(start,end) in enumerate(blocks):
    print("Decoding block %4d/%4d..."%(j,len(blocks)-1))
    bits=[]
    block=data[start:end]
    dvt=np.convolve(block,[-1,-2,-4,4,2,1],'same')
    #plt.plot(range(start,end),block)
    #plt.plot(range(start,end),dvt)
    prevbit=-1
    samples=len(block)
    blkerror=False
    for i in range(len(block)-1):
        if abs(block[i])>5000 and dvt[i]*dvt[i+1]<=0:
            if prevbit==-1:
                #print("Found first bit")
                zeropol=block[i]>0
                prevbit=i
                #plt.axvline(start+i,color='red')
                bits.append(False)
            else:
                tslb=i-prevbit
                if tslb<8:
                    if tslb>3:
                        print("  BITS CHANGING TOO FAST - NOISY RECORDING?:",tslb)
                        print('  tslb=%d at sample %d'%(tslb,start+i))
                        blkerror=True
                        break
                elif tslb<24:
                    pass
                    #print("  Inter-bit transition")
                elif tslb<40:
                    bit=(block[i]>0)!=zeropol
                    #print("%d Bit goes here"%(1 if bit else 0))
                    prevbit=i
                    #plt.axvline(start+i,color='gray')#('green' if bit else 'red'))
                    bits.append(bit)
                    print('\033[1G          %10d/%10d (%2d%%): %6d bits'%(i,samples,100*i/samples,len(bits)),end=' ',flush=True)
                    
                else:
                    print()
                    print("  BITS CHANGING TOO SLOW - PERHAPS TAPE IS REVERSED?")
                    print("  tslb=%d at sample %d"%(tslb,start+i))
                    blkerror=True
                    break

    plt.show()
    if blkerror or (len(bits)%16)!=0:
        print('FAILED TO DECODE. SKIPPING')
        continue
    print()
    if (len(bits)%16)!=0:
        print('INCORRECT LENGTH. SKIPPING')
        continue

    preamble=decode(bits[:16],16)
    block=decode(bits[16:27],11)
    track=decode(bits[27:29],2)
    eof=bits[29]
    deof=bits[30]

    print('    Header:')
    print('        Preamble: %04x'%preamble)
    print('        Block   : %03x'%block)
    print('        Track   : %1x'%track)
    print('        EOF     : %s'%tf(eof))
    print('        DEOF    : %s'%tf(deof))
    print()
    print()
    print()

    if not os.path.isdir(str(track)):
        os.mkdir(str(track))
    with open(os.path.join(str(track),'%04d.bin'%block),'wb') as f:
        for i in range(0,len(bits),16):
            f.write(struct.pack('>H',decode(bits[i:i+16],16)))
