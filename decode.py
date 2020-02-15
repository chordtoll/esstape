#!/usr/bin/env python3
from scipy.io import wavfile
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import struct
import time
import pickle
import esscrc

fs, data = wavfile.read(sys.argv[1])

samples=len(data)

WINDOW=10000

state='IBG'

PROGRESS=True

if not os.path.isfile(os.path.join("blockmaps",sys.argv[1]+'.map')):
    blocks=[]
    print("Locating blocks...")
    i=WINDOW
    sumabs=sum(abs(data[:WINDOW]))  #We find blocks by calculating the average magnitude
    while True:
        if i>=samples:
            break
        if i%10000==0:  #Make a status message every now and then
            if PROGRESS: print('\033[1G    [%s] %10d/%10d (%2d%%): %6d blocks'%(state,i,samples,100*i/samples,len(blocks)),end=' ',flush=True)
        if state=='IBG':
            sumabs+=abs(data[i])
            sumabs-=abs(data[i-WINDOW])
            if sumabs>2000*WINDOW:  #If the average magnitude is above the threshold,
                start=i-WINDOW          #Start of block. Add some extra on the start for good measure.
                state='BLK'
                if PROGRESS: print('\033[1G    [%s] %10d/%10d (%2d%%): %6d blocks'%(state,i,samples,100*i/samples,len(blocks)),end=' ',flush=True)
        else:
            sumabs+=abs(data[i])
            sumabs-=abs(data[i-WINDOW])
            if sumabs<500*WINDOW:  #If the average magnitude is below the threshold,
                end=i                   #End of block.
                blocks.append([start,end])
                print()
                print(start,end,end-start)
                state='IBG'
                if PROGRESS: print('\033[1G    [%s] %10d/%10d (%2d%%): %6d blocks'%(state,i,samples,100*i/samples,len(blocks)),end=' ',flush=True)
        i+=1

    with open(sys.argv[1]+'.map','wb') as f:
        pickle.dump(blocks,f,0)
else:
    with open(sys.argv[1]+'.map','rb') as f:
        blocks=pickle.load(f)


print()           
print()           
print()           

def decode(bits,n): #Converts an array of N bits (LSB first) to an integer
    res=0
    mask=1
    for i in bits:
        if i:
            res|=mask
        mask<<=1
    return res

def tf(b):
    return 'T' if b else 'F' #Yep, it's that simple


DO_GRAPH = False

fails=[]

#Now we decode the blocks
for j,(start,end) in enumerate(blocks):
    print("Decoding block %4d/%4d..."%(j,len(blocks)-1))
    bits=[]
    if DO_GRAPH: end=start+15000
    block=data[start:end]
    dvt=np.convolve(block,[-1,-2,-4,4,2,1],'same') #This is a fancy derivative. Numbers are magic.
    if DO_GRAPH: plt.plot(range(start,end),block)
    prevbit=-1
    samples=len(block)
    blkerror=False
    BIT_THRESH=1000
    zeropol=None
    pl=[]
    plx=[]
    for i in range(20,len(block)-8):
        d1=int(block[i])-int(block[i-8])
        d2=int(block[i])-int(block[i+8])
        pl.append(d1+d2)
        plx.append(start+i)
        #To search for peaks, we check the derivative. If it crosses zero and the sample is larger than nearby samples, we've found a peak.
        if (d1>0)==(d2>0) and abs(d1+d2)>5000 and dvt[i]*dvt[i+1]<=0:   #Check if we're at a peak
            if DO_GRAPH: plt.plot(start+i,block[i],'x')
            if DO_GRAPH: plt.plot(start+i-8,block[i-8],'4')
            if DO_GRAPH: plt.plot(start+i+8,block[i+8],'3')
            if prevbit==-1: #If we're at our first bit,
                zeropol=d1>0  #We know it's a zero. Save the direction for later.
                prevbit=i   #Store the time
                if DO_GRAPH: plt.axvline(start+i,color='blue')   #Graph it! Red line for 0.
                bits.append(False)  #Yep, that's a zero
                break   #We found the zero, time for a new strategy
            '''else:
                tslb=i-prevbit #See how long ago our last bit was
                if tslb<8:  #If it's 0-8, yikes
                    if tslb>3:  #But 0-3 are probably fine
                        print("  BITS CHANGING TOO FAST - NOISY RECORDING?:",tslb)
                        print('  tslb=%d at sample %d'%(tslb,start+i))
                        blkerror=True   #uh oh, sisters
                        break           #i'm outta here
                elif tslb<24: #If it's 8-24, it's our optional pulse
                    pass    #So we can ignore it
                elif tslb<40: #If it's 24-40, we got a bit!
                    bit=(d1>0)!=zeropol   #If its direction is different from our initial zero, it's a one!
                    prevbit=i #Store its time
                    if DO_GRAPH: plt.axvline(start+i,color='green' if bit else 'red')   #Heck yeck we're makin' graphs
                    bits.append(bit) #Put that bit in the bucket
                    if len(bits)>16 and len(bits)<13328:#This is sneaky
                        BIT_THRESH=1000                     #Drop the threshold during the block
                    else:
                        BIT_THRESH=2000                     #Bring it back up at the header/footer
                    if PROGRESS: print('\033[1G          %10d/%10d (%2d%%): %6d bits'%(i,samples,100*i/samples,len(bits)),end=' ',flush=True) #Tell our friend the user that we've got a live one
                    
                else: #More than 40, things are wrong
                    print()
                    print("  BITS CHANGING TOO SLOW - PERHAPS TAPE IS REVERSED?")
                    print("  tslb=%d at sample %d"%(tslb,start+i)) #Fun fact, a TSLB of 48 means the tape is probably backwards
                    blkerror=True
                    break'''
    #if zeropol is None:
    #    print(5/0)
    plt.plot(plx,pl)
    while True:
        if (prevbit+48>(end-start)): break
        dval=[]
        yval=[]
        for i in range(24,40):
            d1=int(block[prevbit+i])-int(block[prevbit+i-8])
            d2=int(block[prevbit+i])-int(block[prevbit+i+8])
            ds=abs(d1+d2)
            dval.append(ds)
            yval.append(block[prevbit+i])
        mpos=prevbit+24+dval.index(max(dval))
        #if DO_GRAPH: plt.plot(start+mpos,block[mpos],'x')
        #if DO_GRAPH: plt.plot(start+prevbit+24,block[prevbit+24],'4')
        #if DO_GRAPH: plt.plot(start+prevbit+40,block[prevbit+40],'3')
        
        bit=((block[mpos]-block[mpos-8])>0)!=zeropol

        if PROGRESS: print('\033[1G          %10d/%10d (%2d%%): %6d bits'%(mpos,samples,100*mpos/samples,len(bits)),end=' ',flush=True) #Tell our friend the user that we've got a live one

        bits.append(bit)
        if (len(bits)==13344): break
        
        if DO_GRAPH: plt.axvline(start+mpos,color='green' if bit else 'red')
        
        #if DO_GRAPH: plt.plot(range(start+prevbit+24,start+prevbit+40),dval)
        prevbit=mpos
    if DO_GRAPH: plt.show()
    
    if blkerror:    #include <yikes.h>
        print('FAILED TO DECODE. SKIPPING') #If the block's no good, throw it out
        fails.append((j,start+i))
        continue

    print()
    print(len(bits))

    if len(bits)==13345: #If we catch a ghost bit,
        bits=bits[:-1]     #Drop it. This is risky, but CRC

    if (len(bits)%16)!=0: #We need an integer number of words
        print('INCORRECT LENGTH. SKIPPING')
        fails.append((j,len(bits)))
        continue

    preamble=decode(bits[:16],16) #Get those tasty fields out of the header
    block=decode(bits[16:27],11)
    track=decode(bits[27:29],2)
    eof=bits[29]
    deof=bits[30]
    crc=decode(bits[13312:13328],16)
    ccrc=esscrc.compute_crc(bits[16:-32])

    print('    Header:') #Tell the user so they can celebrate!
    print('        Preamble: %04x'%preamble)
    print('        Block   : %03x'%block)
    print('        Track   : %1x'%track)
    print('        EOF     : %s'%tf(eof))
    print('        DEOF    : %s'%tf(deof))
    print('        CRC     : %04x'%crc)
    print('    Computed CRC: %04x'%ccrc)
    print()
    print()
    print()

    if crc!=ccrc:
        print("CRC MISMATCH. SKIPPING")
        #input("Press enter to continue...")
        continue

    if not os.path.isdir(str(track)): #Make a folder for the track
        os.mkdir(str(track))
    with open(os.path.join(str(track),'%04d.bin'%block),'wb') as f: #Write the block
        for i in range(0,len(bits),16):
            f.write(struct.pack('>H',decode(bits[i:i+16],16)))

if fails:
    print("Failed to decode blocks:")
    for n,offs in fails:
        print("\t%d @ %d"%(n,offs))
