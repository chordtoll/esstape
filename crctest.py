def crc16(data : bytearray, offset , length,poly):
    if data is None or offset < 0 or offset > len(data)- 1 and offset+length > len(data):
        return 0
    crc = 0x0000
    for i in range(0, length):
        crc ^= data[offset + i] << 8
        for j in range(0,8):
            if (crc & 0x8000) > 0:
                crc =(crc << 1) ^ poly 
            else:
                crc = crc << 1
    return crc & 0xFFFF


