def calculate_fcs(payload):
    """Calculate the LRC of the string.

    The Frame Check Sequence is calculated as follows: For each
    byte of the payload, interpret this as an ascii character. For
    each ascii character, get the two-character hex
    representation. the FCS is the XOR of the leftmost of these
    characters, followed by the XOR of the rightmost of these
    characters.

    """
    left = 0
    right = 0
    for byte in payload:
        hexcode = hex(ord(byte))
        print("Byte = [{}]  Hexcode = [{}]".format(byte, hexcode))
        print("left = [{}]([{}]) right = [{}]([{}])".format(hexcode[2], bin(int(hexcode[2],16)), hexcode[3], bin(int(hexcode[3],16))))
        left ^= int(hexcode[2],16)
        right ^= int(hexcode[3],16)

    print("Left = [{}]([{}])  Right = [{}]([{}])".format(left, bin(left), right, bin(right)))
    result = "{:X}{:X}".format(left, right)
    print("Result = {}".format(result))

if __name__ == '__main__':
#    payload = "@00MS"
    payload = "@00RD01100012"
    calculate_fcs(payload)
              
