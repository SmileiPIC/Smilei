import os, sys

'''
Makes a hexdump from the command line, in Smilei's format.
This code is intended as replacement for xxd command

Usage:
    python hexdump.py source_file dest_file
'''

if __name__ == "__main__":
    # Read input
    basename = os.path.splitext(os.path.basename(sys.argv[1]))[0] + "_py"
    with open(sys.argv[1], 'r') as src:
        src_text = src.read()

    # hexdump, inspired from https://gist.github.com/sbz/1080258
    length = 12
    totlen = len(src_text)
    try: # python 2
        lines = [', '.join(["0x%02x" % ord(x) for x in src_text[c:c+length]]) for c in xrange(0, totlen, length)]
    except: # python 3
        lines = [', '.join(["0x%02x" % ord(x) for x in src_text[c:c+length]]) for c in range(0, totlen, length)]

    # Write out
    with open(sys.argv[2], 'w') as dest:
        dest.write("unsigned char "+basename+"[] = {\n  ")
        dest.write(',\n  '.join(lines))
        dest.write("\n};\n")
        dest.write("unsigned int "+basename+"_len = "+str(totlen)+";\n")
