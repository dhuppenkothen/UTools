import numpy as numpy

def resample(n, step_size):
    output = []
    for i in numpy.arange(0, len(n), step_size):
        total = 0
        print "Bin is " + str(i)

        prev_frac = int(i+1) - i
        prev_bin = int(i)
        print "Fractional part of bin %d is %f"  %(prev_bin, prev_frac)
        total += prev_frac * n[prev_bin]

        if i + step_size < len(n):
            # Fractional part of next bin:
            next_frac = i+step_size - int(i+step_size)
            next_bin = int(i+step_size)
            print "Fractional part of bin %d is %f"  %(next_bin, next_frac)
            total += next_frac * n[next_bin]

        print "Fully included bins: %d to %d" % (int(i+1), int(i+step_size)-1)
        total += sum(n[int(i+1):int(i+step_size)])
        output.append(total)
    return output


if __name__ == "__main__":
    my_array = numpy.array([19,11,13,11,12,10])
    print regrid(my_array, 1.2)
