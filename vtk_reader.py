import array
import numpy as np


def read_vtk(fp, n1, n2, n3, endian='<', dtype=np.float32):
    """ Scans the VTK data files.

    **Inputs**:

     fp -- Data file pointer\n
     n1 -- No. of points in X1 direction\n
     n2 -- No. of points in X2 direction\n
     n3 -- No. of points in X3 direction\n
     endian -- Endianess of the data\n
     dtype -- datatype

    **Output**:

      Dictionary consisting of variable names as keys and its values.

    """
    ks = []
    vtkvar = []
    while True:
        l = fp.readline()
        try:
            l.split()[0]
        except IndexError:
            pass
        else:
            if l.split()[0] == 'SCALARS':
                ks.append(l.split()[1])
            elif l.split()[0] == 'LOOKUP_TABLE':
                A = array.array(dtype)
                fmt = endian + str(n1 * n2 * n3) + dtype
                nb = np.dtype(fmt).itemsize
                A.fromstring(fp.read(nb))
                vtkvar_buf = np.frombuffer(A, dtype=np.dtype(fmt))
                vtkvar.append(np.reshape(vtkvar_buf, (n1, n2, n3)).transpose())
            else:
                pass
        if l == '':
            break

    vtkvardict = dict(zip(ks, vtkvar))
    return vtkvardict
