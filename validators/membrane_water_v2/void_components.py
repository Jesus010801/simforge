"""Linear-time six-connected labeling with periodic seam equivalence."""
import numpy as np
from scipy import ndimage


def components(mask, periodic=(True, True, True)):
    labels, count = ndimage.label(mask)
    parent = np.arange(count + 1)

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for axis, wraps in enumerate(periodic):
        if not wraps:
            continue
        a, b = np.take(labels, 0, axis), np.take(labels, -1, axis)
        pairs = np.unique(np.stack((a.ravel(), b.ravel()), axis=1), axis=0)
        for u, v in pairs:
            if u and v:
                ru, rv = root(u), root(v)
                parent[max(ru, rv)] = min(ru, rv)
    roots = np.array([root(i) for i in range(count + 1)])
    unique, compact = np.unique(roots, return_inverse=True)
    return compact[labels].astype(np.int32), len(unique) - 1


def neighbors(flat, shape):
    x, rem = divmod(flat, shape[1] * shape[2])
    y, z = divmod(rem, shape[2])
    for axis, coord in enumerate((x, y, z)):
        stride = (shape[1] * shape[2], shape[2], 1)[axis]
        for step in (-1, 1):
            yield flat + (((coord + step) % shape[axis]) - coord) * stride
