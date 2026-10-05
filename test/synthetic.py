import numpy as np

def chessboard(gridn, square, shape=None, margin=80):
    r'''Returns (image, corners), with corners in mrgingham's order'''
    n = gridn + 1
    if shape is None:
        shape = (n*square + 2*margin,) * 2
    image = np.full(shape, 255.)
    for i in range(n):
        for j in range(n):
            if (i + j) % 2 == 0:
                image[margin + i*square: margin + (i+1)*square,
                      margin + j*square: margin + (j+1)*square] = 0

    # Sharp boards are only found at some square sizes
    k = np.array((1., 4., 6., 4., 1.)) / 16.
    for axis in (0, 1):
        image = np.apply_along_axis(lambda r: np.convolve(r, k, 'same'), axis, image)
    image = image.astype(np.uint8)

    # corners are on pixel boundaries; pixel centers are at integers
    c = margin + square*np.arange(1, n) - 0.5
    x, y = np.meshgrid(c, c)
    return image, np.column_stack((x.ravel(), y.ravel()))
