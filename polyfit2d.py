#!/usr/bin/env python

import unittest
from numpy.polynomial import polynomial
import numpy as np


def polyfit2d(x, y, f, deg):
    x = np.asarray(x)
    y = np.asarray(y)
    f = np.asarray(f)
    deg = np.asarray(deg)
    vander = polynomial.polyvander2d(x, y, deg)
    vander = vander.reshape((-1,vander.shape[-1]))
    f = f.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, f)[0]
    return c.reshape(deg+1)

class MyTest(unittest.TestCase):

    def setUp(self):
        return self

    def test_1(self):
        x = [-1,2,3]
        y = [ 4,5,6]
        z = [[1,2,3],[4,5,6],[7,8,9]]
        print('x.shape = ',np.asarray(x).shape,', y.shape = ',np.asarray(y).shape,', z.shape = ',np.asarray(z).shape)
        self._test_fit(
            x,
            y,
            z,
            [2,2])

    def test_2(self):
        self._test_fit(
            [-1,2],
            [ 4,5],
            [[1,2],[4,5]],
            [1,1])

    def test_3(self):
        x = [-1,2,3]
        y = [ 4,5]
        z = [[1,2],[4,5],[7,8]]
        print('x.shape = ',np.asarray(x).shape,', y.shape = ',np.asarray(y).shape,', z.shape = ',np.asarray(z).shape)
        self._test_fit(
            x,
            y,
            z,
            [2,1])

    def test_4(self):
        self._test_fit(
            [-1,2,3],
            [ 4,5],
            [[1,2],[4,5],[0,0]],
            [2,1])

    def test_5(self):
        self._test_fit(
            [-1,2,3],
            [ 4,5],
            [[1,2],[4,5],[0,0]],
            [1,1])

    def test_6(self):
        x = np.ndarray(shape=(1700), dtype=np.int16)
        y = np.ndarray(shape=(2100), dtype=np.int16)
        z = np.ones(shape=(2100,1700), dtype=np.int32)
        self._test_fit(
            x,
            y,
            z,
            [3,3])

    def _test_fit(self, x, y, c, deg):
        from numpy.polynomial import polynomial
        import numpy as np
        X = np.array(np.meshgrid(x,y))
        f = polynomial.polyval2d(X[0], X[1], c)
        c1 = polyfit2d(X[0], X[1], f, deg)
        print('X.shape = ',X.shape,', f.shape = ',f.shape,', c1.shape = ',c1.shape)
        np.testing.assert_allclose(c1,
                                np.asarray(c)[:deg[0]+1,:deg[1]+1],
                                atol=1e-12)

unittest.main()
