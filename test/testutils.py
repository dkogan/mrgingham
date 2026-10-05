# Copied from mrcal, git revision 214d13a1ab36, test/testutils.py. Removed
# confirm_covariances_equal(), the only thing that needed mrcal

import sys
import numpy as np
import numpysane as nps
import os
import re
from inspect import currentframe

Nchecks = 0
NchecksFailed = 0

# no line breaks. Useful for test reporting. Yes, this sets global state, but
# we're running a test harness. This is fine
np.set_printoptions(linewidth=1e10, suppress=True)



def percentile_compat(*args, **kwargs):
    r'''Wrapper for np.percentile() to handle their API change

In numpy 1.24 the "interpolation" kwarg was renamed to "method". I need to pass
the right thing to work with both old and new numpy. This function tries the
newer method, and if that fails, uses the old one. The test is only done the
first time.

It is assumed that this is called with the old 'interpolation' key.

    '''

    if not 'interpolation' in kwargs or \
       percentile_compat.which == 'interpolation':
        return np.percentile(*args, **kwargs)

    kwargs_no_interpolation = dict(kwargs)
    del kwargs_no_interpolation['interpolation']

    if percentile_compat.which == 'method':
        return np.percentile(*args, **kwargs_no_interpolation,
                             method = kwargs['interpolation'])

    # Need to detect

    try:
        result = np.percentile(*args, **kwargs_no_interpolation,
                               method = kwargs['interpolation'])
        percentile_compat.which = 'method'
        return result
    except:
        percentile_compat.which = 'interpolation'
        return np.percentile(*args, **kwargs)

percentile_compat.which = None


def test_location():
    r'''Reports string describing current location in the test'''


    filename_this = os.path.split( __file__ )[1]

    frame = currentframe().f_back.f_back

    # I keep popping the stack until I leave the testutils file and I'm not in a
    # function called "check"
    while frame:
        if frame.f_back is None or \
           (not frame.f_code.co_filename.endswith(filename_this) and
            frame.f_code.co_name != "check" ):
            break
        frame = frame.f_back

    testfile = os.path.split(frame.f_code.co_filename)[1]
    try:
        return "{}:{}".format(testfile, frame.f_lineno)
    except:
        return ''


def print_red(x):
    """print the message in red"""
    sys.stdout.write("\x1b[31m" + test_location() + ": " + x + "\x1b[0m\n")


def print_green(x):
    """Print the message in green"""
    sys.stdout.write("\x1b[32m" + test_location() + ": " + x + "\x1b[0m\n")

def print_blue(x):
    """Print the message in blue"""
    sys.stdout.write("\x1b[34m" + test_location() + ": " + x + "\x1b[0m\n")



def relative_scale(a,b,
                   *,
                   smooth_radius = None,
                   eps           = 1e-6):
    if smooth_radius is not None and smooth_radius > 0:
        d = smooth_radius*2 + 1
        f = np.ones((d,),) / d
        a = np.convolve(a, f, mode='same')
        b = np.convolve(b, f, mode='same')
    return (np.abs(a) + np.abs(b)) / 2 + eps

def relative_diff(a,b,
                  *,
                  smooth_radius = None,
                  eps           = 1e-6):
    return (a - b) / relative_scale(a,b,
                                    eps           = eps,
                                    smooth_radius = smooth_radius)

def confirm_equal(x, xref,
                  *,
                  msg='',
                  eps=1e-6,
                  reldiff_eps = 1e-6,
                  reldiff_smooth_radius = None,
                  relative=False,
                  worstcase=False,
                  percentile=None,
                  r=False):
    r'''If x is equal to xref, report test success.

    msg identifies this check. eps sets the RMS equality tolerance. The x,xref
    arguments can be given as many different types. This function tries to do
    the right thing.

    if relative: I look at a relative error:
                 err = (a-b) / ((abs(a)+abs(b))/2 + eps)
                 a,b can be smoothed with a kernel of the given smooth_radius
    else:        I look at absolute error:
                 err = a-b

    if worstcase: I look at the worst-case error
                  error = np.max(np.abs(err))
    elif percentile is not None: I look at the given point in the error distribution
                  error = percentile_compat(np.abs(err), percentile)
    else:         RMS error
                  error = np.sqrt(nps.norm2(err) / len(err))

    if r: we are comparing rodrigues rotations. More than one set of r values
          can represent the same rotation

          Let k be an integer. r = th * vaxis.
          Changing th -> th + k*2pi implies the same rotation
          Changing vaxis -> -vaxis and th -> 2pi-th also implies the same rotation
          I normalize the inputs first by finding the rotation with the smallest th

    '''

    if r:
        if not (x.shape[-1] == 3 and xref.shape[-1] == 3):
            raise Exception("confirm_equal(r=True) only makes sense if x and xref have shape (...,3)")

        def normalize_r(r):
            th = nps.mag(r)
            v = r / nps.dummy(th, -1)
            th %= 2.*np.pi
            # th is now in [0,2pi)
            if th > np.pi:
                th = 2.*np.pi - th
                v  *= -1
            # th is in [0,pi)
            return th * v

        x    = normalize_r(x)
        xref = normalize_r(xref)






    global Nchecks
    global NchecksFailed
    Nchecks = Nchecks + 1

    # strip all trailing whitespace in each line, in case these are strings
    if isinstance(x, str):
        x = re.sub('[ \t]+(\n|$)', '\\1', x)
    if isinstance(xref, str):
        xref = re.sub('[ \t]+(\n|$)', '\\1', xref)

    # convert data to numpy if possible
    try:
        xref = np.array(xref)
    except:
        pass
    try:
        x = np.array(x)
    except:
        pass

    try:  # flatten array if possible
        x = x.ravel()
        xref = xref.ravel()
    except:
        pass

    try:
        N = x.shape[0]
    except:
        N = 1
    try:
        Nref = xref.shape[0]
    except:
        Nref = 1

    if N != Nref:

        # Comparing an array to a scalar reference is allowed
        if Nref == 1:
            xref = np.ones((N,), dtype=float) * xref
            Nref = N
        else:
            print_red(("FAILED{}: mismatched array sizes: N = {} but Nref = {}. Arrays: \n" +
                       "x = {}\n" +
                       "xref = {}").
                      format((': ' + msg) if msg else '',
                             N, Nref,
                             x, xref))
            NchecksFailed = NchecksFailed + 1
            return False

    if N != 0:
        try:  # I I can subtract, get the error that way
            if relative:
                diff = relative_diff(x, xref,
                                     eps           = reldiff_eps,
                                     smooth_radius = reldiff_smooth_radius)
            else:
                diff = x - xref

            if worstcase:
                what = 'worst-case'
                err  = np.max(np.abs(diff))
            elif percentile is not None:
                what = f'{percentile}%-percentile'
                err  = percentile_compat(np.abs(diff), percentile, interpolation='higher')
            else:
                what = 'RMS'
                err  = np.sqrt(nps.norm2(diff) / len(diff))

            if not np.all(np.isfinite(err)):
                print_red(f"FAILED{(': ' + msg) if msg else ''}: Some comparison results are NaN or Inf. {what}. error_x_xref =\n{nps.cat(err,x,xref)}")
                NchecksFailed = NchecksFailed + 1
                return False
            if err > eps:
                print_red(f"FAILED{(': ' + msg) if msg else ''}: {what} error = {err}. x_xref_err =\n{nps.cat(x,xref,diff)}")
                NchecksFailed = NchecksFailed + 1
                return False
        except:  # Can't subtract. Do == instead
            if not np.array_equal(x, xref):
                print_red(f"FAILED{(': ' + msg) if msg else ''}: x_xref =\n{nps.cat(x,xref)}")
                NchecksFailed = NchecksFailed + 1
                return False
    print_green("OK" + (': ' + msg) if msg else '')
    return True


def confirm(x, msg=''):
    r'''If x is true, report test success.

    msg identifies this check'''

    global Nchecks
    global NchecksFailed
    Nchecks = Nchecks + 1

    if not x:
        print_red("FAILED{}".format((': ' + msg) if msg else ''))
        NchecksFailed = NchecksFailed + 1
        return False
    print_green("OK{}".format((': ' + msg) if msg else ''))
    return True


def confirm_raises(f, msg=''):
    r'''If f() raises an exception, report test success.

    msg identifies this check'''

    global Nchecks
    global NchecksFailed
    Nchecks = Nchecks + 1

    try:
        f()
        print_red("FAILED{}".format((': ' + msg) if msg else ''))
        NchecksFailed = NchecksFailed + 1
        return False
    except:
        print_green("OK{}".format((': ' + msg) if msg else ''))
        return True


def confirm_does_not_raise(f, msg=''):
    r'''If f() raises an exception, report test failure.

    msg identifies this check'''

    global Nchecks
    global NchecksFailed
    Nchecks = Nchecks + 1

    try:
        f()
        print_green("OK{}".format((': ' + msg) if msg else ''))
        return True
    except:
        print_red("FAILED{}".format((': ' + msg) if msg else ''))
        NchecksFailed = NchecksFailed + 1
        return False


def finish():
    r'''Finalize the executed tests.

    Prints the test summary. Exits successfully iff all the tests passed.

    '''
    if not Nchecks and not NchecksFailed:
        print_red("No tests defined")
        sys.exit(0)

    if NchecksFailed:
        print_red("Some tests failed: {} out of {}".format(NchecksFailed, Nchecks))
        sys.exit(1)

    print_green("All tests passed: {} total".format(Nchecks))
    sys.exit(0)
