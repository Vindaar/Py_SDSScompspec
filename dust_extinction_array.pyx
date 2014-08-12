# dust_extinction_array.pyx

def dust_extinction_array(flux, flux_error, wave, R, A, npix):
   
    cdef double xx1 = 0.0
    cdef double xx2 = 0.0
    cdef double xx3 = 0.0
    cdef double xx4 = 0.0
    cdef double xx5 = 0.0
    cdef double xx6 = 0.0
    cdef double a = 0.0
    cdef double b = 0.0
    cdef double x = 0.0
    cdef double y = 0.0
    cdef double A_lam = 0.0
    cdef double A_V, R_V
    cdef double wavel
    cdef double fluxl
    cdef double flux_errorl
    cdef int    i
    cdef int    np = npix
    A_V = A
    R_V = R
        
    for i in range(np):
        wavel = wave[i]
        fluxl = flux[i]
        flux_errorl = flux_error[i]
        x = 10000.0 / wavel
        if x <= 1.1:
            # NIR 
            xx1 = pow(x, 1.61)
            a = 0.574 * xx1
            b = -0.527 * xx1
        elif x <= 3.3:
            # Optical 
            y = x - 1.82
            xx1 = y * y * y
            xx2 = xx1 * y
            xx3 = xx2 * y
            xx4 = xx3 * y
            xx5 = xx4 * y
            xx6 = xx5 * y
            a = 1.0 + 0.104 * y - 0.609 * y * y + 0.701 * xx1 + 1.137 * xx2 - 1.718 * xx3 - 0.827 * xx4 + 1.647 * xx5 - 0.505 * xx6
            b = 1.952 * y + 2.908 * y * y - 3.989 * xx1 - 7.985 * xx2 + 11.102 * xx3 + 5.491 * xx4 - 10.805 * xx5 + 3.347 * xx6
        elif x <= 8.0:
            # UV
            xx1 = x - 4.67
            xx2 = xx1 * xx1
            xx3 = x - 4.62
            xx4 = xx2 * xx2
            a = 1.752 - 0.316 * x - 0.104 / (xx2 + 0.341)
            b = -3.09 + 1.825 * x + 1.206 / (xx4 + 0.263)
            if x >= 5.9:
                xx3 = x - 5.9
                xx4 = xx3 * xx3
                xx5 = xx4 * xx3
                a += -0.04473 * xx4 - 0.009779 * xx5
                b += 0.213 * xx4 + 0.1207 * xx5
        else:
            # Far UV
            xx1 = x - 8.0
            xx2 = xx1 * xx1
            xx3 = xx2 * xx1
            a = -1.073 - 0.628 * xx1 + 0.137 * xx2 - 0.07 * xx3
            b = 13.67 + 4.257 * xx1 - 0.42 * xx2 + 0.374 * xx3
        # a + b/R_V not (a+b)/R_V!
        A_lam         = a + b/R_V
        A_lam        *= A_V
        redfac        = 10**(0.4*A_lam)
        fluxl        *= redfac
        flux_errorl  *= redfac
        flux[i]       = fluxl
        flux_error[i] = flux_errorl
    return flux, flux_error
