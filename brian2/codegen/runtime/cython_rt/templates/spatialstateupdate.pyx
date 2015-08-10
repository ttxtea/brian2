{# USES_VARIABLES { N, _invr, Ri, Cm, dt, area, diameter, length,
                        ab_star0, ab_star1, ab_star2,
                        ab_plus0, ab_plus1, ab_plus2,
                        ab_minus0, ab_minus1, ab_minus2,
                        _starts, _ends,
                        Ri, diameter, length, area,
                        _invr, _invr0, _invrn, b_plus, b_minus } #}

{% extends 'common.pyx' %}

{% macro before_run() %}
{{ preamble() }}
def main(_namespace):
    cdef int _idx
    cdef int _vectorisation_idx
    {{ load_namespace | autoindent }}
    cdef int _i
    cdef int _counter
    cdef int _first
    cdef int _last
    cdef double __invr0
    cdef double __invrn

    cdef double _Ri = {{Ri}}[0]  # Ri is a shared variable

    {% if owner.morphology.type == 'soma' %}
    # Correction for soma (a bit of a hack),
    #  so that it has negligible axial resistance
    {{length}}[0] = {{diameter}}[0] * 0.01;
    {% endif %}

    # Inverse axial resistance
    for _i in range(1, N):
        {{_invr}}[_i] = (M_PI / (2 * _Ri) * ({{diameter}}[_i-1] * {{diameter}}[_i]) /
                       ({{length}}[_i-1] + {{length}}[_i]))
    # Note: this would give nan for the soma
    # Cut branches
    for _i in range(_num{{_starts}}):
        {{_invr}}[{{_starts}}[_i]] = 0

    # Linear systems
    # The particular solution
    # a[i,j]=ab[u+i-j,j]   --  u is the number of upper diagonals = 1
    for _i in range(N):
        {{ab_star1}}[_i] = (-({{Cm}}[_i] / dt) - {{_invr}}[_i] / {{area}}[_i])
    for _i in range(1, N):
        {{ab_star0}}[_i] = {{_invr}}[_i] / {{area}}[_i-1]
        {{ab_star2}}[_i-1] = {{_invr}}[_i] / {{area}}[_i]
        {{ab_star1}}[_i-1] -= {{_invr}}[_i] / {{area}}[_i-1]
    for _i in range(N):
        # Homogeneous solutions
        {{ab_plus0}}[_i] = {{ab_star0}}[_i]
        {{ab_minus0}}[_i] = {{ab_star0}}[_i]
        {{ab_plus1}}[_i] = {{ab_star1}}[_i]
        {{ab_minus1}}[_i] = {{ab_star1}}[_i]
        {{ab_plus2}}[_i] = {{ab_star2}}[_i]
        {{ab_minus2}}[_i] = {{ab_star2}}[_i]

    # Set the boundary conditions
    for _counter in range(_num{{_starts}}):
        _first = {{_starts}}[_counter]
        _last = {{_ends}}[_counter]
        # Inverse axial resistances at the ends: r0 and rn
        __invr0 = (M_PI / (2 * _Ri) * {{diameter}}[_first] * {{diameter}}[_first]  /
                              {{length}}[_first])
        __invrn = (M_PI / (2 * _Ri) * {{diameter}}[_last] * {{diameter}}[_last] /
                              {{length}}[_last])
        {{_invr0}}[_counter] = __invr0
        {{_invrn}}[_counter] = __invrn
        # Correction for boundary conditions
        {{ab_star1}}[_first] -= (__invr0 / {{area}}[_first])
        {{ab_star1}}[_last] -= (__invrn / {{area}}[_last])
        {{ab_plus1}}[_first] -= (__invr0 / {{area}}[_first])
        {{ab_plus1}}[_last] -= (__invrn / {{area}}[_last])
        {{ab_minus1}}[_first] -= (__invr0 / {{area}}[_first])
        {{ab_minus1}}[_last] -= (__invrn / {{area}}[_last])
        # RHS for homogeneous solutions
        {{b_plus}}[_last] = -(__invrn / {{area}}[_last])
        {{b_minus}}[_first] = -(__invr0 / {{area}}[_first])
{% endmacro %}

{% block maincode %}
    {# USES_VARIABLES { Cm, dt, v, N,
                        ab_star0, ab_star1, ab_star2, b_plus,
                        ab_plus0, ab_plus1, ab_plus2, b_minus,
                        ab_minus0, ab_minus1, ab_minus2,
                        v_star, u_plus, u_minus, gtot_all
                        _P, _B, _morph_i, _morph_parent_i, _starts, _ends,
                        _invr0, _invrn} #}

    cdef double[:] c = _numpy.zeros(N, dtype=_numpy.double)

    cdef double ai
    cdef double bi
    cdef double _m

    cdef int _i
    cdef int _j

    cdef int _n_segments = _num{{_B}}
    cdef int __morph_i
    cdef int __i_parent
    cdef int __first
    cdef int __last
    cdef double __invr0
    cdef double __invrn
    cdef int col
    cdef int i_pivot
    cdef double pivot_magnitude
    cdef double pivot_element
    cdef double pivot_factor
    cdef double tmp

    _vectorisation_idx = 1

    ## MAIN CODE ######
    {{scalar_code|autoindent}}

    # Tridiagonal solving
    # Pass 1
    for _i in range(N):
        _idx = _i
        _vectorisation_idx = _idx

        {{vector_code|autoindent}}
        {{gtot_all}}[_idx] = _gtot

        {{v_star}}[_i] = -({{Cm}}[_i]/dt*{{v}}[_i])-_I0 # RHS -> v_star (solution)
        bi={{ab_star1}}[_i] - {{gtot_all}}[_i] # main diagonal
        if (_i<N-1):
            c[_i]= {{ab_star0}}[_i+1] # superdiagonal
        if (_i>0):
            ai = {{ab_star2}}[_i-1] # subdiagonal
            _m = 1.0/(bi-ai*c[_i-1])
            c[_i] = c[_i]*_m
            {{v_star}}[_i] = ({{v_star}}[_i] - ai*{{v_star}}[_i-1])*_m
        else:
            c[0]=c[0]/bi
            {{v_star}}[0]={{v_star}}[0]/bi

    for _i in range(N-2, -1, -1):
        {{v_star}}[_i]={{v_star}}[_i] - c[_i]*{{v_star}}[_i+1]
    
    # Pass 2
    for _i in range(N):
        {{u_plus}}[_i] = {{b_plus}}[_i] # RHS -> v_star (solution)
        bi = {{ab_plus1}}[_i]-{{gtot_all}}[_i] # main diagonal
        if (_i<N-1):
            c[_i] = {{ab_plus0}}[_i+1] # superdiagonal
        if (_i>0):
            ai = {{ab_plus2}}[_i-1] # subdiagonal
            _m = 1.0/(bi-ai*c[_i-1])
            c[_i] = c[_i]*_m
            {{u_plus}}[_i] = ({{u_plus}}[_i] - ai*{{u_plus}}[_i-1])*_m
        else:
            c[0]=c[0]/bi
            {{u_plus}}[0] = {{u_plus}}[0]/bi
    for _i in range(N-2, -1, -1):
        {{u_plus}}[_i] = {{u_plus}}[_i] - c[_i]*{{u_plus}}[_i+1]
    
    # Pass 3
    for _i in range(N):
        {{u_minus}}[_i] = {{b_minus}}[_i] # RHS -> v_star (solution)
        bi = {{ab_minus1}}[_i] - {{gtot_all}}[_i] # main diagonal
        if (_i<N-1):
            c[_i] = {{ab_minus0}}[_i+1] # superdiagonal
        if (_i>0):
            ai = {{ab_minus2}}[_i-1] # subdiagonal
            _m = 1.0/(bi-ai*c[_i-1])
            c[_i] = c[_i]*_m
            {{u_minus}}[_i] = ({{u_minus}}[_i] - ai*{{u_minus}}[_i-1])*_m
        else:
            c[0] = c[0]/bi
            {{u_minus}}[0] = {{u_minus}}[0]/bi
    for _i in range(N-2, -1, -1):
        {{u_minus}}[_i] = {{u_minus}}[_i] - c[_i]*{{u_minus}}[_i+1]

    # Prepare matrix for solving the linear system
    for _i in range(_num{{_B}}):
        {{_B}}[_i] = 0.
    for _i in range(_num{{_P}}):
        {{_P}}[_i] = 0.
    for _i in range(_n_segments - 1):
        __morph_i = {{_morph_i}}[_i]
        __i_parent = {{_morph_parent_i}}[_i]
        __first = {{_starts}}[_i]
        __last = {{_ends}}[_i]
        __invr0 = {{_invr0}}[_i]
        __invrn = {{_invrn}}[_i]
        # Towards parent
        if __morph_i == 1: # first branch, sealed end
            {{_P}}[0] = {{u_minus}}[__first] - 1
            {{_P}}[0 + 1] = {{u_plus}}[__first]
            {{_B}}[0] = -{{v_star}}[__first]
        else:
            {{_P}}[__i_parent*_n_segments + __i_parent] += (1 - {{u_minus}}[__first]) * __invr0
            {{_P}}[__i_parent*_n_segments + __morph_i] -= {{u_plus}}[__first] * __invr0
            {{_B}}[__i_parent] += {{v_star}}[__first] * __invr0
        # Towards children
        {{_P}}[__morph_i*_n_segments + __morph_i] = (1 - {{u_plus}}[__last]) * __invrn
        {{_P}}[__morph_i*_n_segments + __i_parent] = -{{u_minus}}[__last] * __invrn
        {{_B}}[__morph_i] = {{v_star}}[__last] * __invrn

    # Solve the linear system (the result will be in _B in the end)
    for _i in range(_n_segments):
        # find pivot element
        i_pivot = _i
        pivot_magnitude = fabs({{_P}}[_i*_n_segments + _i])
        for _j in range(_i+1, _n_segments):
            if fabs({{_P}}[_j*_n_segments + _i]) > pivot_magnitude:
                i_pivot = _j
                pivot_magnitude = fabs({{_P}}[_j*_n_segments + _i])

        if pivot_magnitude == 0.:
            raise ValueError('Matrix is singular')

        # swap rows
        if _i != i_pivot:
            for col in range(_i, _n_segments):
                tmp = {{_P}}[_i*_n_segments + col]
                {{_P}}[_i*_n_segments + col] = {{_P}}[i_pivot*_n_segments + col]
                {{_P}}[i_pivot*_n_segments + col] = tmp

        # Deal with rows below
        for _j in range(_i+1, _n_segments):
            pivot_element = {{_P}}[_j*_n_segments + _i]
            if pivot_element == 0.0:
                continue
            pivot_factor = pivot_element/{{_P}}[_i*_n_segments + _i]
            for _k in range(_i+1, _n_segments):
                {{_P}}[_j*_n_segments + _k] -= {{_P}}[_i*_n_segments + _k]*pivot_factor
            {{_B}}[_j] -= {{_B}}[_i]*pivot_factor
            {{_P}}[_j*_n_segments + _i] = 0

    # Back substitution
    for _i in range(_n_segments-1, -1, -1):
        # substitute all the known values
        for _j in range(_n_segments-1, _i, -1):
            {{_B}}[_i] -= {{_P}}[_i*_n_segments + _j]*{{_B}}[_j]
            {{_P}}[_i*_n_segments + _j] = 0
        # divide by the diagonal element
        {{_B}}[_i] /= {{_P}}[_i*_n_segments + _i]
        {{_P}}[_i*_n_segments + _i] = 1

    # Linear combination
    for _i in range(_n_segments-1):
        __morph_i = {{_morph_i}}[_i]
        __i_parent = {{_morph_parent_i}}[_i]
        __first = {{_starts}}[_i]
        __last = {{_ends}}[_i]
        for _j in range(__first, __last+1):
            if (_j < _num{{v}}):
                {{v}}[_j] = {{v_star}}[_j] + {{_B}}[__i_parent] * {{u_minus}}[_j] + {{_B}}[__morph_i] * {{u_plus}}[_j]

{% endblock %}
