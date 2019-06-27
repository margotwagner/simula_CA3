TITLE CA3_I_Leak.mod   leak current

COMMENT
    I_Leak
ENDCOMMENT


UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	    (S) = (siemens)
}


? interface
NEURON {
        SUFFIX CA3_I_Leak
        NONSPECIFIC_CURRENT il
        RANGE glbar, gl, el
}

PARAMETER {
        glbar = 0.0001 (S/cm2)
        el = 0 (mV)
}


STATE {
}


ASSIGNED {
        v (mV)
        celsius (degC)
        gl (S/cm2)
        il (mA/cm2)
}


BREAKPOINT {
        SOLVE states METHOD cnexp
        gl = glbar
        il = gl * (v - el) 
}


UNITSON
