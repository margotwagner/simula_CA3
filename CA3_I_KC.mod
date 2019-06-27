TITLE ikc.mod       ikc

COMMENT
    Ca activated potassium current I_{K-C}.
ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	    (S) = (millisiemens)
}


? interface
NEURON {
        SUFFIX ikc
        USEION k READ ek WRITE ikc
        USEION caa READ ica
        RANGE gkcbar, gkc, c 
        GLOBAL cinf, ctau
        THREADSAFE : assigned GLOBALs will be per thread
}


PARAMETER {
        gkcbar = 0.015 (S/cm2)      <0, 1e9>
}


STATE {
       c 
}


ASSIGNED {
        v (mV)
        cai (1)
}


BREAKPOINT {
        SOLVE states METHOD cnexp
        gkc = gkcbar * c
        ikc = gkc * xi(cai) * (v - 
}


INITIAL {
}


DERIVATIVE states {
        c' = (cinf - c) / ctau
        cai' = -0.13*ica - 0.075*cai
}


PROCEDURE {
}

FUNCTION xi(y) {
        if (y / 250.0 > 1.0) {
                xi = 1.0
        } else {
                xi = y / 250.0
        }
}


UNITSON
