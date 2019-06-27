TITLE CA3_I_KC.mod  potassium

COMMENT
    Ca activated potassium current I_{K-C}.
ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	    (S) = (siemens)
}


? interface
NEURON {
        SUFFIX CA3_I_KC
        USEION k READ ek WRITE ik
        RANGE gkcbar, gkc, c 
        GLOBAL cinf, ctau
        THREADSAFE : assigned GLOBALs will be per thread
}


PARAMETER {
        gkcbar = 0.015 (S/cm2)      <0, 1e9>
        :cai(1)
}


STATE {
       c 
}


ASSIGNED {
        v (mV)
        celsius (degC)
        ek (mV)
        gkc (S/cm2)
        ik (mA/cm2)
        cai (1)
        cinf
        ctau (ms)
}


BREAKPOINT {
        SOLVE states METHOD cnexp
        gkc = gkcbar * c
        ik = gkc * xi(cai) * (v - ek) 
}


INITIAL {
        rates(v)
        c = cinf 
}


DERIVATIVE states {
        rates(v)
        c' = (cinf - c) / ctau
}


PROCEDURE rates(v(mV)) {
        LOCAL alpha, beta
        TABLE cinf, ctau DEPEND celcius FROM -100 TO 100 WITH 200 
        UNITSOFF

        if (v > 50) {
            alpha = 2 * exp((6.5 - v)/27)
            beta = 0.0
        } else {
            alpha = (exp((v - 10)/11) - exp((v - 6.5)/27)) / 18.975
            beta = 2 * exp((6.5 - v)/27) - alpha 
        }

        cinf = alpha / (alpha + beta)
        ctau = 1 / (alpha + beta)
}

FUNCTION xi(y) {
        if (y / 250.0 > 1.0) {
                xi = 1.0
        } else {
                xi = y / 250.0
        }
}


UNITSON
