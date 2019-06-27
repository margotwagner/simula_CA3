TITLE CA3_I_Ca.mod   calcium channel
 
COMMENT
 The inward fast calcium current I_Ca into the dendrite in the CA3 model. 
  ("Intrinsic and Network Rhythmogensis in a Reduced Traub Model for CA3 Neurons" J.Comp.Neuro. (Neth.) 1:39-60 (1994).)
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}

? interface
NEURON {
        SUFFIX CA3_I_Ca
        USEION ca READ eca WRITE ica
        RANGE gcabar, s, gca
        GLOBAL sinf, stau
        THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gcabar = 0.01 (S/cm2)	
}
 
STATE {
        s
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        eca (mV)
        gca (S/cm2)
        ica (mA/cm2)
        sinf 
        stau (ms) 
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gca = gcabar*s*s
	ica = gca*(v - eca)
}
 
 
INITIAL {
	rates(v)
	s = sinf
}

? states
DERIVATIVE states {  
        rates(v)
        s' =  (sinf-s)/stau

}
 

? rates

PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
        TABLE sinf, stau DEPEND celsius FROM -100 TO 100 WITH 200
UNITSOFF
                :"s" calcium activation system
        alpha = 1.6 /  (1 + exp(-0.072*(v-65)))
        beta = 0.2 * vtrap(v-51.1,5)
        
        sum = alpha + beta
        stau = 1/(sum)
        sinf = alpha/sum
}
        
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON

: test syntax in terminal using 'modlunit file.mod'
: compile in terminal using 'nrnivmodl file.mod'