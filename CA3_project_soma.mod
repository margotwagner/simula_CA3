TITLE CA3_project_dend.mod   calcium, potassium (AHP), potassium (C), synaptic and leak channels
 
COMMENT
 This is the single soma compartment Hodgkin-Huxley treatment for the set of sodium, 
  potassium, and leakage channels found in the CA3.
  ("Intrinsic and Network Rhythmogensis in a Reduced Traub Model for CA3 Neurons" J.Comp.Neuro. (Neth.) 1:39-60 (1994).)
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}

? interface
NEURON {
        SUFFIX CA3_project_soma
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, m, h, n, gl, el, gna, gk
        GLOBAL minf, hinf, ninf, mtau, htau, ntau
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .03 (S/cm2)	
        gkbar = .015 (S/cm2)	
        gl = .0001 (S/cm2)	
        el = 0 (mV)
}
 
STATE {
        m h n
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

	gna (S/cm2)
	gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf
	mtau (ms) htau (ms) ntau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*minf*minf*h
	ina = gna*(v - ena)
        gk = gkbar*n
	ik = gk*(v - ek)      
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}
 

? rates

PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
        TABLE minf, mtau, hinf, htau, ninf, ntau DEPEND celsius FROM -100 TO 100 WITH 200
UNITSOFF
        alpha = .32 * vtrap(13.1-v,4) :(13.1 - v)/(exp((13.1-v)/4)-1) :"m" sodium activation system
        beta =  .28 * vtrap(v-40.1,5) :(v - 40.1)/(exp((v-40.1)/5)-1)
        sum = alpha + beta
	mtau = 1/(sum)
        minf = alpha/sum
        
                :"h" sodium inactivation system
        alpha = 0.128*exp((17-v)/18)
        beta = 4/(1+exp((40 - v)/5))
        sum = alpha + beta
	htau = 1/(sum)
        hinf = alpha/sum
                :"n" potassium - dr activation system
        alpha = .016 * vtrap(35.1-v,5) :(35.1 - v)/(exp((35.1-v)/5)-1)
        beta = 0.25*exp(0.5 - .025*v)
        sum = alpha + beta
        ntau = 1/(sum)
        ninf = alpha/sum
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