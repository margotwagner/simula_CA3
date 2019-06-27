TITLE CA3_K_AHP.mod   potassium (AHP)

COMMENT
 This is the potassium afterhyperpolarization current, which has slower, calcium-activated dynamics. Found in the CA3.
  ("Intrinsic and Network Rhythmogensis in a Reduced Traub Model for CA3 Neurons" J.Comp.Neuro. (Neth.) 1:39-60 (1994).)
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
    (S) = (siemens)
}

? interface
NEURON {
        SUFFIX CA3_K_AHP
        USEION k READ ek WRITE ik
        RANGE gkbar, q, gk
        GLOBAL qinf, qtau
}
 
PARAMETER {	
        gkbar = .0008 (S/cm2)
        :cai (1)
}
 
 
ASSIGNED {
        celsius (degC)
        gk (S/cm2)
        v (mV)
        ek (mV)
        ik (mA/cm2)
        qinf
        qtau (ms)
        cai (1)
}

STATE {
    q
}
 
? currents
BREAKPOINT {
    SOLVE states METHOD cnexp
        gk = gkbar*q
    ik = gk*(v - ek)      
}
 
 
INITIAL {
	rates(cai)
	q = qinf
}

? states
DERIVATIVE states {  
        rates(cai)
        q' =  (qinf-q)/qtau
}
 

? rates

PROCEDURE rates(cai) {  :Computes rate and other constants at intracellular calcium level cai.
                      :Call once from HOC to initialize.
        LOCAL  alpha, beta, sum
        TABLE qinf, qtau DEPEND celsius FROM -100 TO 100 WITH 200
UNITSOFF
                :"q" potassium afterhyperpolarization activation system
                
        if (0.2e-4*cai < 0.01) {
            alpha = 0.2e-4*cai
            } else {
            alpha = 0.01
            }
            
        beta = 0.001
        sum = alpha + beta
        qtau = 1/(sum)
        qinf = alpha/sum
}

UNITSON

: test syntax in terminal using 'modlunit file.mod'
: compile in terminal using 'nrnivmodl file.mod'