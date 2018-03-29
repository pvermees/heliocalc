// Using Meesters and Dunai (2005)
function get_age(U,Th,Sm,He,Ft){
    U = Math.max(0,U);
    Th = Math.max(0,Th);
    Sm = Math.max(0,Sm);
    He = Math.max(0,He);
    var A147 = A('Sm',147);
    var P = (8*Ft[3]*lambda(238)*137.818/138.818 + 7*Ft[2]*lambda(235)/138.818)*U + 
	6*Ft[1]*lambda(232)*Th + A147*Ft[0]*lambda(147)*Sm;
    var L = ((8*Ft[3]*lambda(238)*lambda(238)*137.818/138.818 +
	      7*Ft[2]*lambda(235)*lambda(235)/138.818)*U + 
             6*Ft[1]*Th*lambda(232)*lambda(232) +
	     A147*Ft[0]*Sm*lambda(147)*lambda(147)) / P;
    var t = Math.log(1 + L*He/P)/L;
    return t;
}
