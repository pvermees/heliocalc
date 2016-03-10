function getNaturalIsotope(element){
    switch (element){
    case 'U': return 238;
    case 'Th': return 232;
    case 'Sm': return 147;
    default: return 0;
    }
}

function density(mineral){
    switch (mineral){
    case 'apatite': return 3.20;
    case 'zircon' : return 4.65;
    case 'titanite' : return 3.53;
    case 'monazite' : return 5.26;
    case 'xenotime' : return 4.75;
    case 'rutile' : return 4.25;
    case 'magnetite' : return 5.18;
    case 'haematite' : return 5.26;
    case 'goethite' : return 4.28;
    case 'barite' : return 4.5;
    default: return 0;
    }
}

// get Molar Mass of an element or mineral
function MM(stuff){
    switch (stuff){
    case 'O' : return 15.999;
    case 'F' : return 18.9984;
    case 'P' : return 30.973762;
    case 'Cl': return 35.45;
    case 'Ca': return 40.078;
    case 'Si': return 28.0855;
    case 'Mn': return 54.9380443;
    case 'Sr': return 87.621;
    case 'Y' : return 88.905842;
    case 'Zr': return 91.224;
    case 'Sm': return 150.36;
    case 'Hf': return 178.492;
    case 'Th': return 232.03806;
    case 'U': return 238.02891;
    case 'apatite': return 5*MM('Ca')+3*MM('P')+12*MM('O')+MM('F');
    case 'zircon' : return MM('Zr')+MM('Si')+4*MM('O');
    default: return 0;
    }
}

// gramme of stoichiometric element per gramme of mineral
function ppm(element){
    var mineral = (element === 'Zr' | element === 'Si') ? 'zircon' : 'apatite';
    var out = mpm(element)*MM(element)/MM(mineral);
    return out*1e6;
}

// mol of mineral per mol of a stoichiometric element
function mpm(element){
    switch (element){
    case 'Ca': return 5;
    case 'P': return 3;
    case 'Si': return 1;
    case 'Zr': return 1;
    default: return 0;
    }
}

function A(isotope){
    switch (isotope){
    case 28 : return 0.9223;
    case 29 : return 0.0467;
    case 30 : return 0.031;
    case 31 : return 1;
    case 35 : return 0.7577;
    case 37 : return 0.2423;
    case 40 : return 0.96941;
    case 42 : return 0.00647;
    case 43 : return 0.00135;
    case 44 : return 0.02086;
    case 46 : return 0.00004;
    case 47 : return 0.00187;
    case 55 : return 1;
    case 84 : return 0.0056;
    case 86 : return 0.0986;
    case 87 : return 0.07;
    case 88 : return 0.8258;
    case 89 : return 1;
    case 90 : return 0.5145;
    case 91 : return 0.1122;
    case 92 : return 0.1715;
    case 94 : return 0.1738;
    case 96 : return 0.028;
    case 144: return 0.0307;
    case 147: return 0.1499;
    case 148: return 0.1124;
    case 149: return 0.1382;
    case 150: return 0.0738;
    case 152: return 0.2675;
    case 154: return 0.2275;
    case 174: return 0.00162;
    case 176: return 0.05206
    case 177: return 0.18606;
    case 178: return 0.27297;
    case 179: return 0.13629;
    case 180: return 0.351;
    case 232: return 1;
    case 235: return 1/138.818;
    case 238: return 137.818/138.818;
    default: return 0;
    }
}

function lambda(isotope){
    switch (isotope){
    case 238: return 0.000155125;
    case 235: return 0.0009845841;
    case 232: return 0.00004933431;
    case 147: return 0.000006539;
    default: return 0;
    }
}

// get the alpha-travelling distances for Sm147, Th232, U235 and U238
function getR(mineral){
    switch (mineral){
    case "apatite": return [5.93,22.25,21.8,18.81];
    case "zircon": return [4.76,18.43,18.05,15.55];
    case "titanite": return [5.47,20.68,20.25,17.46];
    case "monazite": return [4.98,19.11,18.74,16.18];
    case "xenotime": return [4.68,17.99,17.63,15.20];
    case "rutile": return [4.77,18.14,17.76,15.30];
    case "magnetite": return [4.51,16.49,16.16,13.97];
    case "haematite": return [4.39,16.04,15.72,13.59];
    case "goethite": return [4.95,18.38,18,15.54];
    case "barite": return [5.54,21.50,21.05,18.14];
    default: return [0,0,0,0];
    }
}

function append(smp,x,sx){
    var out = smp;
    if (x > 0 & sx > 0) {
	var extradig = 0;
	var lx = Math.max(1,Math.ceil(Math.log10(x)));
	var lsx = Math.floor(Math.log10(sx));
	var numdig = 0;
	if (lsx < 1) { numdig = lx - lsx + extradig; }
	out.push(x.toFixed(numdig),sx.toFixed(numdig));
    } else if (x == 0 & sx === null) {
	out.push('');
    } else if (x == 0 & sx == 0) {
	out.push('','');
    } else if (x != 0 & sx === null) {
	out.push(x.toFixed(3)); 
    } else if (sx <= 0) {
	out.push(x,sx);
    }
    return out;
}

// uses Newton's method to solve equation 6.9 of Galbraith (2005)
function getZeta(zvz,mu){
    var zeta = 0;
    var dzeta = 0;
    var f, fprime, x;
    for (var i=0; i<10; i++){
	x = getwuzumu(zvz,zeta,mu);
	f = x.sumwu2zumu2 - x.sumwu;
	dzeta = f/x.fprime;
	zeta = zeta - dzeta;
	if (Math.abs(dzeta/zeta) < 0.001) { 
	    return zeta;
	}
    }
    // the algorithm didn't converge!
    return null;
}

function getMu(zvz,zeta){
    var mu = zvz[0][0]; // won't be used
    var x = getwuzumu(zvz,zeta,mu);
    return (x.sumwuzu/x.sumwu);
}

function getwuzumu(zvz,zeta,mu){
    var wu, a = 0, b = 0, c = 0, d = 0;
    for (var i = 0; i<zvz.length; i++){
	wu = 1/(zvz[i][1]+zeta);
	a += wu*zvz[i][0];
	b += wu;
	c += Math.pow(wu*(zvz[i][0]-mu),2);
	d += 1/(zeta+zvz[i][1]) - 0.5*Math.pow(zvz[i][0]-mu,2)/Math.pow(zeta+zvz[0][1],3);
    }
    return {sumwuzu: a, sumwu: b, sumwu2zumu2: c, fprime: d};
}

function age(U,Th,Sm,He,Ft){
    var A147 = A('Sm',147);
    var P = (8*Ft[3]*lambda(238)*137.818/138.818 + 7*Ft[2]*lambda(235)/138.818)*U + 
	6*Ft[1]*lambda(232)*Th + A147*Ft[0]*lambda(147)*Sm;
    var L = ((8*Ft[3]*lambda(238)*lambda(238)*137.818/138.818 + 7*Ft[2]*lambda(235)*lambda(235)/138.818)*U + 
             6*Ft[1]*Th*lambda(232)*lambda(232) + A147*Ft[0]*Sm*lambda(147)*lambda(147)) / P;
    var t = Math.log(1 + L*He/P)/L;
    return t;
}

function err_age(data,i,doFt){
    var J = data.Jt(i,doFt);
    var E = data.Sigma_t(i);
    var t_err = propagate(J,E);
    return t_err;
}

function propagate(J,E){
    var n = J.length, JEi = 0, out = 0;
    for (var i=0; i<n; i++){
	JEi = 0;
	for (var j=0; j<n; j++){ JEi += J[j]*E[i][j]; }
	out += J[i]*JEi;
    }
    return Math.sqrt(out);
}

function getCSMXAsMSAx(X,x,S,s){
    return ppm(S)*(MM(X)*A(s))/(MM(S)*A(x));
}

function dX8dt(t){
    return 8 * Math.exp(lambda(238) * t) * lambda(238) * 137.818 / 138.818;
}

function X8(t){
    return 8 * (Math.exp(lambda(238) * t) - 1) * 137.818 / 138.818;
}

function dX5dt(t){
    return 7 * Math.exp(lambda(235) * t) * lambda(235) * 137.818 / 138.818;
}

function X5(t){
    return 7 * (Math.exp(lambda(235) * t) - 1) / 138.818;
}

function dX2dt(t){
    return 6 * Math.exp(lambda(232) * t) * lambda(232);
}

function X2(t){
    return 6 * (Math.exp(lambda(232) * t) - 1);
}

function dX7dt(t){
    return 0.1499 * Math.exp(lambda(147) * t) * lambda(147);
}

function X7(t){
    return 0.1499 * (Math.exp(lambda(147) * t) - 1);
}

function JB2(IX,IU){
    var CS = ppm(IU.S);
    var xesm = IX.get_xesm();
    var sesm = IU.get_sesm();
    var CSMXAsMSAx = getCSMXAsMSAx(IX.X,IX.x,IU.S,IU.s);
    var dCXdxesm = CSMXAsMSAx / sesm;
    var dCXdsesm = -CSMXAsMSAx * xesm / Math.pow(sesm,2);
    var dCXdyxmsmsp = dCXdxesm*IX.dxesm_dyxmsmsp();
    var dCXdsxmsmsp = dCXdsesm*IU.dsesm_dsxmsmsp();
    var dCXdxsmss   = dCXdsesm*IU.dsesm_dxsmss();
    return [dCXdyxmsmsp,dCXdsxmsmsp,dCXdxsmss];
}

function SigmaB2(IX,IU){
    var var_yxmsmsp = Math.pow(IX.err_yxmsmsp,2);
    var var_sxmsmsp = Math.pow(IU.err_sxmsmsp,2);
    var var_xsmss = IU.var_xsmss;
    var EB = [[var_yxmsmsp,       0          ,       0        ],
	      [0          ,   var_sxmsmsp    ,       0        ],
	      [0,                 0,                 var_xsmss]];
    return EB;
}
