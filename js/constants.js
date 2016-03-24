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
    case 'H' : return 1.008;
    case 'O' : return 15.999;
    case 'F' : return 18.9984;
    case 'P' : return 30.973762;
    case 'Cl': return 35.45;
    case 'Ca': return 40.078;
    case 'Ti': return 47.867;
    case 'Si': return 28.0855;
    case 'S' : return 32.06;
    case 'Mn': return 54.9380443;
    case 'Fe': return 55.845;
    case 'Sr': return 87.621;
    case 'Y' : return 88.905842;
    case 'Zr': return 91.224;
    case 'Ba': return 137.327;
    case 'La': return 138.90547;
    case 'Ce': return 140.116;
    case 'Sm': return 150.36;
    case 'Hf': return 178.492;
    case 'Th': return 232.03806;
    case 'U': return 238.02891;
    case 'apatite': return 5*MM('Ca')+3*MM('P')+12*MM('O')+MM('F');
    case 'zircon' : return MM('Zr')+MM('Si')+4*MM('O');
    case 'titanite' : return MM('Ca')+MM('Ti')+MM('Si')+5*MM('O');
    case 'monazite' : return MM('Ce')/2+MM('La')/2+MM('P')+4*MM('O');
    case 'xenotime' : return MM('Y')+MM('P')+4*MM('O');
    case 'rutile' : return MM('Ti')+2*MM('O');
    case 'magnetite' : return 3*MM('Fe')+4*MM('O');
    case 'haematite' : return 2*MM('Fe')+3*MM('O');
    case 'goethite' : return MM('Fe')+2*MM('O')+MM('H');
    case 'barite' : return MM('Ba')+MM('S')+4*MM('O');
    default: return 0;
    }
}

// gramme of stoichiometric element per gramme of mineral
function gpg(mineral,element){
    var out = mpm(element,mineral)*MM(element)/MM(mineral);
    return out;
}

function ppm(mineral,element){
    return gpg(mineral,element)*1e6;
}

// moles of a stoichiometric element per mol of mineral
function mpm(element,mineral){
    switch (element){
    case 'Ca':
	switch (mineral){
	case 'apatite': return 5;
	case 'titanite': return 1;
	default: return 0;
	}
    case 'P':
	switch (mineral){
	case 'apatite': return 3;
	case 'monazite': return 1;
	case 'xenotime': return 1;
	default: return 0;
	}
    case 'Si':
	switch (mineral){
	case 'zircon': return 1;
	case 'titanite': return 1;
	default: return 0;
	}
    case 'Zr':
	switch (mineral){
	case 'zircon': return 1;
	default: return 0;
	}
    case 'Fe':
	switch (mineral){
	case 'magnetite': return 3;
	case 'haematite': return 2;
	case 'goethite': return 1;
	default: return 0;
	}
    case 'Ba':
	switch (mineral){
	case 'barite': return 1;
	default: return 0;
	}
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
