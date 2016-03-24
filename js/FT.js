function FT(R,l,w,h,habit) {
    out = []; 	// Sm147, Th232, U235, U238, 
    var F = 0;
    if (l <= 0) { return [1,1,1,1]; }
    else if (w <= 0) { w = l; }
    if (h <= 0) { h = w; }
    for (var j=0; j<4; j++) {
	switch (habit){
	case "ellipsoid":
	    F = FT_ellipsoid(R[j],l,w,h); break;
	case "cylinder":
	    F = FT_cylinder(R[j],l,w,h); break;
	case "cylinder (inf)":
	    F = FT_inf_cylinder(R[j],l,w,h); break;
	case "tetragonal (0T)":
	    F = FT_tetragonal(R[j],l,w,h,0); break;
	case "tetragonal (1T)":
	    F = FT_tetragonal(R[j],l,w,h,1); break;
	case "tetragonal (2T)":
	    F = FT_tetragonal(R[j],l,w,h,2); break;
	case "tetragonal (inf)":
	    F = FT_inf_tetragonal(R[j],l,w,h); break;
	case "hexagonal (0T)":
	    F = FT_hexagonal(R[j],l,w,h,0); break;
	case "hexagonal (1T)":
	    F = FT_hexagonal(R[j],l,w,h,1); break;
	case "hexagonal (2T)":
	    F = FT_hexagonal(R[j],l,w,h,2); break;
	case "hexagonal (inf)":
	    F = FT_inf_hexagonal(R[j],l,w,h); break;
	}
	out.push(F);
    }
    return out;
}

function FT_ellipsoid(R,c,b,a){
    var V = Math.PI*a*b*c*4/3;
    var p = 1.6075;
    var S = 4*Math.PI*Math.pow((Math.pow(a,p)*Math.pow(b,p) + 
				Math.pow(b,p)*Math.pow(c,p) + 
				Math.pow(c,p)*Math.pow(a,p))/3,1/p);
    var Rs = 3*V/S;
    var out = 1 - 0.75*R/Rs + (1/16 + 0.1686*Math.pow(1-a/Rs,2))*Math.pow(R/Rs,3);
    return out;
}

function FT_cylinder(R,l,w,h){
    var r = Math.sqrt(w*h/Math.PI);
    var out = 1 - (r+l)*R/(2*r*l) + 0.2122*R*R/(r*l) + 0.0153*Math.pow(R/r,3);
    return out;
}

function FT_inf_cylinder(R,l,w,h){
    var r = Math.sqrt(w*h/Math.PI);
    var out = 1 - 0.5*R/r + 0.0153*Math.pow(R/r,3);
    return out;
}

function FT_tetragonal(R,c,b,a,Np){
    var Rs, out;
    if (Np == 0) {
	Rs = 1.5/(1/a + 1/b + 1/c);
	out = 1 - 0.75*R/Rs + 0.2122*(a+b+c)*R*R/(a*b*c) - 0.00995*Math.pow(R,3)/(a*b*c);
    } else {
	var V = VOL_tetragonal(c,b,a,Np);
	var S = 2*(a*b + b*c + a*c) - Np*((a*a + b*b)/2 + (2-Math.sqrt(2))*a*b);
	Rs = 3*V/S;
	out = 1 - 0.75*R/Rs + (0.2095*(a+b+c) - (0.096-0.013*(a*a+b*b)/(c*c))*(a+b)*Np)*R*R/V;
    }
    return out;
}

function FT_inf_tetragonal(R,c,b,a){
    Rs = 1.5/(1/a + 1/b);
    return 1 - 0.75*R/Rs + 0.2122*R*R/(a*b);
}

function FT_hexagonal(R,H,W,L,Np){
    var V = VOL_hexagonal(H,W,L,Np);
    var S = 2*H*(W + L/Math.sqrt(3)) + 2*L*(W-L/(2*Math.sqrt(3))) -
	Np*(Math.sqrt(3)*W*W/4 + (2 - Math.sqrt(2))*W*L + (Math.sqrt(2)-1)*L*L/(2*Math.sqrt(3)));
    var Rs = 3*V/S;
    var out = 1 - 3*R/(4*Rs) + (
        (0.2093 - 0.0465*Np)*(W+L/Math.sqrt(3)) + 
	(0.1062 + 0.2234*R/(R + 6*(W*Math.sqrt(3)-L))) *
	(H - Np*(W*Math.sqrt(3)/2 + L)/4)
                               )*R*R/V;
    return out;
}

function FT_inf_hexagonal(R,H,W,L){
    var V = VOL_hexagonal(H,W,L,0);
    var S = 2*H*(W + L/Math.sqrt(3)) + 2*L*(W-L/(2*Math.sqrt(3)));
    var Rs = 3*V/S;
    var out = 1 - 3*R/(4*Rs) + ((0.1062 + 0.2234*R/(R + 6*(W*Math.sqrt(3)-L))) * H )*R*R/V;
    return out;
}

function VOL(l,w,h,habit){
    if (l <= 0) { return 0; }
    else if (w <= 0) { w = l; }
    if (h <= 0) { h = w; }
    switch (habit){
    case "ellipsoid":
	return VOL_ellipsoid(l/2,w/2,h/2);
    case "cylinder":
	return VOL_cylinder(l,w,h);
    case "cylinder (inf)":
	return VOL_cylinder(l,w,h);
    case "tetragonal (0T)":
	return VOL_tetragonal(l,w,h,0);
    case "tetragonal (1T)":
	return VOL_tetragonal(l,w,h,1);
    case "tetragonal (2T)":
	return VOL_tetragonal(l,w,h,2);
    case "tetragonal (inf)":
	return VOL_tetragonal(l,w,h,0);
    case "hexagonal (0T)":
	return VOL_hexagonal(l,w,h,0);
    case "hexagonal (1T)":
	return VOL_hexagonal(l,w,h,1);
    case "hexagonal (2T)":
	return VOL_hexagonal(l,w,h,2);
    case "hexagonal (inf)":
	return VOL_hexagonal(l,w,h,0);
    }
    return l*w*h;
}

function VOL_ellipsoid(c,b,a){
    return Math.PI*a*b*c*4/3;
}

function VOL_cylinder(l,w,h){
    var r = Math.sqrt(w*h/4);
    return Math.PI*r*r*l;
}

function VOL_tetragonal(c,b,a,Np){
    var ap = Math.min(a,b);
    var bp = Math.max(a,b);
    return a*b*c - (bp*bp + ap*ap/3)*Np*ap/4;
}

// using the parameter names of Ketcham et al. (2011)
function VOL_hexagonal(H,W,L,Np){
    var dV = 0;
    var diff = L - Math.sqrt(3)*W/2;
    if (diff > 0) { dV = Math.pow(diff,3)/(6*Math.sqrt(3)); }
    var V = H*L*(W-L/(2*Math.sqrt(3))) - Np*(Math.sqrt(3)*W*W*L/8 - dV)
    return V
}
