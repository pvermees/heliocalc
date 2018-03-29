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

function append(smp,x,sx){
    var out = smp;
    if (x > 0 & sx > 0) {
	var extradig = 0;
	var lx = Math.max(1,Math.ceil(Math.log10(x)));
	var lsx = Math.floor(Math.log10(sx));
	var numdig = 0;
	if (lsx < 1) { numdig = lx - lsx + extradig; }
	out.push(x.toFixed(numdig),sx.toFixed(numdig));
    } else if (sx === null) {
	out.push(x.toFixed(3)); 
    } else {
	out.push('','');
    }
    return out;
}
