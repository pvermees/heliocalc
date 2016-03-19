get_zxmss = function(dat,ne,ni,de,di){
    var n = get_xss(dat,ne,ni);
    var d = get_xss(dat,de,di);
    var average = n[0]/d[0];
    var variance = average*average*(n[1]/(n[0]*n[0]) + d[1]/(d[0]*d[0]));
    return [average,variance];
}

// used by get_k(dat,i)
get_cov_s7mss_8smss = function(dat,i){
    var S = dat.get_S();
    var s = dat.get_s(S);
    var A = get_xss(dat,S,s,i);
    var B = get_xss(dat,'Sm',147,i);
    var C = get_xss(dat,'U',238,i);
    return get_cov(A[0],Math.sqrt(A[1]),B[0],C[0]);
}

get_xss = function(dat,X,x){
    var snames = dat.names();
    var sum = 0, avg = 0, n = 0, ndmss = 0, mu;
    var zvz = [], wu = [];
    var z = 0, vz = 0, variance = 0;
    var zeta = 0, mu = 0, oldzeta = 0, oldmu = 0;
    for (var i = 0; i<snames.length; i++){
	if (dat.isStand(i)) {
	    z = Math.log(dat.get_signal(X,x,i));
	    vz = Math.pow(dat.get_error(X,x,i)/dat.get_signal(X,x,i),2);
	    zvz.push([z,vz]);
	}
    }
    if (zvz.length < 2){ // when there is only one standard
	avg = Math.exp(z);
	variance = avg*avg*vz;
	return [avg,variance]; 
    }
    for (var j = 0; j<10; j++){
	mu = get_mu(zvz,zeta);
	zeta = get_zeta(zvz,mu);
	if (zeta === null) {
	    mu = get_mu(zvz,0);
	    break;
	} else if ((Math.abs((mu-oldmu)/mu) < 0.001) & (Math.abs((zeta-oldzeta)/zeta) < 0.001)) {
	    break;
	} else {
	    oldmu = mu;
	    oldzeta = zeta;
	}
    }
    avg = Math.exp(mu);
    if (zeta === null){ // if the algorithm didn't converge, use the sample variance
	var n = zvz.length;
	for (var i=0; i<n; i++){
	    variance += Math.pow(zvz[i][0]-mu,2)/(n*(n-1));
	}
	variance *= avg*avg;
    } else {
	var wuzumu = get_wuzumu(zvz,zeta,mu);
	variance = avg*avg/wuzumu.sumwu;
    }
    return [avg,variance]; 
}

// uses Newton's method to solve equation 6.9 of Galbraith (2005)
function get_zeta(zvz,mu){
    var zeta = 0;
    var dzeta = 0;
    var f, fprime, x;
    for (var i=0; i<10; i++){
	x = get_wuzumu(zvz,zeta,mu);
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

function get_mu(zvz,zeta){
    var mu = zvz[0][0]; // won't be used
    var x = get_wuzumu(zvz,zeta,mu);
    return (x.sumwuzu/x.sumwu);
}

function get_wuzumu(zvz,zeta,mu){
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
