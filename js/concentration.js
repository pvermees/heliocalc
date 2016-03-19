get_conc = function(dat,X,i){
    var x = dat.get_x(X);
    var CX = get_CX(dat,X,x,i);
    var err_CX = get_err_CX(dat,X,x,i);;
    return [CX,err_CX];
}

get_CX = function(dat,X,x,i){
    if (dat.standardson()){
	var S = dat.get_S();
	var s = dat.get_s();
	var CS = ppm(dat.mineral(),S);
	var xsesm = get_xsesm(dat,X,x,S,s,i);
	return (MM(X)*A(s)*xsesm*CS)/(MM(S)*A(x));
    } else {
	var xesm = get_xesm(dat,X,x,i);
	var mass = dat.get_mass(i)[0];
	var CX = (mass > 0) ? xesm*MM(X)/(A(x)*mass) : 0 ;
	return CX;
    }
}

get_xsesm = function(dat,X,x,S,s,i){
    var xsesm = 0;
    if (dat.inStand(X) & !dat.inSpike(X)){
	var xsmsm = dat.get_signal(X,x,i)/dat.get_signal(S,s,i);
	var xstss = dat.get_zxtss(X,x,S,s);
	var sxmss = dat.get_zxmss(S,s,X,x)[0];
	xsesm = xsmsm*xstss*sxmss;
    } else if (dat.inSpike(X)){
	var xesm = dat.get_xesm(X,x,i);
	var sesm = dat.get_zesm(S,s,i);
	xsesm = xesm/sesm;
    } else {
	// do nothing
    }
    return xsesm;
}
