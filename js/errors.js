propagate = function(Jac,Cov){
    var J = math.matrix(Jac);
    var E = math.matrix(Cov);
    var Jt =  math.transpose(J);
    var EJ = math.multiply(E,J);
    return math.multiply(Jt,EJ).valueOf();
}

get_err_yxmsmsp = function(dat,X,x,Y,y,i){
    return Math.sqrt(get_var_yxmsmsp(dat,X,x,Y,y,i));
}

// equation 2
get_var_yxmsmsp = function(dat,Y,y,X,x,i){
    var n = dat.get_signal(Y,y,i); // numerator
    var d = dat.get_signal(X,x,i); // denominator
    var err_d = dat.get_error(X,x,i);
    var err_n = dat.get_error(Y,y,i);
    var variance = Math.pow(n/d,2)*(Math.pow(err_n/n,2)+Math.pow(err_d/d,2));
    return variance;
}

// equation 13
get_err_xesm = function(dat,X,x,i){
    var y = dat.get_y(X);
    var err_yxmsmsp = get_err_yxmsmsp(dat,X,y,X,x,i);
    var dxesm_dyxmsmsp = get_dxesm_dyxmsmsp(dat,X,x,y,i);
    return Math.abs(err_yxmsmsp * dxesm_dyxmsmsp);
}

// equation 15
get_err_zesm = function(dat,Z,z,i){
    var J = JA(dat,Z,z,i);
    var E = EA(dat,Z,z,i);
    var V = propagate(J,E);
    return Math.sqrt(V);
}

get_err_CX = function(dat,X,x,i){
    var J = [], E = [], V = 0;
    if (dat.inStand(X) & !dat.inSpike(X)){
	J = JB(dat,X,x,i);
	E = EB(dat,X,x,i);
    } else if (X === 'U' & dat.standardson()){
	J = JC(dat,i);
	E = EA(dat,dat.get_S(),dat.get_s(),i);
    } else if (dat.inSpike(X) & dat.standardson()){
	J = JD(dat,X,x,i);
	E = ED(dat,X,x,i);
    } else if (dat.inSpike(X)){
	var err_xesm = get_err_xesm(dat,X,x,i);
	var mass = dat.get_mass(i)[0];
	var err = (mass > 0) ? (MM(X)*err_xesm)/(A(x)*mass) : 0 ;
	return err;
    }
    V = propagate(J,E);
    return Math.sqrt(V);
}

get_err_age_radialplotter = function(dat,t,Ft,i){
    var J = JH(dat,t,Ft,i);
    var E = EH(dat,i);
    var V = propagate(J,E);
    return Math.sqrt(V);
}

get_covmat_age_eU = function(dat,t,Ft,i){
    var J = JG(dat,t,Ft,i);
    var E = EG(dat,i);
    return propagate(J,E);
}

// covmat is a 2x2 covariance matrix
cov2cor = function(covmat){
    var cov = covmat[0][1];
    var ss = Math.sqrt(covmat[0][0]*covmat[1][1]);
    var out = (ss > 0) ? cov/ss : 0 ;
    return out;
}
