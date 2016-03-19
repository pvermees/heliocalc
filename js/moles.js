get_moles = function(dat,X,i){
    var x = dat.get_x(X);
    if (X === 'He'){
	var m = 0, err_m = 0;
	m = dat.get_signal('He',4,i) / 22.414;
	err_m = dat.get_error('He',4,i) / 22.414;
	return [m,err_m];
    }
    if (dat.inSpike(X)){
	var Ax = A(x);
	var xesm = get_xesm(dat,X,x,i)/Ax;
	var err_xesm = get_err_xesm(dat,X,x,i)/Ax;
	return [xesm,err_xesm];
    }
    if (dat.inStand(X) & dat.inSpike('U')){
	var Ax = A(x);
	var zesm = get_zesm(dat,X,x,i)/Ax;
	var err_zesm = get_err_zesm(dat,X,x,i)/Ax;
	return [zesm,err_zesm];
    }
    return [0,0];
}

// equation 8
get_xesm = function(dat,X,x,i){
    var y = dat.get_y(X);
    var ytsp = dat.get_ytsp(X,i);
    var yxtsp = dat.get_yxtsp(X);
    var xytsp = 1/yxtsp;
    var yxtsm = A(y)/A(x);
    var yxmsmsp = dat.get_signal(X,y,i)/dat.get_signal(X,x,i);
    return ytsp*xytsp*(yxtsp-yxmsmsp)/(yxmsmsp-yxtsm);
}

get_zesm = function(dat,Z,z,i){
    var U8esm = get_xesm(dat,'U',238,i);
    var utsp = dat.get_ytsp('U',i);
    var U8utsp = 1/dat.get_yxtsp('U');
    var z8msmsp = dat.get_signal(Z,z,i)/dat.get_signal('U',238,i);
    var z8tss = dat.get_zxtss(Z,z,'U',238);
    var U8zmss = dat.get_zxmss('U',238,Z,z)[0];
    return (U8esm+utsp*U8utsp)*z8msmsp*z8tss*U8zmss;
}
