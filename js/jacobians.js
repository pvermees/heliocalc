JA = function(dat,Z,z,i){
    var dzesm_du8msmsp = get_dzesm_du8msmsp(dat,Z,z,i);
    var dzesm_dz8msmsp = get_dzesm_dz8msmsp(dat,Z,z,i);
    var dzesm_d8zmss = get_dzesm_d8zmss(dat,Z,z,i);
    return [[dzesm_du8msmsp],[dzesm_dz8msmsp],[dzesm_d8zmss]];
}

JB = function(dat,X,x,i){
    var dCX_dxsmsm = get_dCX_dxsmsm(dat,X,x,i);
    var dCX_dsxmss = get_dCX_dsxmss(dat,X,x,i);
    return [[dCX_dxsmsm],[dCX_dsxmss]];
}

JC = function(dat,i){
    var dCU_du8msmsp = get_dCU_du8msmsp(dat,i);
    var dCU_ds8msmsp = get_dCU_ds8msmsp(dat,i);
    var dCU_d8smss = get_dCU_d8smss(dat,i);
    return [[dCU_du8msmsp],[dCU_ds8msmsp],[dCU_d8smss]];
}

JD = function(dat,X,x,i){
    var dCX_du8msmsp = get_dCX_du8msmsp(dat,X,x,i);
    var dCX_ds8msmsp = get_dCX_ds8msmsp(dat,X,x,i);
    var dCX_d8smss = get_dCX_d8smss(dat,X,x,i);
    var dCX_dyxmsmsp = get_dCX_dyxmsmsp(dat,X,x,i);
    return [[dCX_du8msmsp],[dCX_ds8msmsp],[dCX_d8smss],[dCX_dyxmsmsp]];
}

JG = function(dat,t,Ft,i){
    var dt_du8msmsp = get_dt_du8msmsp(dat,t,Ft,i);
    var deU_du8msmsp = get_deU_du8msmsp(dat,i);
    var deU_ds8msmsp = get_deU_ds8msmsp(dat,i);
    var deU_d8smss = get_deU_d8smss(dat,i);
    var dt_dv2msmsp = get_dt_dv2msmsp(dat,t,Ft,i);
    var deU_dv2msmsp = get_deU_dv2msmsp(dat,i);
    var dt_dw7msmsp = get_dt_dw7msmsp(dat,t,Ft,i);
    var deU_dw7msmsp = get_deU_dw7msmsp(dat,i);
    var deU_d7smsmsp = get_deU_d7smsmsp(dat,i);
    var deU_ds7mss = get_deU_ds7mss(dat,i);
    var dt_d78msmsp = get_dt_d78msmsp(dat,t,Ft,i);
    var dt_d87mss = get_dt_d87mss(dat,t,Ft,i);
    var dt_d4esm = get_dt_d4esm(dat,t,Ft,i);
    var J = [[dt_du8msmsp,deU_du8msmsp],
	     [0,deU_ds8msmsp],
	     [0,deU_d8smss],
	     [dt_dv2msmsp,deU_dv2msmsp],
	     [dt_dw7msmsp,deU_dw7msmsp],
	     [0,deU_d7smsmsp],
	     [0,deU_ds7mss],
	     [dt_d78msmsp,0],
	     [dt_d87mss,0],
	     [dt_d4esm,0]];
    return J;
}

// used for age uncertainty of radialplotter input
JH = function(dat,t,Ft,i){
    var dt_du8msmsp = get_dt_du8msmsp(dat,t,Ft,i);
    var dt_dv2msmsp = get_dt_dv2msmsp(dat,t,Ft,i);
    var dt_dw7msmsp = get_dt_dw7msmsp(dat,t,Ft,i);
    var dt_d78msmsp = get_dt_d78msmsp(dat,t,Ft,i);
    var dt_d4esm = get_dt_d4esm(dat,t,Ft,i);
    var J = [[dt_du8msmsp],[dt_dv2msmsp],[dt_dw7msmsp],[dt_d78msmsp],[dt_d4esm]];
    return J;
}

get_dt_du8msmsp = function(dat,t,Ft,i){
    var dD_du8msmsp = get_dD_du8msmsp(dat,t,Ft,i); // equation 85
    var dD_dt = get_dD_dt(dat,t,Ft,i); // equation 80
    return -dD_du8msmsp/dD_dt;
}

get_dt_dv2msmsp = function(dat,t,Ft,i){
    var dD_dv2msmsp = get_dD_dv2msmsp(dat,t,Ft,i); // equation 88
    var dD_dt = get_dD_dt(dat,t,Ft,i); // equation 80
    return -dD_dv2msmsp/dD_dt;
}

get_dt_dw7msmsp = function(dat,t,Ft,i){
    var dD_dw7msmsp = get_dD_dw7msmsp(dat,t,Ft,i); // equation 89
    var dD_dt = get_dD_dt(dat,t,Ft,i); // equation 80
    return -dD_dw7msmsp/dD_dt;
}

get_dt_d78msmsp = function(dat,t,Ft,i){
    var dD_d78msmsp = get_dD_d78msmsp(dat,t,Ft,i); // equation 86
    var dD_dt = get_dD_dt(dat,t,Ft,i); // equation 80
    return -dD_d78msmsp/dD_dt;
}

get_dt_d87mss = function(dat,t,Ft,i){
    var dD_d87mss = get_dD_d87mss(dat,t,Ft,i); // equation 87
    var dD_dt = get_dD_dt(dat,t,Ft,i); // equation 80
    return -dD_d87mss/dD_dt;
}

get_dt_d4esm = function(dat,t,Ft,i){
    var dD_d4esm = -1; // equation 90
    var dD_dt = get_dD_dt(dat,t,Ft,i); // equation 80
    return -dD_d4esm/dD_dt;
}

// equation 14
get_dxesm_dyxmsmsp = function(dat,X,x,y,i){
    var ytsp = dat.get_ytsp(X,i);
    var yxtsm = A(y)/A(x);
    var xytsp = 1/dat.get_yxtsp(X);
    var yxmsmsp = dat.get_signal(X,y,i)/dat.get_signal(X,x,i);
    return ytsp*(yxtsm*xytsp-1)/Math.pow(yxmsmsp-yxtsm,2);
}

// equation 32
get_dCX_dxsmsm = function(dat,X,x,i){
    var S = dat.get_S();
    var s = dat.get_s();
    var xstss = dat.get_zxtss(X,x,S,s,i);
    var sxmss = dat.get_zxmss(S,s,X,x,i)[0];
    return (MM(X)*A(s)*dat.CS(S)*xstss*sxmss)/(MM(S)*A(x));
}

// equation 33
get_dCX_dsxmss = function(dat,X,x,i){
    var S = dat.get_S();
    var s = dat.get_s();
    var xstss = dat.get_zxtss(X,x,S,s,i);
    var xsmsm = dat.get_signal(X,x,i)/dat.get_signal(S,s,i);
    return (MM(X)*A(s)*dat.CS(S)*xstss*xsmsm)/(MM(S)*A(x));
}

// equation 58
get_dzesm_du8msmsp = function(dat,Z,z,i){
    var u = dat.get_y('U');
    var dzesm_d8esm = get_dzesm_d8esm(dat,Z,z,i); // equation 59
    var dxesm_dyxmsmsp = get_dxesm_dyxmsmsp(dat,'U',238,u,i); // equation 14
    return dzesm_d8esm*dxesm_dyxmsmsp;
}

// equation 59
get_dzesm_d8esm = function(dat,Z,z,i){
    var z8msmsp = dat.get_signal(Z,z,i)/dat.get_signal('U',238,i);
    var z8tss = dat.get_zxtss(Z,z,'U',238);
    var U8zmss = 1/dat.get_zxmss(Z,z,'U',238)[0];
    return z8msmsp*z8tss*U8zmss;
}

// equation 60
get_dzesm_dz8msmsp = function(dat,Z,z,i){
    var U8esm = dat.get_xesm('U',238,i);
    var utsp = dat.get_ytsp('U',i);
    var U8utsp = 1/dat.get_yxtsp('U',i);
    var z8tss = dat.get_zxtss(Z,z,'U',238);
    var U8zmss = 1/dat.get_zxmss(Z,z,'U',238)[0];
    return (U8esm+utsp*U8utsp)*U8zmss*z8tss;
}

// equation 61
get_dzesm_d8zmss = function(dat,Z,z,i){
    var U8esm = dat.get_xesm('U',238,i);
    var utsp = dat.get_ytsp('U',i);
    var U8utsp = 1/dat.get_yxtsp('U',i);
    var z8msmsp = dat.get_signal(Z,z,i)/dat.get_signal('U',238,i);
    var z8tss = dat.get_zxtss(Z,z,'U',238);
    return (U8esm+utsp*U8utsp)*z8msmsp*z8tss;
}

// equation 62
get_dCU_du8msmsp = function(dat,i){
    var dCU_d8esm = get_dCU_d8esm(dat,i); // equation 63
    var u = dat.get_y('U');
    var d8esm_du8msmsp = get_dxesm_dyxmsmsp(dat,'U',238,u,i); // equation 14
    return dCU_d8esm*d8esm_du8msmsp;
}

// equation 63
get_dCU_d8esm = function(dat,i){
    if (dat.standardson()){
	var S = dat.get_S();
	var s = dat.get_s();
	var sesm = dat.get_zesm(S,s,i); // equation 12
	return (MM('U')*A(s)*dat.CS(S))/(MM(S)*A(238)*sesm);
    } else {
	var mass = dat.get_mass(i)[0];
	var out = (mass > 0) ? MM('U')/(A(238)*mass) : 0 ;
	return out;
    }
}

// equation 64
get_dCU_ds8msmsp = function(dat,i){
    if (dat.standardson()){
	var S = dat.get_S();
	var s = dat.get_s();
	var dCU_dsesm = get_dCU_dsesm(dat,S,s,i); // equation 65
	var dsesm_ds8msmsp = get_dzesm_dz8msmsp(dat,S,s,i); // equation 14
	return dCU_dsesm*dsesm_ds8msmsp;
    } else {
	return 0;
    }
}

// equation 65
get_dCU_dsesm = function(dat,S,s,i){
    var U8esm = dat.get_xesm('U',238,i); // equation 8
    var sesm = dat.get_zesm(S,s,i); // equation 12
    return -(MM('U')*A(s)*U8esm*dat.CS(S))/(MM(S)*A(238)*sesm*sesm);
}

// equation 66
get_dCU_d8smss = function(dat,i){
    if (dat.standardson()){
	var S = dat.get_S();
	var s = dat.get_s();
	var dCU_dsesm = get_dCU_dsesm(dat,S,s,i); // equation 65
	var dsesm_d8smss = get_dzesm_d8zmss(dat,S,s,i); // equation 61
	return dCU_dsesm*dsesm_d8smss;
    } else {
	return 0;
    }
}

// equation 67
get_dCX_du8msmsp = function (dat,X,x,i){
    if (dat.standardson()){
	var S = dat.get_S();
	var s = dat.get_s();
	var dCX_dsesm = get_dCX_dsesm(dat,S,s,X,x,i); // equation 68
	var dsesm_du8msmsp = get_dzesm_du8msmsp(dat,S,s,i); // equation 58
	return dCX_dsesm*dsesm_du8msmsp;
    } else {
	return 0;
    }
}

// equation 68
get_dCX_dsesm = function(dat,S,s,X,x,i){
    var xesm = dat.get_xesm(X,x,i); // equation 8
    var sesm = dat.get_zesm(S,s,i); // equation 12
    return -(MM(X)*A(s)*xesm*dat.CS(S))/(MM(S)*A(x)*sesm*sesm);
}

// equation 69
get_dCX_ds8msmsp = function(dat,X,x,i){
    if (dat.standardson()){
	var S = dat.get_S();
	var s = dat.get_s();
	var dCX_dsesm = get_dCX_dsesm(dat,S,s,X,x,i); // equation 68
	var dsesm_ds8msmsp = get_dzesm_dz8msmsp(dat,S,s,i); // equation 60
	return dCX_dsesm*dsesm_ds8msmsp;
    } else {
	return 0;
    }
}

// equation 70
get_dCX_d8smss = function(dat,X,x,i){
    if (dat.standardson()){
	var S = dat.get_S();
	var s = dat.get_s();
	var dCX_dsesm = get_dCX_dsesm(dat,S,s,X,x,i); // equation 68
	var dsesm_d8smss = get_dzesm_d8zmss(dat,S,s,i); // equation 61
	return dCX_dsesm*dsesm_d8smss;
    } else {
	return 0;
    }
}

// equation 71
get_dCX_dyxmsmsp = function(dat,X,x,i){
    var dCX_dxesm = get_dCX_dxesm(dat,X,x,i); // equation 72
    var y = dat.get_y(X);
    var dxesm_dyxmsmsp = get_dxesm_dyxmsmsp(dat,X,x,y,i); // equation 14
    return dCX_dxesm*dxesm_dyxmsmsp;
}

// equation 72
get_dCX_dxesm = function(dat,X,x,i){
    if (dat.standardson()){
	var S = dat.get_S();
	var s = dat.get_s();
	var sesm = get_zesm(dat,S,s,i); // equation 12
	return (MM('U')*A(s)*dat.CS(S))/(MM(S)*A(238)*sesm);
    } else {
	var mass = dat.get_mass(i)[0];
	var out = (mass > 0) ? MM(X)/(A(x)*mass) : 0 ;
	return out;
    }
}

// equation 73
get_deU_du8msmsp = function(dat,i){
    var dCU_du8msmsp =  get_dCU_du8msmsp(dat,i); // equation 62
    var dCTh_du8msmsp = get_dCX_du8msmsp(dat,'Th',232,i); // equation 67
    var dCSm_du8msmsp = (dat.inSpike('Sm')) ? get_dCX_du8msmsp(dat,'Sm',147,i) : 0 ; // equation 67
    return dCU_du8msmsp + 0.235*dCTh_du8msmsp + 0.00121*dCSm_du8msmsp;
}

// equation 74
get_deU_ds8msmsp = function(dat,i){
    var dCU_ds8msmsp =  get_dCU_ds8msmsp(dat,i); // equation 64
    var dCTh_ds8msmsp = get_dCX_ds8msmsp(dat,'Th',232,i); // equation 69
    var dCSm_ds8msmsp = (dat.inSpike('Sm')) ? get_dCX_ds8msmsp(dat,'Sm',147,i) : 0 ; // equation 69
    return dCU_ds8msmsp + 0.235*dCTh_ds8msmsp + 0.00121*dCSm_ds8msmsp;
}

// equation 75
get_deU_d8smss = function(dat,i){
    var dCU_d8smss =  get_dCU_d8smss(dat,i); // equation 66
    var dCTh_d8smss = get_dCX_d8smss(dat,'Th',232,i); // equation 70
    var dCSm_d8smss = (dat.inSpike('Sm')) ? get_dCX_d8smss(dat,'Sm',147,i) : 0 ; // equation 70
    return dCU_d8smss + 0.235*dCTh_d8smss + 0.00121*dCSm_d8smss;
}

// function 76
get_deU_dv2msmsp = function(dat,i){
    var dCTh_dv2msmsp = get_dCX_dyxmsmsp(dat,'Th',232,i); // equation 71
    return 0.235*dCTh_dv2msmsp;
}

// equation 77
get_deU_dw7msmsp = function(dat,i){
    var deU_dw7msmsp = dat.inSpike('Sm') ? get_dCX_dyxmsmsp(dat,'Sm',147,i) : 0 ; // equation 71
    return 0.00121*deU_dw7msmsp;
}

// equation 78
get_deU_d7smsmsp = function(dat,i){
    var dSm_d7smsmsp = (dat.inStand('Sm') & !dat.inSpike('Sm')) ? get_dCX_dxsmsm(dat,'Sm',147,i) : 0; // equation 32
    return 0.00121*dSm_d7smsmsp;
}

// equation 79
get_deU_ds7mss = function(dat,i){
    var dCSm_ds7mss = (dat.inStand('Sm') & !dat.inSpike('Sm')) ? get_dCX_dsxmss(dat,'Sm',147,i) : 0; // equation 33
    return 0.00121*dCSm_ds7mss;
}

// equation 80
get_dD_dt = function(dat,t,Ft,i){
    var dX8_dt = 8*lambda(238)*Math.exp(lambda(238)*t)*137.818/138.818; // equation 81
    var dX5_dt = 7*lambda(235)*Math.exp(lambda(235)*t)/138.818; // equation 82
    var dX2_dt = 6*lambda(232)*Math.exp(lambda(232)*t); // equation 83
    var dX7_dt = 0.1499*lambda(147)*Math.exp(lambda(147)*t); // equation 84
    var U8esm = get_xesm(dat,'U',238,i); // equation 8
    var T2esm = get_xesm(dat,'Th',232,i); // equation 8
    var S7esm = 0;
    if (dat.inSpike('Sm')){
	S7esm = get_xesm(dat,'Sm',147,i); // equation 8
    } else if (dat.inStand('Sm')){
	S7esm = get_zesm(dat,'Sm',147,i); // equation 12
    }
    return (dX8_dt*Ft[3]+dX5_dt*Ft[2])*U8esm + dX2_dt*Ft[1]*T2esm + dX7_dt*Ft[0]*S7esm;
}

// equation 85
get_dD_du8msmsp = function(dat,t,Ft,i){
    var u = dat.get_y('U');
    var d8esm_du8msmsp = get_dxesm_dyxmsmsp(dat,'U',238,u,i); // equation 14
    var d7esm_du8msmsp = (dat.inStand('Sm') & !dat.inSpike('Sm')) ? get_dzesm_du8msmsp(dat,'Sm',147,i) : 0 ; // equation 58
    return (X8(t)*Ft[3]+X5(t)*Ft[2])*d8esm_du8msmsp + X7(t)*Ft[0]*d7esm_du8msmsp;
}

// equation 86
get_dD_d78msmsp = function(dat,t,Ft,i){
    var d7esm_d78msmsp = (dat.inStand('Sm') & !dat.inSpike('Sm')) ? get_dzesm_dz8msmsp(dat,'Sm',147,i) : 0 ; // equation 60
    return X7(t)*Ft[0]*d7esm_d78msmsp;
}

// equation 87
get_dD_d87mss = function(dat,t,Ft,i){
    var d7esm_d87mss = (dat.inStand('Sm') & !dat.inSpike('Sm')) ? get_dzesm_d8zmss(dat,'Sm',147,i) : 0 ; // equation 61
    return X7(t)*Ft[0]*d7esm_d87mss;
}

// equation 88
get_dD_dv2msmsp = function(dat,t,Ft,i){
    var v = dat.get_y('Th');
    var d2esm_dv2msmsp = get_dxesm_dyxmsmsp(dat,'Th',232,v,i); // equation 14
    return X2(t)*Ft[1]*d2esm_dv2msmsp;
}

// equation 89
get_dD_dw7msmsp = function(dat,t,Ft,i){
    var d7esm_dw7msmsp = (dat.inSpike('Sm')) ? get_dxesm_dyxmsmsp(dat,'Sm',147,dat.get_y('Sm'),i) : 0 ; // equation 14
    return X7(t)*Ft[0]*d7esm_dw7msmsp;
}
