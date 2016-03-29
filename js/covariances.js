EA = function(dat,Z,z,i){
    var a = get_a(dat,i);
    var b = get_b(dat,Z,z,i);
    var c = get_c(dat,Z,z,i);
    var d = get_d(dat,Z,z,i);
    var out = [[a,b,0],
	       [b,c,0],
	       [0,0,d]];
    return out;
}

EB = function(dat,X,x,i){
    var S = dat.get_S();
    var s = dat.get_s();
    var err_xsmsm = Math.pow(get_err_yxmsmsp(dat,X,x,S,s,i),2);
    var err_sxmss = get_zxmss(dat,S,s,X,x)[1];
    var out = [[err_xsmsm,0],
	       [0,err_sxmss]];
    return out;
}

ED = function(dat,X,x,i){
    var S = dat.get_S();
    var s = dat.get_s();
    var a = get_a(dat,i);
    var b = get_b(dat,S,s,i);
    var c = get_c(dat,S,s,i);
    var d = get_d(dat,S,s,i);
    var y = dat.get_y(X);
    var s2yxmsmsp = get_var_yxmsmsp(dat,X,y,X,x,i);
    var out = [[a,b,0,0],
	       [b,c,0,0],
	       [0,0,d,0],
	       [0,0,0,s2yxmsmsp]];
    return out;
}

EG = function(dat,i){
    var S, s;
    var a = get_a(dat,i);
    var bs = 0;
    var cs = 0;
    var ds = 0;
    if (dat.standardson()){
	S = dat.get_S();
	s = dat.get_s();
	bs = get_b(dat,S,s,i);
	cs = get_c(dat,S,s,i);
	ds = get_d(dat,S,s,i);
    }
    var b7 = 0;
    var c7 = 0;
    var d7 = 0;
    if (dat.doSm()){
	b7 = get_b(dat,'Sm',147,i);
	c7 = get_c(dat,'Sm',147,i);
	d7 = get_d(dat,'Sm',147,i);
    }
    var e = get_e(dat,i);
    var f = get_f(dat,i);
    var g = 0;
    var h = 0;
    var ii = 0;
    var j = 0;
    var k = 0;
    if (dat.standardson()){
	g = get_g(dat,S,s,i);
	h = get_h(dat,S,s,i);
	ii = get_i(dat,S,s,i);
	j = get_j(dat,S,s,i);
	k = get_k(dat,i);
    }
    var l = get_l(dat,i);
    var out = [[a ,bs,0 ,0,0,0 ,0,b7,0 ,0],
	       [bs,cs,0 ,0,0,ii,0,0 ,0 ,0],
	       [0 ,0 ,ds,0,0,0 ,k,0 ,0 ,0],
	       [0 ,0 ,0 ,e,0,0 ,0,0 ,0 ,0],
	       [0 ,0 ,0 ,0,f,h ,0,0 ,0 ,0],
	       [0 ,ii,0 ,0,h,g ,0,0 ,0 ,0],
	       [0 ,0 ,k ,0,0,0 ,j,0 ,0 ,0],
	       [b7,0 ,0 ,0,0,0 ,0,c7,0 ,0],
	       [0 ,0 ,0 ,0,0,0 ,0,0 ,d7,0],
	       [0 ,0 ,0 ,0,0,0 ,0,0 ,0 ,l]];
    return out;
}

EH = function(dat,i){
    var S, s;
    var a = get_a(dat,i);
    if (dat.standardson()){
	S = dat.get_S();
	s = dat.get_s();
    }
    var b7 = 0;
    var c7 = 0;
    if (dat.doSm()){
	b7 = get_b(dat,'Sm',147,i);
	c7 = get_c(dat,'Sm',147,i);
    }
    var e = get_e(dat,i);
    var f = get_f(dat,i);
    var l = get_l(dat,i);
    var out = [[a ,0,0,b7,0],
	       [0 ,e,0,0 ,0],
	       [0 ,0,f,0 ,0],
	       [b7,0,0,c7,0],
	       [0 ,0,0,0 ,l]];
    return out;
}

// equation 18
function get_a(dat,i){
    var u = dat.get_y('U');
    return get_var_yxmsmsp(dat,'U',u,'U',238,i); // equation 2
}

// equation 19
function get_b(dat,Z,z,i){
    var u = dat.get_y('U');
    var Amsmsp = dat.get_signal('U',u,i);
    var Bmsmsp = dat.get_signal(Z,z,i);
    var Cmsmsp = dat.get_signal('U',238,i);
    var err_Cmsmsp = dat.get_error('U',238,i);
    return Amsmsp*Bmsmsp*Math.pow(err_Cmsmsp,2)/Math.pow(Cmsmsp,4);
}

// equation 20
function get_c(dat,Z,z,i){
    return get_var_yxmsmsp(dat,Z,z,'U',238,i); // equation 2
}

// equation 21
function get_d(dat,Z,z,i){
    return get_zxmss(dat,'U',238,Z,z,i)[1];
}

// equation 41
function get_e(dat,i){
    var v = dat.get_y('Th');
    return get_var_yxmsmsp(dat,'Th',v,'Th',232,i); // equation 2
}

// equation 42
function get_f(dat,i){
    if (dat.inSpike('Sm')){
	var w = dat.get_y('Sm');
	return get_var_yxmsmsp(dat,'Sm',w,'Sm',147,i); // equation 2
    } else {
	return 0;
    }
}

// equation 43
function get_g(dat,S,s,i){
    return get_var_yxmsmsp(dat,'Sm',147,S,s,i); // equation 2
}

// equation 44
function get_h(dat,S,s,i){
    if (dat.inSpike('Sm')){
	var A = dat.get_signal('Sm',147,i);
	var A_err = dat.get_error('Sm',147,i);
	var B = dat.get_signal(S,s,i);
	var w = dat.get_y('Sm');
	var C = dat.get_signal(S,w,i);
	return get_cov(A,A_err,B,C); // equation 49
    } else {
	return 0;
    }
}

// equation 45
function get_i(dat,S,s,i){
    var A = dat.get_signal(S,s,i);
    var A_err = dat.get_error(S,s,i);
    var B = dat.get_signal('U',238,i);
    var C = dat.get_signal('Sm',147,i);
    return get_cov(A,A_err,B,C); // equation 49
}

// equation 46
function get_j(dat,S,s,i){
    return get_zxmss(dat,S,s,'Sm',147)[1];
}

// equation 47
function get_k(dat,i){
    return get_cov_s7mss_8smss(dat,i);
}

// equation 54
function get_l(dat,i){
    var HesHe = get_moles(dat,'He',i);
    return Math.pow(HesHe[1],2);
}

// equation 48
function get_cov(A,A_err,B,C){
    var log_cov = Math.log(C)+2*Math.log(A_err)-Math.log(B)-2*Math.log(A);
    return -Math.exp(log_cov);
}
