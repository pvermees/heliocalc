var heliodata = function(data){

    this.data = data;

    this.names = function(){ return this.data.smp.sname; }

    this.mineral = function(){ return this.data.settings.mineral; }

    this.dens = function(){ return density(this.mineral()); }

    this.prefix = function(){ return this.data.prefix; }

    this.standardson = function(){ return this.standardsonoff() === 'on'; }

    this.standardsonoff = function(){ return this.data.settings.standardson; }

    this.gets= function(S){ return this.data.std[S].isotope; }

    this.getSpikeIsotope = function(element){ return this.data.spk[element].isotope; }

    this.R = function(){ return getR(this.mineral()); } 

    this.habit = function(){ return this.data.settings.habit; }

    this.length = function(i){ 
	var l = this.data.smp.length[i];
	if ($.isNumeric(l)) {
	    return Number(l);
	} else {
	    return 0;
	}
    }

    this.width = function(i){ 
	var w = this.data.smp.width[i];
	if ($.isNumeric(w)) {
	    return Number(w);
	} else {
	    return 0;
	}
    }

    this.height = function(i){ 
	var h = this.data.smp.height[i];
	if ($.isNumeric(h)) {
	    return Number(h);
	} else {
	    return 0;
	}
    }

    // get the internal standard
    this.getS = function(){
	switch (this.mineral()) {
	case "apatite": return 'Ca';
	case "zircon": return 'Si';
	default: return 'Ca';
	}
    }

    this.doSm = function(){
	var c1 = (this.standardson() & this.data.std.Sm.isotope !== null);
	var c2 = (this.data.spk.Sm.isotope !== null);
	return (c1 | c2) ;
    }

    this.inSpike = function(element){ return this.presentin('spk',element); }

    this.inStand = function(element){
	if (element === 'U') { return true; }
	return this.presentin('std',element);
    }

    this.isStand = function(i){
	return (this.data.smp.sname[i].indexOf(this.prefix()) >= 0);
    }

    this.presentin = function(spkstd,element){
	var out = false;
	$.each(this.data[spkstd], function(k, v) {
	    if (k === element & v.isotope !== null) { out = true; }
	});
	return(out);
    }

    this.moles = function(X,i){
	var m = 0, err_m = 0;
	if (X === 'He'){
	    m = this.signal('He',4,i) / 22.414;
	    err_m = this.error('He',4,i) / 22.414;
	} else {
	    var x = getNaturalIsotope(X);
	    var S = this.getS();
	    var s = this.gets(S);
	    var xesm = 0, err_xesm = 0;
	    var I;
	    if (this.inSpike(X)){
		I = this.info(X,x,S,s,i)
		xesm = I.get_xesm();
		err_xesm = I.get_errxesm();
	    } else if (this.inStand(X) & this.inSpike('U')){
		I = this.info('U',238,X,x,i);
		xesm = I.get_sesm();
		err_xesm = I.get_errsesm();
	    }
	    m = xesm/A(x);
	    err_m = err_xesm/A(x);
	}
	return [m,err_m];
    }

    this.signal = function(element,isotope,i){
	if (element==='He'){
	    return this.data.smp.HEdat[element][isotope].signal[i];
	} else {
	    return this.data.smp.ICPdat[element][isotope].signal[i];
	}
    }

    this.error = function(element,isotope,i){
	if (element==='He'){
	    return this.data.smp.HEdat[element][isotope].err[i];
	} else {
	    return this.data.smp.ICPdat[element][isotope].err[i];
	}
    }

    // returns and object with all the measured and true ratios
    this.info = function(X,x,S,s,i){
	var I = new info(this,X,x,S,s,i);
	return I;
    }

    // equation 21
    this.ratio_err = function(ne,ni,de,di,i){
	var n = this.signal(ne,ni,i);
	var d = this.signal(de,di,i);
	var err_n = this.error(ne,ni,i);
	var err_d = this.error(de,di,i);
	var err = (n/d)*Math.sqrt(Math.pow(err_n/n,2) + Math.pow(err_d/d,2));
	return err;
    }

    // returns mean and variance of standard solution ratios
    this.avgStdRatio = function(ne,ni,de,di){
	var snames = this.names();
	var sum = 0, avg = 0, n = 0, ndmss = 0, mu;
	var zvz = [], wu = [];
	var z = 0, vz = 0, variance = 0;
	var zeta = 0, mu = 0, oldzeta = 0, oldmu = 0;
	for (var i = 0; i<snames.length; i++){
	    if (this.isStand(i)) {
		z = Math.log(this.signal(ne,ni,i)) - Math.log(this.signal(de,di,i));
		vz = Math.pow(this.error(ne,ni,i)/this.signal(ne,ni,i),2) + 
		    Math.pow(this.error(de,di,i)/this.signal(de,di,i),2);
		zvz.push([z,vz]);
	    }
	}
	if (zvz.length < 2){ // when there is only one standard
	    avg = Math.exp(z);
	    variance = avg*avg*vz;
	    return [avg,variance]; 
	}
	for (var j = 0; j<10; j++){
	    mu = getMu(zvz,zeta);
	    zeta = getZeta(zvz,mu);
	    if (zeta === null) {
		mu = getMu(zvz,0);
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
	    var x = getwuzumu(zvz,zeta,mu);
	    variance = avg*avg/x.sumwu;
	}
	return [avg,variance]; 
    }

    // get array of alpha-retention factors and uncertainties for ith grain
    this.getFt = function(i){
	var l = this.length(i);
	var w = this.width(i);
	var h = this.height(i);
	var R = this.R();
	return FT(R,l,w,h,this.habit());
    }

    // get the volume of the ith grain, in um^3
    this.getV = function(i){
	var V = 0, err_V = 0;
	if (this.standardson()){
	    var S = this.getS();
	    var s = this.gets(S);
	    var IU = this.info('U',238,S,s,i);
	    V = IU.v*IU.get_sesm()/A(IU.s);
	    err_V = IU.v*IU.get_errsesm()/A(IU.s);
	} else {
	    var l = this.length(i);
	    var w = this.width(i);
	    var h = this.height(i);
	    V = VOL(l,w,h,this.habit());
	}
	return [V, err_V];
    }

    // get the mass of the ith grain, in ug
    this.getMass = function(i){
	var V = this.getV(i); // in um^3
	var m = V[0] * this.dens() / 1e6; // mass in ug
	var err_m = V[1] * this.dens() / 1e6; // mass in ug
	return [m, err_m];
    }

    this.getC = function(X,x,i){
	var CX = 0, err_CX = 0;
	var S = this.getS();
	var s = this.gets(S);
	var spiked = this.inSpike(X);
	var J = [0,0];
	var E = [[0,0],[0,0]];
	if (!this.standardson()) {
	    var I = this.info(X,x,S,s,i);
	    var xesm = I.get_xesm();
	    var mass = this.getMass(i);
	    if (mass[0] > 0) {
		var MmAx = MM(X)/(mass[0]*A(I.x));
		var CX = xesm*MmAx;
		var err_CX = I.get_errxesm()*MmAx;
	    }
	    return [CX,err_CX];
	}
	if (this.inSpike(X) & this.inStand(X)){        // equation 7
	    var I = this.info(X,x,S,s,i);
	    var xesm = I.get_xesm();
	    var sesm = I.get_sesm();
	    CX = getCSMXAsMSAx(X,x,S,s)*xesm/sesm;
	    J = I.JB();
	    E = I.SigmaB();                   // equation 22
	} else if (!this.inSpike() & this.inStand(X)){ // equation 2
	    var I = this.info(X,x,S,s,i);
	    CX = getCSMXAsMSAx(X,x,S,s)*I.xsmsmsp*I.xstss*I.sxmss[0];
	    J = I.JA();
	    E = I.SigmaA();                     // equation 18
	} else if (this.inSpike(X) & this.inStand('U')){ // bespoke method for Th
	    var IX = this.info(X,x,S,s,i);
	    var IU = this.info('U',238,S,s,i);
	    var xesm = IX.get_xesm();
	    var sesm = IU.get_sesm();
	    CX = getCSMXAsMSAx(X,x,S,s)*xesm/sesm;
	    J = JB2(IX,IU);
	    E = SigmaB2(IX,IU);
	}
	err_CX = propagate(J,E);
	return [CX,err_CX];
    }

    this.get_eU = function(i){
	var eU = this.getC('U',238,i)[0] + 0.235 * this.getC('Th',232,i)[0];
	var J = this.JeU(i);
	var E = this.Sigma_eU(i);
	var err_eU = propagate(J,E);
	return [eU,err_eU];
    }

    this.appendConcentrations = function(smparray,i){
	var me = this;
	$.each(this.data.std, function(k, v) {
	    if ((me.mineral() === "apatite" & k!=="Ca" & k!=="P" & k!=="Sm" & v.isotope!==null) |
		(me.mineral() === "zircon" & k!=="Si" & k!=="Zr" & k!=="Hf" & k!=="Sm" & v.isotope!==null)){
		isotope = me.data.std[k].isotope;
		C = me.getC(k,isotope,i);
		smparray = append(smparray,C[0],C[1]);
	    }
	});
	return smparray;
    }

    this.appendConcHeaders = function(headarray){
	var me = this;
	$.each(this.data.std, function(k, v) {
	    if ((me.mineral() === "apatite" & k!=="Ca" & k!=="P" & k!=="Sm" & v.isotope!==null) |
		(me.mineral() === "zircon" & k!=="Si" & k!=="Zr" & k!=="Hf" & v.isotope!==null)){
		headarray.push(k+" [ppm]","&sigma;("+k+")");
	    }
	});
	return headarray;
    }

    this.JeU = function(i){
	var S = this.getS();
	var s = this.gets(S);
	var IU = this.info('U',238,S,s,i);
	var ITh = this.info('Th',232,S,s,i);
	var xesmU = IU.get_xesm();
	var xesmTh = ITh.get_xesm();
	var sesm = IU.get_sesm();
	var CU = this.getC('U',238,i)[0];
	var CTh = this.getC('Th',232,i)[0];
	var eU = CU + 0.235 * CTh;
	var deUdCU = 1;
	var deUdCTh = 0.235;
	var dCUdxesmU = IU.dCX_dxesm();
	var dCThdxesmU = -(CTh/sesm)*IU.dsesm_dxesm();
	var dCThdxesmTh = CTh/xesmTh;
	var dCThdsesm = -CTh/sesm;
	var deUdxesmU = deUdCU*dCUdxesmU + 0.235*deUdCTh*dCThdxesmU;
	var deUdyxmsmspU = deUdxesmU*IU.dxesm_dyxmsmsp();
	var deUdyxmsmspTh = deUdCTh*dCThdxesmTh*ITh.dxesm_dyxmsmsp();
	var deUdsxmsmsp = (deUdCU*IU.dCX_dsesm()+deUdCTh*dCThdsesm)*IU.dsesm_dsxmsmsp();
	var deUdxsmss = (deUdCU*IU.dCX_dsesm()+deUdCTh*dCThdsesm)*IU.dsesm_dxsmss();
	var J = [deUdyxmsmspU, deUdyxmsmspTh, deUdsxmsmsp, deUdxsmss];
	return J;
    }

    this.Sigma_eU = function(i){
	var S = this.getS();
	var s = this.gets(S);
	var IU = this.info('U',238,S,s,i);
	var ITh = this.info('Th',232,S,s,i);
	var var_yxmsmspU = Math.pow(IU.err_yxmsmsp,2);
	var var_yxmsmspTh = Math.pow(ITh.err_yxmsmsp,2);
	var var_sxmsmsp = Math.pow(IU.err_sxmsmsp,2);
	var var_xsmss = IU.var_xsmss;
	var E = [[var_yxmsmspU,0,0,0],
		 [0,var_yxmsmspTh,0,0],
		 [0,0,var_sxmsmsp,0],
		 [0,0,0,var_xsmss]];
	return E;
    }

    this.Jt = function(i,doFt){
	var S = this.getS();
	var s = this.gets(S);
	var IU = this.info('U',238,S,s,i);
	var ITh = this.info('Th',232,S,s,i);
	var ISm = this.info('Sm',147,S,s,i);
	var IUSm = this.info('U',238,'Sm',147,i);
	var U = this.moles('U',i)[0];
	var Th = this.moles('Th',i)[0];
	var Sm = this.moles('Sm',i)[0];
	var He = this.moles('He',i)[0];
	var Ft = doFt ? this.getFt(i) : [1,1,1,1];
	var t = age(U,Th,Sm,He,Ft);
	var V = this.getV(i);
	var R = this.R();
	var dUdyxmsmspU = IU.dxesm_dyxmsmsp()/A(238);
	var dSmdyxmsmspU = IUSm.dsesm_dyxmsmsp()/A(232);
	var dSmdsxmsmspU = IUSm.dsesm_dsxmsmsp()/A(147);
	var dSmdxsmssU = IUSm.dsesm_dxsmss()/A(147);
	var dDdt = (dX8dt(t)*Ft[3] + dX5dt(t)*Ft[2])*U + dX2dt(t)*Ft[1]*Th + dX7dt(t)*Ft[0]*Sm;
	var dDdyxmsmspU = (X8(t)*Ft[3] + X5(t)*Ft[2])*dUdyxmsmspU + X7(t)*Ft[0]*dSmdyxmsmspU;
	var dDdyxmsmspTh = X2(t)*Ft[1]*ITh.dxesm_dyxmsmsp();
	var dDdyxmsmspSm = X7(t)*Ft[0]*ISm.dxesm_dyxmsmsp();
	var dDdsxmsmspU = X7(t)*Ft[0]*dSmdsxmsmspU;
	var dDdxsmssU = X7(t)*Ft[0]*dSmdxsmssU;
	var dDdHe = 0;
	var h = -dDdyxmsmspU/dDdt;
	var i = -dDdyxmsmspTh/dDdt;
	var j = -dDdyxmsmspSm/dDdt;
	var k = -dDdsxmsmspU/dDdt;
	var l = -dDdxsmssU/dDdt;
	var m = -dDdHe/dDdt;
	var J = [h,i,j,k,l,m];
	return J;
    }

    this.Sigma_t = function(i){
	var S = this.getS();
	var s = this.gets(S);
	var IU = this.info('U',238,S,s,i);
	var ITh = this.info('Th',232,S,s,i);
	var ISm = this.info('Sm',147,S,s,i);
	var a = Math.pow(IU.err_yxmsmsp,2);
	var b = IU.cov_yxmsmsp_sxmsmsp; 
	var c = Math.pow(ITh.err_yxmsmsp,2);
	var d = Math.pow(ISm.err_yxmsmsp,2);
	var e = Math.pow(IU.err_sxmsmsp,2);
	var f = IU.var_xsmss;
	var g = Math.pow(this.error('He',4,i),2);
	var E = [[a,0,0,b,0,0],
		 [0,c,0,0,0,0],
		 [0,0,d,0,0,0],
		 [b,0,0,e,0,0],
		 [0,0,0,0,f,0],
		 [0,0,0,0,0,g]];
	return E;
    }

    this.stdTabs2json = function(){ this.tabs2json(false); }

    this.spkTabs2json = function(){ this.tabs2json(true); }

    this.tabs2json = function(spk){
	var theclass = spk ? '.spk' : '.std';
	var key = spk ? 'spk' : 'std';
	var selector1, selector2, selector3, conc, ratio, isotope;
	var data = this.data;
	$.each(this.data[key], function(k, v) {
	    selector1 = 'input[type=text][element=' + k + ']' + theclass;
	    if (spk) { selector2 = 'input[type=text][element=' + k + '].ratio'; }
	    selector3 = 'input[type=radio][element=' + k + ']' + theclass + ":checked";
	    conc = $(selector1).val();
	    ratio = $(selector2).val();
	    isotope = $(selector3).attr("id");
	    if ($(selector3).hasClass("NA")){
		data[key][k]['conc'] = null;
		data[key][k]['isotope'] = null;
		data[key][k]['ratio'] = null;
	    } else {
		data[key][k]['conc'] = Number(conc);
		data[key][k]['isotope'] = Number(isotope);
		data[key][k]['ratio'] = Number(ratio);
	    }
	});
    }

    this.settings2json = function(){
	this.data.settings.standardson = $("#onoff").val();
	this.data.settings.mineral = $("#mineral").val();
	this.data.settings.habit = $("#habit").val();
	this.data.prefix = $("#prefix").val();
    }

    this.handson2json = function(){
	var data = $('#smptable').handsontable('getData');
	var nlstd = this.nuclideList('std',true);
	var nlspk = this.nuclideList('spk');
	var signal, err, sname;
	var smp = this.data.smp;
	for (var i=0; i < data.length; i++){ // loop through rows
	    var j = 0;
	    sname = data[i][j++];
	    if (sname != null){
		smp.sname[i] = sname;
		smp.spkvol[i] = data[i][j++];
		smp.length[i] = data[i][j++];
		smp.width[i] = data[i][j++];
		if (!this.standardson()){
		    smp.height[i] = data[i][j++];
		} else {
		    smp.height[i] = 0;
		}
		smp.HEdat.He["4"].signal[i] = data[i][j++];
		smp.HEdat.He["4"].err[i] = data[i][j++];
		smp.ICPdat['U'][238].signal[i] = data[i][j++];
		smp.ICPdat['U'][238].err[i] = data[i][j++];
		smp.ICPdat['Th'][232].signal[i] = data[i][j++];
		smp.ICPdat['Th'][232].err[i] = data[i][j++];
		if (this.inSpike('Sm')){
		    smp.ICPdat['Sm'][147].signal[i] = data[i][j++];
		    smp.ICPdat['Sm'][147].err[i] = data[i][j++];          
		}
		if (this.standardson()){
		    for (var k=0; k < nlstd.length; k++){ // loop through all the standard nuclides
			smp.ICPdat[nlstd[k][1]][nlstd[k][0]].signal[i] = data[i][j++];
			smp.ICPdat[nlstd[k][1]][nlstd[k][0]].err[i] = data[i][j++];
		    }
		}
		for (var k=0; k < nlspk.length; k++){ // loop through all the spike nuclides
		    signal = data[i][j++];
		    err = data[i][j++];
		    smp.ICPdat[nlspk[k][1]][nlspk[k][0]].signal[i] = signal;
		    smp.ICPdat[nlspk[k][1]][nlspk[k][0]].err[i] = err;
		}
	    }
	}
    }

    this.csv2json = function(file){
	var lines = file.split('\n');
	var line1 = lines[0].split(',');
	var line2 = lines[1].split(',');
	var keys = this.getKeys(line1,line2);
	var json = this.initJSON(keys);
	for (var i=2; i<lines.length; i++){
	    line = lines[i].split(',');
	    if (line.length >= keys.length){
		for (var j=0; j<keys.length; j++){
		    json = this.setVal(json,keys[j],line[j]);
		}
	    }
	}
	return json;
    }

    this.json2handson = function(){
	var out = [], row;
	var mysmp = this.data.smp;
	var skip = true;
	var issmp;
	for (var i=0; i < mysmp.sname.length; i++){ // loop through samples
	    issmp = (mysmp.sname[i].indexOf(mydata.prefix) < 0); // smp or std?
	    row = [];
	    row.push(mysmp.sname[i]);
	    row.push(mysmp.spkvol[i]);
	    row.push(mysmp.length[i]);
	    row.push(mysmp.width[i]);
	    if (!this.standardson()){
		row.push(mysmp.height[i]);
	    }
	    row.push(this.signal('He',4,i));
	    row.push(this.error('He',4,i));
	    row.push(this.signal('U',238,i));
	    row.push(this.error('U',238,i));
	    row.push(this.signal('Th',232,i));
	    row.push(this.error('Th',232,i));
	    if (this.inSpike('Sm')){
		row.push(this.signal('Sm',147,i));
		row.push(this.error('Sm',147,i));
	    }
	    var me = this;
	    if (this.standardson()){
		$.each(this.data.std, function(k, v) { // loop through std isotopes
		    if (v.isotope !== null & !me.inSpike(k)){ // is element selected in settings?
			row.push(me.signal(k,v.isotope,i));
			row.push(me.error(k,v.isotope,i));
		    }
		});
	    }
	    $.each(this.data.spk, function(k, v) { // loop through spk isotopes
		if (v.isotope !== null){ // is element selected in settings?
		    row.push(me.signal(k,v.isotope,i));
		    row.push(me.error(k,v.isotope,i));
		}
	    });
	    out.push(row);
	}
	return out;
    }

    this.initJSON = function(keys){
	var out = {};
	for (var i=0; i<keys.length; i++){
	    out = this.setVal(out,keys[i],null);
	}
	return out;
    }

    // recursive function to set value based on key array
    this.setVal = function(json,key,val){
	var out = json;
	if (key.length == 1){
	    if (val === null){
		out[key[0]] = []; // initialise
	    } else {
		out[key[0]].push(val);
	    }
	} else {
	    if (!out[key[0]]){
		out[key[0]] = {}; // add key
	    }
	    out[key[0]] = this.setVal(out[key[0]],key.slice(1,key.length),val);
	}
	return out;
    }


    this.crosscheck = function(stdjson,ICPjson){
	var out = stdjson;
	$.each(stdjson, function(k, v) { // loop through std isotopes
	    if (k in ICPjson){
		// everything all right, do nothing
	    } else {
		out[k].isotope = null;
	    }
	});
	return out;
    }

    // get array of arrays with JSON keys from first two lines of .csv import file
    this.getKeys = function(line1,line2){
	var key1,key2, keys = [];
	for (var i=0; i<line1.length; i++){
	    key2 = null;
	    if (line1[i] === "name"){
		key1 = ["sname"];
	    } else if (line1[i] === "spike"){
		key1 = ["spkvol"];
	    } else if (line1[i] === "length"){
		key1 = ["length"];
	    } else if (line1[i] === "width"){
		key1 = ["width"];
	    } else if (line1[i] === "height"){
		key1 = ["height"];
	    } else if (line1[i] === "He" & line2[i] === "4"){
		key1 = ["HEdat",line1[i],line2[i],"signal"];
		key2 = ["HEdat",line1[i],line2[i++],"err"];
	    } else {
		key1 = ["ICPdat",line1[i],line2[i],"signal"];
		key2 = ["ICPdat",line1[i],line2[i++],"err"];
	    }
	    keys.push(key1);
	    if (key2 !== null){ keys.push(key2); }
	}
	return keys;
    }

    this.nuclideList = function(stdorspk,skipSpike){ // returns 2D Array of isotopes and elements
	var out = [], skip = false;
	json = this.data[stdorspk];
	var me = this;
	$.each(json, function(k, v) {
	    skip = (me.inSpike(k) & skipSpike);
	    if (!(v.isotope === null) & !skip){
		out.push([v.isotope,k]);
	    }
	});
	return out;
    }

    this.save2localStorage = function(){
	localStorage.setItem( 'data' , JSON.stringify(this.data) );
    }

}
