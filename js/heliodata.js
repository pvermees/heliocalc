var heliodata = function(data){

    this.data = data;

    this.names = function(){ return this.data.smp.sname; }

    this.mineral = function(){ return this.data.settings.mineral; }

    this.dens = function(){ return density(this.mineral()); }

    this.prefix = function(){ return this.data.prefix; }

    this.standardson = function(){ return this.standardsonoff() === 'on'; }

    this.standardsonoff = function(){ return this.data.settings.standardson; }

    this.get_s= function(){ return this.data.std[this.get_S()].isotope; }

    this.get_y = function(element){ return this.data.spk[element].isotope; }

    this.R = function(){ return getR(this.mineral()); } 

    this.habit = function(){ return this.data.settings.habit; }

    // get the internal standard
    this.get_S = function(){
	var S;
	switch (this.mineral()){
	case 'apatite': 
	    if (this.inStand('Ca')) { S = 'Ca'; }
	    if (this.inStand('P')) { S = 'P'; }
	    break;
	case 'zircon' :
	    if (this.inStand('Si')) { S = 'Si'; }
	    if (this.inStand('Zr')) { S = 'Zr'; }
	    break;
	case 'titanite' :
	    if (this.inStand('Ti')) { S = 'Ti'; }
	    if (this.inStand('Si')) { S = 'Si'; }
	    break;
	case 'monazite' :
	    if (this.inStand('P')) { S = 'P'; }
	    break;
	case 'xenotime' :
	    if (this.inStand('Y')) { S = 'Y'; }
	    if (this.inStand('P')) { S = 'P'; }
	    break;
	case 'rutile' :
	    if (this.inStand('Ti')) { S = 'Ti'; }
	    break;
	case 'magnetite' :
	case 'haematite' :
	case 'goethite' :
	    if (this.inStand('Fe')) { S = 'Fe'; }
	    break;
	case 'barite' : 
	    if (this.inStand('Ba')) { S = 'Ba'; }
	    if (this.inStand('S')) { S = 'S'; }
	    break;
	default: S = 'none';
	}
	if (S === 'none'){
	    this.data.settings.standardson = $("#onoff").val('off');
	    refresh();
	} else {
	    return S;
	}
    }

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

    this.get_x = function(element){
	if (element === 'U'){
	    return 238;
	}
	if (element === 'Th') {
	    return 232;
	}
	if (element === 'Sm') {
	    return 147;
	}
	if (element === 'He') {
	    return 4;
	}
	if (this.inStand(element)){
	    return this.data.std[element].isotope;
	}
	return 0;
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

    this.get_signal = function(element,isotope,i){
	if (element==='He'){
	    return Number(this.data.smp.HEdat[element][isotope].signal[i]);
	} else {
	    return Number(this.data.smp.ICPdat[element][isotope].signal[i]);
	}
    }

    this.get_error = function(element,isotope,i){
	if (element==='He'){
	    return Number(this.data.smp.HEdat[element][isotope].err[i]);
	} else {
	    return Number(this.data.smp.ICPdat[element][isotope].err[i]);
	}
    }

    this.CS = function(S){ return ppm(this.mineral(),S); }

    this.get_xesm = function(X,x,i){ return get_xesm(this,X,x,i); }

    this.get_xsesm = function(X,x,i){ return get_xsesm(this,X,x,i); }

    this.get_yxtsp = function(element){ return this.data.spk[element].ratio; }

    this.get_ytsp = function(element,i){ return this.data.smp.spkvol[i] * this.data.spk[element].conc; }

    this.get_xesm = function(X,x,i){ return get_xesm(this,X,x,i); }

    this.get_zesm = function(Z,z,i){ return get_zesm(this,Z,z,i); }

    this.get_zxtss = function(Z,z,X,x){
	var Xconc = 1, Yconc = 1;
	if (X !== 'U') { Xconc = this.data.std[X].conc; }
	if (Z !== 'U') { Zconc = this.data.std[Z].conc; }
	return (A(z)*MM(X)*Zconc)/(A(x)*MM(Z)*Xconc);
    }

    // returns mean and variance of standard solution ratios
    this.get_zxmss = function(Z,z,X,x){	return get_zxmss(this,Z,z,X,x); }

    this.get_moles = function(X,x,i){ return get_moles(this,X,x,i); }

    // get array of alpha-retention factors and uncertainties for ith grain
    this.get_Ft = function(i){
	var l = this.length(i);
	var w = this.width(i);
	var h = this.height(i);
	var R = this.R();
	return FT(R,l,w,h,this.habit());
    }

    // get the volume of the ith grain, in um^3
    this.get_vol = function(i){
    }

    // get the mass of the ith grain, in ug
    this.get_mass = function(i){
	var m = 0, err_m = 0;
	if (this.standardson()){
	    var S = this.get_S();
	    var s = this.get_s();
	    var MPM = mpm(S,this.mineral());
	    var As = A(s);
	    var sesm = get_zesm(this,S,s,i);
	    var err_sesm = get_err_zesm(this,S,s,i);
	    var Mm = MM(this.mineral());
	    m = sesm*Mm/(MPM*As*1000000);
	    err_m = err_sesm*Mm/(MPM*As*1000000);
	} else {
	    var V = VOL(this.length(i),this.width(i),this.height(i),this.habit()); // in um^3
	    m = V * this.dens() / 1e6; // mass in ug
	}
	return [m, err_m];
    }

    this.get_conc = function(X,i){ return get_conc(this,X,i); }

    this.appendConcentrations = function(smparray,i){
	var me = this;
	$.each(this.data.std, function(k, v) {
	    if ((me.mineral() === "apatite" & k!=="Ca" & k!=="P" & k!=="Sm" & v.isotope!==null) |
		(me.mineral() === "zircon" & k!=="Si" & k!=="Zr" & k!=="Hf" & k!=="Sm" & v.isotope!==null)){
		isotope = me.data.std[k].isotope;
		C = me.get_conc(k,i);
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
		smp.height[i] = data[i][j++];
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
		    json = this.set_val(json,keys[j],line[j]);
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
	    row.push(mysmp.height[i]);
	    row.push(this.get_signal('He',4,i));
	    row.push(this.get_error('He',4,i));
	    row.push(this.get_signal('U',238,i));
	    row.push(this.get_error('U',238,i));
	    row.push(this.get_signal('Th',232,i));
	    row.push(this.get_error('Th',232,i));
	    if (this.inSpike('Sm')){
		row.push(this.get_signal('Sm',147,i));
		row.push(this.get_error('Sm',147,i));
	    }
	    var me = this;
	    if (this.standardson()){
		$.each(this.data.std, function(k, v) { // loop through std isotopes
		    if (v.isotope !== null & !me.inSpike(k)){ // is element selected in settings?
			row.push(me.get_signal(k,v.isotope,i));
			row.push(me.get_error(k,v.isotope,i));
		    }
		});
	    }
	    $.each(this.data.spk, function(k, v) { // loop through spk isotopes
		if (v.isotope !== null){ // is element selected in settings?
		    row.push(me.get_signal(k,v.isotope,i));
		    row.push(me.get_error(k,v.isotope,i));
		}
	    });
	    out.push(row);
	}
	return out;
    }

    this.initJSON = function(keys){
	var out = {};
	for (var i=0; i<keys.length; i++){
	    out = this.set_val(out,keys[i],null);
	}
	return out;
    }

    // recursive function to set value based on key array
    this.set_val = function(json,key,val){
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
	    out[key[0]] = this.set_val(out[key[0]],key.slice(1,key.length),val);
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
