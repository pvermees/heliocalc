$(function(){

    //////////// function definitions //////////////////
    
    // refresh the screen based on the current state of mydata
    function refresh(){
	setSettings();
	setStdTabs();
	setSpkTabs();
	showorhide();
	var smpheaders = get_smpheaders();
	var resheaders = get_resheaders();
	setTable(smpheaders,mydata.json2handson(),'#smptable');
	setTable(resheaders,crunch(),'#restable'); // crunch does all the actual work
	localStorage.setItem('data',JSON.stringify(mydata.data));
    }
    
    function crunch(){
	var out = [];
	var names = mydata.names();
	var issmp, sname, UsU, ThsTh, SmsSm, HesHe;
	var t, Ft, tstar, mass, Uppm, Thppm, Smppm, eU;
	var covmat_age_eU, err_t, err_tstar, err_eU;
	var smp;
	for (var i=0; i<names.length; i++){ // loop through all the ICPMS data
	    issmp = !mydata.standardson() | (names[i].indexOf(mydata.prefix()) < 0); // smp or std?
	    if (issmp){
		// compute output values
		sname = names[i];
		UsU = mydata.get_moles('U',i);
		ThsTh= mydata.get_moles('Th',i);
		SmsSm = mydata.get_moles('Sm',i);
		HesHe = mydata.get_moles('He',i);
		t = get_age(UsU[0],ThsTh[0],SmsSm[0],HesHe[0],[1,1,1,1],i);
		Ft = mydata.get_Ft(i);
		tstar = get_age(UsU[0],ThsTh[0],SmsSm[0],HesHe[0],Ft,i);
		mass = mydata.get_mass(i);
		Uppm = mydata.get_conc('U',i);
		Thppm = mydata.get_conc('Th',i);
		eU = Uppm[0] + 0.235*Thppm[0];
		if (mydata.doSm()){
		    Smppm = mydata.get_conc('Sm',i);
		    eU += 0.00121*Smppm[0];
		}
		covmat_age_eU = get_covmat_age_eU(mydata,t,[1,1,1,1],i); // uncorrected
		err_t = Math.sqrt(covmat_age_eU[0][0]);
		covmat_age_eU = get_covmat_age_eU(mydata,tstar,Ft,i); // corrected
		err_tstar = Math.sqrt(covmat_age_eU[0][0]);
		err_eU = Math.sqrt(covmat_age_eU[1][1]);
		// add output values to array
		smp = [sname];
		smp = append(smp,UsU[0],UsU[1]);
		smp = append(smp,Ft[2],null);
		smp = append(smp,Ft[3],null);
		smp = append(smp,ThsTh[0],ThsTh[1]);
		smp = append(smp,Ft[1],null);
		if (mydata.doSm()){
		    smp = append(smp,SmsSm[0],SmsSm[1]);
		    smp = append(smp,Ft[0],null);
		}
		smp = append(smp,HesHe[0],HesHe[1]);
		smp = append(smp,t,err_t);
		smp = append(smp,tstar,err_tstar);
		if (mydata.standardson()) { smp = append(smp,mass[0],mass[1]); }
		else { smp = append(smp,mass[0],null); }
		smp = append(smp,Uppm[0],Uppm[1]);
		smp = append(smp,Thppm[0],Thppm[1]);
		if (mydata.doSm()){ smp = append(smp,Smppm[0],Smppm[1]); }
		smp = append(smp,eU,err_eU);
		smp = append(smp,cov2cor(covmat_age_eU),null);
		if (mydata.standardson()){ smp = mydata.appendConcentrations(smp,i); }
		out.push(smp);
	    }
	}
	return out;
    }
    
    function showorhide(){
	if (mydata.standardson()){
	    $("#standardspec").show();
	    $("#prefixbox").show();
	} else {
	    $("#standardspec").hide();
	    $("#prefixbox").hide();
	}
    }
    
    function get_resheaders(){ // results
	var out = ["Name","U [pmol]","&sigma;(U)","F<sub>t</sub>(<sup>235</sup>U)"];
	out.push("F<sub>t</sub>(<sup>238</sup>U)");
	out.push("Th [pmol]","&sigma;(Th)","F<sub>t</sub>(<sup>232</sup>Th)");
	if (mydata.doSm()){ out.push("Sm [pmol]","&sigma;(Sm)","F<sub>t</sub>(<sup>147</sup>Sm)"); }
	out.push("He [pmol]","&sigma;(He)");
	out.push("t[Ma]","&sigma;(t)","t*[Ma]","&sigma;(t*)");
	out.push("mass [&mu;g]");
	if (mydata.standardson()) { out.push("&sigma;(mass)") }
	out.push("U [ppm]","&sigma;(U)","Th [ppm]","&sigma;(Th)");
	if (mydata.doSm()) {out.push("Sm [ppm]","&sigma;(Sm)");}
	out.push("eU","&sigma;(eU)","&rho;(t*,eU)");
	if (mydata.standardson()){ out = mydata.appendConcHeaders(out);	}
	return out;
    }
    
    function get_smpheaders(){
	var out = ["Name","spike [ml]"];
	out.push("length [&mu;m]","width [&mu;m]","height [&mu;m]");
	out.push("<sup>4</sup>He [ncc]","&sigma;(<sup>4</sup>He)");
	out.push("<sup>238</sup>U","&sigma;(<sup>238</sup>U)");
	out.push("<sup>232</sup>Th","&sigma;(<sup>232</sup>Th)");
	if (mydata.inSpike('Sm')){ out.push("<sup>147</sup>Th","&sigma;(<sup>147</sup>Th)"); }
	if (mydata.standardson()){ out = addHeaders(out,mydata.nuclideList('std',true)); }
	out = addHeaders(out,mydata.nuclideList('spk'));
	return out;
    }

    function input2json(){
	mydata.settings2json();
	mydata.stdTabs2json();
	mydata.spkTabs2json();
	mydata.handson2json();
    }
    
    function save2localStorage(){
	input2json();
	mydata.save2localStorage();
    }
    
    function export2helioplot(fname,doSm){
	var res = $("#restable").handsontable("getData");
	var U,sU,Th,sTh,Sm,sSm,C,He,sHe;
	var out=fname.replace(/.csv/,'')+'\n';
	for (var i=0; i<res.length; i++){
	    U = res[i][1]*res[i][4];
	    sU = res[i][2]*res[i][4];
	    Th = res[i][5]*res[i][7];
	    sTh = res[i][6]*res[i][7];
	    if (doSm) {
		Sm = res[i][8]*res[i][10];
		sSm = res[i][9]*res[i][10];
		He = res[i][11];
		sHe = res[i][12];
	    } else {
		Sm = '';
		sSm = '';
		He = res[i][8];
		sHe = res[i][9];
	    }
	    C = Sm;
	    if (U>0 & He>0){
		out += [U,sU,Th,sTh,Sm,sSm,He,sHe,C].join(',');
		out += '\n';
	    }
	}
	return encodeURIComponent(out);
    }

    function export2radialplotter(fname,doSm){
	var tstar, err_tstar, issmp, UsU, ThsTh, SmsSm, HesHe, tstar;
	var out=fname.replace(/.csv/,'')+',O\n';
	var names = mydata.names();
	for (var i=0; i<names.length; i++){ // loop through all the ICPMS data
	    issmp = !mydata.standardson() | (names[i].indexOf(mydata.prefix()) < 0); // smp or std?
	    if (issmp){
		UsU = mydata.get_moles('U',i);
		ThsTh= mydata.get_moles('Th',i);
		SmsSm = mydata.get_moles('Sm',i);
		HesHe = mydata.get_moles('He',i);
		Ft = mydata.get_Ft(i);
		tstar = get_age(UsU[0],ThsTh[0],SmsSm[0],HesHe[0],Ft,i);
		err_tstar = get_err_age_radialplotter(mydata,tstar,Ft,i);	
		if (tstar>0){
		    out += tstar;
		    out += ',';
		    out += err_tstar;
		    out += '\n';
		}
	    }
	}
	return encodeURIComponent(out);
    }
    
    function export2csv(){
	var header = $("#restable").handsontable('getColHeader');
	var dat = $("#restable").handsontable('getData');
	var out = header.join(',');
	out = out.replace(/&sigma;/g,'s.e.');
	out = out.replace(/&mu;/g,'u');
	out = out.replace(/<sup>|<\/sup>|<sub>|<\/sub>/g,'');
	for (var i=0; i<dat.length; i++){
	    out += '\n';
	    out += dat[i].join(',');
	}
	return encodeURIComponent(out);
    }
        
    function addHeaders(headers,nl){
	var label = "";
	for (var i=0 ; i < nl.length ; i++){
	    label = "<sup>" + nl[i][0] + "</sup>" + nl[i][1];
	    headers.push(label);
	    headers.push("&sigma;(" + label + ")");
	}
	return headers;
    }
    
    function setSettings(){
	$('#onoff').val(mydata.standardsonoff()).selectmenu('refresh');
	$('#mineral').val(mydata.mineral()).selectmenu('refresh');
	$('#habit').val(mydata.habit()).selectmenu('refresh');
	$('#prefix').val(mydata.prefix());
    }

    function setStdTabs(){
	setTabs(false);
    }

    function setSpkTabs(){
	setTabs(true);
    }

    function setTabs(spk){
	var theclass = spk ? '.spk' : '.std';
	var json = spk ? mydata.data.spk : mydata.data.std;
	var selector1, selector2, selector3, selector4;
	$.each(json, function(k, v) {
	    selector1 = 'input[type=text][element=' + k + ']' + theclass;
	    if (spk) {
		selector2 = 'input[type=text][element=' + k + '].ratio'; 
		selector3 = 'sup[element=' + k + '].numerator';
	    }
	    if (v.isotope != null){
		$(selector1).val(v.conc);
		if (spk) {
		    $(selector2).val(v.ratio);  
		    $(selector3).text(v.isotope);
		}
		selector4 = 'input[id=' + v.isotope + ']';
	    } else {
		$(selector1).val('');
		if (spk) {
		    $(selector2).val('');
		    $(selector3).text('x');
		}
		selector4 = 'input' + theclass + '[element=' + k + '].NA';
	    }
	    $(selector4).prop('checked',true);
	});
    }

    function setTable(headers,data,theID){
	$(theID).handsontable({
	    colHeaders: headers,
	    data: data,
	    contextMenu: true,
	    observeChanges: true,
	    manualColumnResize: true,
	});
    }
    
    function loadDefaults(){
	if ($("#onoff").val() === "on"){ readJSON('./data/standardson.json'); } 
	else { readJSON('./data/standardsoff.json'); }
    }
    
    function run(){
	save2localStorage();
	refresh();
    }
    
    $("#onoff").selectmenu({ change : function(event,ui){
	var isoff = $("#onoff").val() === "off";
	var wason = mydata.standardson();
	if ((isoff & wason) | (!isoff & !wason)){
            loadDefaults();
	    run();
	}
    }});
    
    $(".setting").selectmenu({ width : 'auto'});

    $(".isotope").click(function(){ // click on selected isotopes
	var v = null, myclass = null;
	if ($(this).hasClass('std')){ myclass = 'std'; }
	if ($(this).hasClass('spk')){ myclass = 'spk'; }
	if (!$(this).hasClass('NA')){ v = this.id; }
	var element = $(this).attr('element');
	mydata.data[myclass][element]['isotope'] = v;
	if (!mydata.data.smp.ICPdat[element].hasOwnProperty(v)){ // no data available for isotope
	    mydata.data.smp.ICPdat[element][v] = {signal:{},err:{}}; // create JSON object for new isotope
	}
	refresh();
    });
    
    $("#DEFAULTS").button().click(function( event ) {
	loadDefaults();
    });
    
    $("#CLEAR").button().click(function( event ) {
	localStorage.clear();
	$('.handsontable').each(function(){
	    $(this).handsontable('clear');
	});
    });
    
    $("#RUN").button().click(function( event ) {
	run();
    });
    
    $("#OPEN").on('change', function(e){
	var file = e.target.files[0];
	var reader = new FileReader();
	reader.onload = function(e){
	    mydata.data = JSON.parse(this.result);
	    refresh();
	}
	reader.readAsText(file);
    });

    $("#IMPORT").on('change', function(e){
	var file = e.target.files[0];
	var reader = new FileReader();
	reader.onload = function(e){
	    mydata.data.smp = mydata.csv2json(this.result);
	    mydata.data.std = mydata.crosscheck(mydata.data.std,mydata.data.smp.ICPdat);
	    setTable(get_smpheaders(),mydata.json2handson(),'#smptable');
	    setStdTabs();
	    localStorage.setItem('data',JSON.stringify(mydata.data));
	}
	reader.readAsText(file);
    });
    
    $(".conc").change(function(e){
	mydata.stdTabs2json();
	mydata.spkTabs2json();
    });
    
    $("#SAVE").button().click(function( event ) {
	var fname = prompt("Please enter a file name", "data.json");
	if (fname != null){
	    save2localStorage();
	    $('#fname').attr("href","data:text/plain," + JSON.stringify(mydata.data));
	    $('#fname').attr("download",fname);
	    $('#fname')[0].click();
	}
    });
    
    $("#EXPORT").button().click(function( event ) {
	var fname = prompt("Please enter a file name", "results.csv");
	if (fname != null){
	    $('#fname').attr("href","data:text/plain," + export2csv());
	    $('#fname').attr("download",fname);
	    $('#fname')[0].click();
	}
    });
    
    $("#HELIOPLOT").button().click(function( event ) {
	var fname = prompt("Please enter a file name", "helioplot.csv");
	if (fname != null){
	    $('#fname').attr("href","data:text/plain," + export2helioplot(fname,mydata.doSm()));
	    $('#fname').attr("download",fname);
	    $('#fname')[0].click();
	}
    });

    $("#RADIALPLOTTER").button().click(function( event ) {
	var fname = prompt("Please enter a file name", "radialplotter.csv");
	if (fname != null){
	    $('#fname').attr("href","data:text/plain," + export2radialplotter(fname,mydata.doSm()));
	    $('#fname').attr("download",fname);
	    $('#fname')[0].click();
	}
    });
    
    $("#HELP").button().click(function( event ) {
	if ($("#HELP").text() === 'HELP'){
	    $("#restable").hide();
	    $("#outputlabel").text('Help:');
	    $("#documentation").show();
	    $("#documentation").load('help.html');
	    $("#HELP").text('HIDE');
	} else {
	    $("#documentation").hide();
	    $("#restable").show();
	    $("#outputlabel").text('Output:');
	    $("#HELP").text('HELP');
	}
	$("#HELP").css({ 'padding-left': '25px', 'padding-right': '25px', 
			 'padding-top': '10px', 'padding-bottom': '10px' });
    });
    
    //////////// define and format html items //////////////////
    
    $("#standardspec").tabs();
    $("#spikespec").tabs();
    $("#documentation").width(1200);
    $(".tabbed").width(600);
    $(".conc").width(70);
    $(".ratio").width(70);
    $("#prefix").width(50);
    $("#smptable").handsontable({
	colHeaders: [],
	minCols: 1,
	minRows: 500,
	manualColumnResize: true,
	stretchH: 'all',
	height: 257,
	contextMenu: true
    });
    $("#restable").handsontable({
	colHeaders: [],
	minCols: 1,
	minRows: 500,
	manualColumnResize: true,
	stretchH: 'all',
	height: 257,
	contextMenu: false
    });
    
    //////////// populate with data //////////////////
    
    function readJSON(json){
	$.getJSON(json, function(data){
	    localStorage.clear();
	    mydata = new heliodata(data);
	    refresh();
	});
    }
    
    var data = JSON.parse(localStorage.getItem('data'));
    mydata = new heliodata(data);
    if (mydata.data === null){
	loadDefaults();
    } else {
	refresh();
    }

    $("#smptable").resizable();
    $("#restable").resizable();
    
});

//////////// save when leaving or refreshing //////////////////

$(window).bind('unload', function(){
    run();
});

$(window).bind('reload', function(){
    run();
});
