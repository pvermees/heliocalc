var info = function(data,X,x,S,s,i){
    this.mydata = data;
    this.X = X;
    this.x = x;
    this.S = S;
    this.s = s;
    this.cov_yxmsmsp_sxmsmsp = 0;
    this.mineral = this.mydata.mineral();
    this.l = Number(this.mydata.length(i));
    this.w = Number(this.mydata.width(i));
    this.h = Number(this.mydata.height(i));
    this.v = MM(this.mydata.mineral())/this.mydata.dens()/mpm(S); // um^3/pmol
    this.y = 0;
    this.xmsmsp = 0;
    this.ytsp = 0;
    this.yxtsp = 0;
    this.xytsp = 0;
    this.yxtsm = 0;
    this.ymsmsp = 0;
    this.yxmsmsp = 0;
    this.err_yxmsmsp = 0;
    this.err_sxmsmsp = 0;
    this.err_xmss = 0;
    this.sxmsmsp = 0;
    this.sxmsmsp = 0;
    this.xsmss = 0;
    this.var_xsmss = 0;
    this.stss = 0;
    this.xtss = 0;
    this.sxtss = 0;
    this.dVdsesm = this.v/A(this.s);
    if (this.mydata.inSpike(X) | this.mydata.inStand(X)){
	this.xmsmsp = this.mydata.signal(X,x,i);
	this.err_xmsmsp = this.mydata.error(X,x,i);
    }
    if (this.mydata.inSpike(X)){
	this.y = this.mydata.data.spk[X].isotope;
	this.ytsp = this.mydata.data.spk[X].conc*this.mydata.data.smp.spkvol[i];
	this.yxtsp = this.mydata.data.spk[X].ratio;
	this.xytsp = 1/this.yxtsp
	this.yxtsm = A(this.y)/A(x);
	this.ymsmsp = this.mydata.signal(X,this.y,i);
	this.yxmsmsp = this.ymsmsp/this.xmsmsp;
	this.err_yxmsmsp = this.mydata.ratio_err(X,this.y,X,x,i);
    }
    if (this.mydata.inStand(S)){
	this.smsmsp = this.mydata.signal(S,s,i);
	this.sxmsmsp = this.smsmsp/this.xmsmsp;
	this.xsmsmsp = 1/this.sxmsmsp;
	this.err_xsmsmsp = this.mydata.ratio_err(X,x,S,s,i);
    }
    if (this.mydata.inStand(X) & this.mydata.inStand(S)){
	var asr = this.mydata.avgStdRatio(X,x,S,s);
	this.xsmss = asr[0];
	this.var_xsmss = asr[1];
	this.err_sxmsmsp = this.mydata.ratio_err(S,s,X,x,i);
	this.stss = this.mydata.data.std[S].conc*A(s)/MM(S);
	this.xtss = A(x)/MM(X);
	if (X !== 'U') { this.xtss = this.xtss*this.mydata.data.std[X].conc; }
	this.sxtss = this.stss/this.xtss;
	this.xstss = 1/this.sxtss;
	this.sxmss = this.mydata.avgStdRatio(S,s,X,x);
    }
    if (this.mydata.inSpike(X) & this.mydata.inStand(X) & this.mydata.inStand(S)){
	this.cov_yxmsmsp_sxmsmsp = this.ymsmsp*this.smsmsp*      // equation 25
        Math.pow(this.err_xmsmsp,2)/Math.pow(this.xmsmsp,4)
    }

    this.get_xesm = function(){ // equation 8
	return (this.ytsp*this.xytsp)*(this.yxtsp-this.yxmsmsp)/(this.yxmsmsp-this.yxtsm);
    }
    
    this.get_errxesm = function(){
	return this.dxesm_dyxmsmsp() * this.err_yxmsmsp;
    }
    
    this.get_sesm = function(){ // equation 12
	var xesm = this.get_xesm();
	return  (xesm+this.ytsp/this.yxtsp)*this.sxmsmsp*this.sxtss*this.xsmss;
    }
    
    this.get_errsesm = function(){
	var EB = this.SigmaB();
	var J = this.Jsesm();
	return propagate(J,EB);
    }
        
    this.dxesm_dyxmsmsp = function(){
	var out = 0;
	if (this.mydata.inSpike(this.X)){
	    out = (this.ytsp*this.xytsp)*
		(1+(this.yxmsmsp-this.yxtsp)/(this.yxtsm-this.yxmsmsp))/
		(this.yxmsmsp-this.yxtsm);
	}
	return out;
    }
    
    this.Jsesm = function(){
	return [this.dsesm_dyxmsmsp(),this.dsesm_dsxmsmsp(),this.dsesm_dxsmss()];
    }
    
    this.dsesm_dxesm = function(){
	return this.sxmsmsp*this.sxtss*this.xsmss;
    }
    
    this.dsesm_dxsmss = function(){
	return (this.get_xesm()+this.ytsp*this.xytsp)*this.sxmsmsp*this.sxtss;
    }
    
    this.dsesm_dsxmsmsp = function(){
	return (this.get_xesm()+this.ytsp*this.xytsp)*this.sxtss*this.xsmss;
    }
    
    this.dsesm_dyxmsmsp = function(){
	return this.sxmsmsp*this.sxtss*this.xsmss*this.dxesm_dyxmsmsp();
    }

    this.SigmaA = function(){
	var EB = [[Math.pow(this.err_xsmsmsp,2) , 0                  ],
		  [0                      , Math.pow(this.sxmss[1],2)]];
	return EB;
    }
    
    // equation 24
    this.SigmaB = function(){
	var var_yxmsmsp = Math.pow(this.err_yxmsmsp,2);
	var var_sxmsmsp = Math.pow(this.err_sxmsmsp,2);
	var var_xsmss = this.var_xsmss;
	var EB = [[var_yxmsmsp,       this.cov_yxmsmsp_sxmsmsp,0        ],
		  [this.cov_yxmsmsp_sxmsmsp,var_sxmsmsp,       0        ],
		  [0,                 0,                 var_xsmss]];
	return EB;
    }

    this.dCX_dxesm = function(){
	var xesm = this.get_xesm();
	var sesm = this.get_sesm();
	return getCSMXAsMSAx(this.X,this.x,this.S,this.s) * 
            (1 - this.dsesm_dxesm()*xesm/sesm) / sesm;
    }
    
    this.dCX_dsesm = function(){
	var xesm = this.get_xesm();
	var sesm = this.get_sesm();
	return getCSMXAsMSAx(this.X,this.x,this.S,this.s) * xesm/Math.pow(sesm,2);
    }

    this.JA = function(){
	var CSMXAsMSAxxstss = getCSMXAsMSAx(this.X,this.x,this.S,this.s)*this.xstss;
	var dCXdxsmsm = CSMXAsMSAxxstss*this.sxmss[0];
	var dCXdsxmss = CSMXAsMSAxxstss*this.xsmsmsp;
	return [dCXdxsmsm,dCXdsxmss];
    }
    
    this.JB = function(){
	var dCXdxesm = this.dCX_dxesm();
	var dCXdsesm = this.dCX_dsesm();
	var dCXdyxmsmsp = dCXdxesm*this.dxesm_dyxmsmsp();
	var dCXdsxmsmsp = dCXdsesm*this.dsesm_dsxmsmsp();
	var dCXdxsmss   = dCXdsesm*this.dsesm_dxsmss();
	return [dCXdyxmsmsp,dCXdsxmsmsp,dCXdxsmss];
    }

}
