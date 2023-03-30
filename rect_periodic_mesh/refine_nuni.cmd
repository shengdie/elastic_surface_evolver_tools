// refinement.3.cmd
// Local refinement of an equilateral mesh, with freeedgebend correction.

// Set form factors from ref_coords.
set_form_factors := {
    local vp1;
    local vp2;
    local vp3;
    foreach facet ff do {
        if max(ff.edge, wrap) == 0 then {
            set ff.form_factors[1]
                (ff.vertex[2].ref_coord[1] - ff.vertex[1].ref_coord[1])^2
            + (ff.vertex[2].ref_coord[2] - ff.vertex[1].ref_coord[2])^2
            + (ff.vertex[2].ref_coord[3] - ff.vertex[1].ref_coord[3])^2;
            set ff.form_factors[2]
                (ff.vertex[2].ref_coord[1] - ff.vertex[1].ref_coord[1])
            *(ff.vertex[3].ref_coord[1] - ff.vertex[1].ref_coord[1])
            + (ff.vertex[2].ref_coord[2] - ff.vertex[1].ref_coord[2])
            *(ff.vertex[3].ref_coord[2] - ff.vertex[1].ref_coord[2])
            + (ff.vertex[2].ref_coord[3] - ff.vertex[1].ref_coord[3])
            *(ff.vertex[3].ref_coord[3] - ff.vertex[1].ref_coord[3]);
            set ff.form_factors[3]
                (ff.vertex[3].ref_coord[1] - ff.vertex[1].ref_coord[1])^2
            + (ff.vertex[3].ref_coord[2] - ff.vertex[1].ref_coord[2])^2
            + (ff.vertex[3].ref_coord[3] - ff.vertex[1].ref_coord[3])^2;
        } else {
            vp1 := (ff.vertex[1].ref_coord[2] == 0) * wd;
            vp2 := (ff.vertex[2].ref_coord[2] == 0) * wd;
            vp3 := (ff.vertex[3].ref_coord[2] == 0) * wd;
            set ff.form_factors[1]
                (ff.vertex[2].ref_coord[1] - ff.vertex[1].ref_coord[1])^2
            + (ff.vertex[2].ref_coord[2] - ff.vertex[1].ref_coord[2] + vp2 - vp1)^2
            + (ff.vertex[2].ref_coord[3] - ff.vertex[1].ref_coord[3])^2;
            set ff.form_factors[2]
                (ff.vertex[2].ref_coord[1] - ff.vertex[1].ref_coord[1])
            *(ff.vertex[3].ref_coord[1] - ff.vertex[1].ref_coord[1])
            + (ff.vertex[2].ref_coord[2] - ff.vertex[1].ref_coord[2] + vp2 - vp1)
            *(ff.vertex[3].ref_coord[2] - ff.vertex[1].ref_coord[2] + vp3 - vp1)
            + (ff.vertex[2].ref_coord[3] - ff.vertex[1].ref_coord[3])
            *(ff.vertex[3].ref_coord[3] - ff.vertex[1].ref_coord[3]);
            set ff.form_factors[3]
                (ff.vertex[3].ref_coord[1] - ff.vertex[1].ref_coord[1])^2
            + (ff.vertex[3].ref_coord[2] - ff.vertex[1].ref_coord[2] + vp3 - vp1)^2
            + (ff.vertex[3].ref_coord[3] - ff.vertex[1].ref_coord[3])^2;
        };
    }
}

function real get_a(integer num, real ldelta) {
    return 1/pi/num * sqrt(- ldelta * wd^2 * (1 + ldelta) ) * 1.01;
}

set_aa := {
    a0 := 1/pi/num0 * sqrt(- delta * wd^2 * (1 + delta));
    al := 1/pi/numl * sqrt(- delta * wd^2 * (1 + delta) );
    recalc;
}

procedure set_a(real lnum0, real lnuml, real ldelta) {
    num0 := lnum0;
    numl := lnuml;
    if lnum0 == lnuml then {
        a0 := (a0 < 0 ? -1 : 1) * get_a(lnuml, ldelta);
    } else {
        al :=  (al < 0 ? -1 : 1) * get_a(lnuml, ldelta);
        a0 :=  (a0 < 0 ? -1 : 1) * get_a(lnum0, ldelta);

        //af := (abs(al) - abs(a0)) / lngth;
        //as := lngth * abs(a0) / (abs(al) - abs(a0));
        //tti := abs(a0) * num0 * pi/ww * 1.05;
        //tta := abs(al) * numl * pi/ww * 1.05;
        //tti := 2 * sqrt( - ldelta / (1+ldelta));
        //tta := 2 * sqrt( - ldelta / (1+ldelta));
        //ttt := 2 * sqrt( - ldelta / (1+ldelta));
    };
}

function real get_ksub(real lnum) {
    return (2 * lnum * pi / ww)^4 * bend.modulus/4; //remember the K/2 z^2
}

procedure set_ksub(real lnum0, real lnuml) {
    num0 := lnum0;
    numl := lnuml;

    if lnum0 == lnuml then { 
        print "uniform\n";
        wl0 := ww / lnum0;
        wll := ww / lnuml;
        subenergy.modulus := get_ksub(num0);
        //a0 := get_amx(numl, delta);
    } else {
        wl0 := ww / lnum0;
        wll := ww / lnuml;
        kl := get_ksub(numl); 
        k0 := get_ksub(num0);
        kf := (kl - k0) / lngth;

        ls := wl0;
        
        if is_defined("righttension") then {
            lf := (wll-wl0) / (lngth * (1+tensionG)); 
        } else {
            lf := (wll - wl0) / lngth;
        };
        // lambda = lf (x + ls)
        // lf := (wll - wl0) / lngth; 
        // ls := (lngth * wl0) / (wll - wl0)
        // E_s = 1/2 * K z^2, K = 1/2 B / (wave_len)^4 * (2 *pi)^4  =  B/ lf^4 / (x+ls)^4 * 16 *pi^4/2 -> sube.modulus := B / lf^4 *8 *pi^4
        // E_s = 1/2 * K z^2, K = 1/2 B / (wave_len)^4 * (2 *pi)^4 = B / (l0 + lf * x)^4 * 16 * pi^4/2 -> sube.modulus := B *8 *pi^4
        subenergy.modulus := bend.modulus * 4 * pi^4; // B=bend.modulus/2
        
    };
    //set_a(lnum0, lnuml, delta);
    recalc;
}

procedure set_ktimes(real ltimes, real lnuml) { //set the d wave_len / dx = ltimes * d old_wav_len /dx, and set numl to lnuml
    wll := ww / lnuml;
    numl := lnuml;
    lf := ltimes * lf;
    ls := (wll - lngth * lf) / lf;
    wl0 := lf * ls;
    num0 := ww / wl0;
    //set_a(num0, numl, delta);
    //set_ksub(num0, numl);
    subenergy.modulus := bend.modulus / lf^4 * 4 *pi^4;
    recalc;
}

procedure set_delta(real dtt) {
    delta := dtt;
    upEndY := wd * (1 + delta);
    if is_defined("righttension") then upEndY := upEndY * sqrt(1 - facet[1].poisson_ratio * tensionG * 2);
    ww := upEndY;
    ttt := 2 * sqrt( - delta / (1+delta));
    //set_a(num0, numl, delta);
    //a0 := get_amx(num0, delta);
    recalc;
}

procedure set_adelta(real dtt) {
    foreach vertex vv do {
        vv.y := vv.y *  (1+dtt) /(1+delta);
        if delta < 0 then vv.z := vv.z * (1+dtt) /(1+delta) * sqrt(dtt/delta);
    };
    foreach vertex vv where on_boundary 1 || on_boundary 2 do {
        vv.p1 := vv.p1 *  (1+dtt) /(1+delta);
        if delta < 0 then vv.p2 := vv.p2 * (1+dtt) /(1+delta) * sqrt(dtt/delta);
    };
    
    delta := dtt;
    upEndY := wd * (1+delta);
    if is_defined("righttension") then upEndY := upEndY * sqrt(1 - facet[1].poisson_ratio * tensionG * 2);
    ww := upEndY;
    ttt := 2 * sqrt( - delta / (1+delta));
    
    //set_a(num0, numl, delta);
    recalc;
}

procedure set_thickness(real th) {
    local pr;

    thicknezz := th;
    // Assume all facets have the same Poisson ratio
    pr := facets[1].poisson_ratio;
    bend.modulus := thicknezz**2 / 6 / (1 - pr**2);
    //gbend.modulus := -bend.modulus * (1 - pr)/2;
    if is_defined("gbend") then exec "gbend.modulus := -thicknezz**2 /12/(1 + facets[1].poisson_ratio);";
    //edgewidth := 1/3;  // Right triangles along edge
    //edgewidth := sqrt(3)/4;  // Equilateral triangles along edge
    //freeedgebend.modulus := bend.modulus/4 * (1 - pr**2) * edgewidth * gridsize;
    set_ksub(num0, numl);
    recalc;
}

procedure set_hdwl(real lhwl, real num) {
    local lwll;
    lwll := ww / num;
    set_thickness(lwll * lhwl);
    set_ksub(num, num);
    recalc;
}

procedure set_thd(real lnum0, real lnuml, real ltimes) {   // set thresold delta, so that = 0.9 * delta_fold and is ltimes times of delta_critical
    num0 := lnum0;
    numl := lnuml;
    set_adelta(-0.9/numl^2);
    set_thickness(sqrt(3/2 * abs(delta) / ltimes) / numl / pi * wd);
    // set_thickness(sqrt(12 * (1 - facets[1].poisson_ratio^2) / ltimes/ sqrt(2) / numl^2 / (2 * pi * numl/ ww)^2));
    set_ksub(num0,numl);
    recalc;
}

all_thd := {
    local dc_local;
    local lratio;
    k0 := get_ksub(num0);
    dc_local := 2 * sqrt(bend.modulus * get_ksub(numl)) * (1-facet[1].poisson_ratio^2);
    lratio := (0.5901888630769108 * (numl - num0) * wd * abs(delta)**(1 / 4) / (bend.modulus * k0)**(1 / 8) / num0 / numl / lngth)**4;
    printf "delta_c: %g, delta / delta_c: %g\n", dc_local, -delta / dc_local;
    printf "lratio: %g, sqrt ratio: %g\n", 1/lratio, 1/sqrt(lratio);
}

//procedure set_tension(real ten) {
//    tensionG := ten;
//    righttension.modulus := tensionG;
//    fLen := sqrt(1+2*tensionG) * lngth;
//    recalc;    
//    printf "Setted to new tension: %g\n", tensionG;
//}

procedure refinemarked(){

    local newmark;
    local redges;
    define redges real[3];
    local rels;
    local tp; //temp
    define rels real[3];  //reference length of triangle
    local inx;
    local iny;
    
    local hypid;
    local baseid;
    local heightid;
    local otherhypid;
    local onedge;

    set edge rheight 0;

    do {
        newmark := 0;
        foreach facet ff where sum(ff.edges, divide) do {
            // If it is an equilateral triangle
            if abs(ff.form_factors[1] - ff.form_factors[3]) < 0.01*ff.form_factors[1] then {
                // If 2/3 of the sides are to be divided, divide all
                if sum(ff.edges, divide) == 2 then {
                    set ff.edges divide 1;
                    newmark := 1;
		    //printf "new in triangle facet %d\n",ff.id;
                };
            } else {
		redges[1] := ff.edges[1].id;  // redges keep the id of three edges, e1 = v1 -> v2; e2 = v2 -> v3; e3 = v3 -> v1;
                redges[2] := ff.edges[2].id;
                redges[3] := ff.edges[3].id;
                rels[1] := ff.form_factors[1];  // rels keep the ref_length^2, fac1 = (v2-v1)^2 = e1^2; fac2 = (v3-v1)*(v2-v1) = e1*(-e3); fac3 = (v3-v1)^2  = e3^2;
                rels[2] := ff.form_factors[3] + ff.form_factors[1] - 2 * ff.form_factors[2];
                rels[3] := ff.form_factors[3]; 

		for(inx := 1; inx < 3; inx++) {
		   for(iny := 1; iny < 4-inx; iny++) {
                      if rels[iny] > rels[iny+1] then {
                         tp := rels[iny+1];			 
                         rels[iny+1] := rels[iny];
			 rels[iny] := tp;
			 tp := redges[iny+1];
			 redges[iny+1] := redges[iny];
			 redges[iny] := tp;
			 
		      };
		      //print iny;
		   };
		};

		baseid := redges[1]; //printf "baseid: %d\n",baseid;
		heightid := redges[2]; //printf "heightid: %d\n",heightid;
		hypid := redges[3]; //printf "hypid: %d\n",hypid;

		edges[heightid].rheight := 1;

		//Make sure the hypotenuse
		if not edges[hypid].divide then
		{
		   set edges[hypid] divide 1; // divide the hypotenuse
		   newmark := 1;
		   //printf "new in right triangle %d, hypotenuse %d\n",ff.id,hypid;
                };      

		//Divide the height if it's on border
		if edges[heightid].border then	
		{  
		   if not edges[heightid].divide then
		   {
			set edges[heightid] divide 1; // divide the height
			newmark := 1;
			//printf "new in right triangle %d, height %d\n",ff.id,heightid;
		   };
		}  //otherwise, don't divide it, divide the hypotenuse of the pair triangle
		else
		{
		   if edges[heightid].divide then
		   {
		   	set edges[heightid] divide 0;
			set edges[heightid] color blue;
			//printf "delete divide of height not on border of right triangle %d, edge %d\n",ff.id,heightid;
		   };

		   foreach edges[heightid].facet ef where id != ff.id do {
			local eid;
			eid := ef.id;  //ef's edges order is different from facet[ef.id], which should be a bug of evolver
			redges[1] := facet[eid].edges[1].id;  // redges keep the id of three edges, e1 = v1 -> v2; e2 = v2 -> v3; e3 = v3 -> v1;
                	redges[2] := facet[eid].edges[2].id;	      
                	redges[3] := facet[eid].edges[3].id;
                	rels[1] := facet[eid].form_factors[1];  // rels keep the ref_length^2, fac1 = (v2-v1)^2 = e1^2; fac2 = (v3-v1)*(v2-v1) = e1*(-e3); fac3 = (v3-v1)^2  = e3^2;
                	rels[2] := facet[eid].form_factors[3] + facet[eid].form_factors[1] - 2 * facet[eid].form_factors[2];
                	rels[3] := facet[eid].form_factors[3];

                        //print rels[1]; print rels[2]; print rels[3];
			/*
			printf "%g,%d\n",rels[1],redges[1];
			printf "%g,%d\n",rels[2],redges[2];
			printf "%g,%d\n",rels[3],redges[3];
			print facet[eid].edge[1].id;
			print facet[eid].edge[2].id;
			print facet[eid].edge[3].id;
			*/

			for(inx:=1;inx<3;inx++) {
			   if rels[inx] > rels[inx+1] then {
			      tp := rels[inx+1];
			      rels[inx+1] := rels[inx];
			      rels[inx] := tp;
			      tp := redges[inx+1];
			      redges[inx+1] := redges[inx];
			      redges[inx] := tp;
			   };			
			};
			/*
                        printf "%g,%d\n",rels[1],redges[1];
                        printf "%g,%d\n",rels[2],redges[2];
                        printf "%g,%d\n",rels[3],redges[3];
			*/
			//print rels[1]; print rels[2]; print rels[3];
			otherhypid := redges[3];
			if not edges[otherhypid].divide then {
                           set edges[otherhypid] divide 1;
                           newmark := 1;
			   //printf "new in pair triangle %d of rtri %d, otherhypid %d\n",ef.id,ff.id,otherhypid;
			   set edges[otherhypid] color red;
			   //pause;
                    	};
                   };
		};

		// If two bases are touching, or base is on perimeter, don't divide it
          	if edges[baseid].divide then {
                    onedge := 1;
                    // Get other face
                    foreach edges[baseid].facets ef where id != ff.id do {
                        onedge := 0;
                        // Is it a right triangle?
                        if abs(ef.form_factors[1] - ef.form_factors[3]) > 0.01*ef.form_factors[1] then
                            // Is baseid also a base here?
                            if edges[baseid].length == min(ef.edges, length) && edges[baseid].divide then
                               { set edges[baseid] divide 0; /*printf "delete one divide because touch of right triangle %d, edge %d\n",ff.id,baseid;*/ }
                    };
                    if onedge && edges[baseid].divide then
                        { set edges[baseid] divide 0; /*printf " delete divide because onedge base of right triangle %d, edge %d\n",ff.id, baseid;*/ }
                };
		
            };
        };
    } while newmark;

    set vertex old_vid id;
    set edge old_eid id;

    //'r';
    refine edges where divide;
    set edge divide 0;
    //gridsize := gridsize/2;

    foreach vertex vv where old_vid == 0 do {
        
        if max(vv.edge where old_eid != 0, wrap) > 0 then {

            vv.ref_coord[1] := avg(vv.edge ee where old_eid != 0, sum(ee.vertex where old_vid != 0, ref_coord[1]));
            vv.ref_coord[2] :=
                avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0, ref_coord[2])) + 1/count(vv.edge ee where old_eid != 0, 1) * wd;
            vv.ref_coord[3] :=
                avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0, ref_coord[3]));
        } else {
        
            vv.ref_coord[1] :=
                avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0,
                        ref_coord[1]));
            vv.ref_coord[2] :=
                avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0,
                        ref_coord[2]));
            vv.ref_coord[3] :=
                avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0,
                        ref_coord[3]));
        };
    };
    //exec "foreach vertex vv where old_vid == 0 do { if ! (vv.on_boundary 1 || vv.on_boundary 2) then { if is_defined(\"egap\") then  { if edgegap then { if vv.ref_coord[1] == lngth then set vv egap;} else {set vv egap;};};};}; ";

    // Only swap edges that were newly added or were heights of right triangles

    do {
        flush_counts;
        equiangulate edge where rheight or old_eid==0;
    } while equi_count;

    set_form_factors;
    // Edge correction needs to be updated:
    set_thickness(thicknezz);
}

fix_bound := {
    foreach edge ee where border do {
        fix ee;
        fix ee.vertex;
    };
}

unfix_bound := {
    foreach edge ee where border do {
        unfix ee;
        unfix ee.vertex where ! (ref_coord[1]==0 && ref_coord[2] == 0); //where ! (on_boundary 1 || on_boundary 2) || ref_coord[1] == lngth;
    };
}

r :::= {
    //set edge divide 1 where length > (1 + 2/sqrt(3))/2 * gridsize;
    set edge divide 1;
    refinemarked();
    gridsize := gridsize/2;
    if is_defined("fa") then exec "fa := fa / 4";
    recalc;
    refine_times++;
}

procedure set_bwave(integer lnum0, integer lnuml) {
    //local ww;
    //ww := wd * (1 + delta);
    foreach vertex vv where ref_coord[1] == 0 do {
        vv.z := 1/pi/lnum0 * sqrt(- delta * wd^2 * (1 + delta) ) * sin(lnum0 * pi *2 / ww * vv.y);
    };
    foreach vertex vv where ref_coord[1] == lngth do {
        vv.z := 1/pi/lnuml * sqrt(- delta * wd^2 * (1 + delta) ) * sin(lnuml * pi *2 / ww * vv.y);
    };
}

procedure set_shift(real shift) {
    foreach vertex vv where ref_coord[1] == 0 do {
        vv.z := 1/pi/num0 * sqrt(- delta * wd^2 * (1 + delta)) * sin(num0 * pi *2 / ww * vv.y);
    };
    foreach vertex vv where ref_coord[1] == lngth do {
        vv.z := 1/pi/numl * sqrt(- delta * wd^2 * (1 + delta) ) * sin(numl * pi *2 / ww * (vv.y + shift * wll));
    };
}

set_twave := {
    if is_defined("rEndX") then {
        exec "foreach vertex vv do {vv.z := 1/pi/num0 * sqrt(- delta * wd^2 * (1 + delta)) * sin(num0 * pi *2 / ww * vv.y + shift/rEndX * x);};"
    };
}

set_prewlen := {
    if is_defined("rEndX") then {
        exec "subenergy.modulus := (2 * num0 * pi / ww)^4 * bend.modulus/4 * (1 + ww^2 /4/ num0^2/rEndX^2)^2;";
    };  
}

set_wave := {
    //local ww;
    local aa;
    //ww := wd * (1 + delta);
    aa := 1/pi/numl * sqrt(- delta * wd^2 * (1 + delta) );
    foreach vertex vv do {
        vv.z := aa * sin(numl * pi * 2 / ww * vv.y);
    };
    if valid_boundary(1) then set_boundary(numl,numl);
}

procedure set_boundary(real lnum0, real lnuml) {
    local aa;
    aa := 1/pi/lnum0 * sqrt(- delta * wd^2 * (1 + delta) );
    foreach vertex vv where on_boundary 1 do {
        set vv p2 aa * sin(lnum0 * pi * 2 / ww * vv.p1);
    };
    aa := 1/pi/lnuml * sqrt(- delta * wd^2 * (1 + delta) );
    foreach vertex vv where on_boundary 2 do {
        set vv p2 aa * sin(lnuml * pi * 2 / ww * vv.p1);
    };
}

set_cwave := {
    //local ww;
    local aa;
    local lnm;
    //ww := wd * (1 + delta);
    foreach vertex vv do {
        //lnm := ((kf * (vv.x + ks)) * 4 / bend.modulus)^(1/4) * (wd * (1+delta));
        lnm := ww / (lf * (vv.x+ls));
        aa := 2/lnm * sqrt(- delta * wd^2 * (1 + delta) );
        vv.z := aa * sin(lnm / ww * vv.y);
    };
}

set_awave := {
    //local ww;
    local aa;
    local lnm;
    //ww := wd * (1 + delta);
    foreach vertex vv do {
        lnm := ((kf * vv.x + k0) * 4 / bend.modulus)^(1/4) * ww;
        //lnm := ww / (lf * (vv.x+ls));
        aa := 2/lnm * sqrt(- delta * wd^2 * (1 + delta) );
        vv.z := aa / 2 * (sin(lnm / ww * vv.y) + sin(lnm / ww * (ww - vv.y)));
    };
}

function real ssqrt(real val) {
    if val < 0 then
        return 0;
    return sqrt(val);
}

correct_tension := {
    local actuallen;
    actuallen := sqrt(1+tensionG*2) * lngth * 1.001;
    printf "current len: %g, expected: %g\n", rEndX, actuallen;
    // rEndX := sqrt(1+tensionG*2) * lngth ;
    foreach vertex vv do vv.x := actuallen / rEndX * vv.x;
    rEndX := actuallen;
}

fix_len := {
    local actuallen;
    actuallen := sqrt(1+tensionG*2) * lngth;
    printf "current len: %g, expected: %g\n", rEndX, actuallen;
    // rEndX := sqrt(1+tensionG*2) * lngth ;
    foreach vertex vv do vv.x := actuallen / rEndX * vv.x;
    rEndX := actuallen;
    fix rEndX;
}

show_len := {
    printf "current len: %g, expected: %g\n", rEndX, sqrt(1+tensionG*2) * lngth;
}

procedure refine_range(real xi, real xf) {
    foreach vertex vv where vv.x >= xi and vv.x <= xf do {
        set vv.edges divide 1;
    };
    if sum(edges, divide) > 0 then {
        refinemarked();
        max_refine_times++;
    };
}

procedure mark_out(real zz) {
    foreach vertex vv do {
        if vv.z > zz then
            set vv.edges divide 1;
    };
    if sum(edges, divide) > 0 then {
        refinemarked();
        //return 1;
    };
}

procedure mark_outofpalne(real rat) {
    local maxz;
    foreach vertex vv where x==0 && y==0 do {
        maxz := vv.z;   
    };
    foreach vertex vv where abs(z) > maxz*rat do {   
        set vv.edges divide 1;
    };
    if sum(edges, divide) > 0 then {
        refinemarked();
        //return 1;
    };
}

procedure mark_stretch(real frac) {
    // Mark edges where jump in stretch density > frac * max energy density
    local maxen;
    maxen := max(facets, stretch/area);
    // For some reason, max(ee.facets, stretch) will report erroneous results,
    // while max(ee.facets, facets[id].stretch) works properly....
    set edge ee divide 1 where max(ee.facets, facets[id].stretch/area)
        - min(ee.facets, facets[id].stretch/area) > maxen * frac;
}





