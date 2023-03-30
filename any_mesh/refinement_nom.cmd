// refinement.3.cmd
// Local refinement of an equilateral mesh, with freeedgebend correction.

// Set form factors from ref_coords.
set_form_factors := {
    foreach facet ff do {
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
    }
}
set_form_factors;

function real get_amax(integer num, real ldelta) {
    return 2/pi/num * sqrt(- ldelta * wd^2 * (1 + ldelta) ) * 1.01;
}

/*
function real cryz(real e2, real e3, real tt, real el){
    local cp := e2 * tt + abs(e3)
    cp < 0 ? 0 : (e2 * tt + abs(e3))^2
    //return e3 < 0 ? (- e2 * t2 + e3 * t1)/el : e2 * t2 + e3 * t1) : (e3 < 0 : - e2 * t2 - e3 * t1)
    //y1 * z2 < minimum(z1, -z1);
}
*/

procedure set_amax(integer lminnum, integer lmaxnum, real ldelta) {
    minnum := lminnum;
    maxnum := lmaxnum;
    if lminnum == lmaxnum then {
        amax := (amax < 0 ? -1 : 1) * get_amax(lmaxnum, ldelta);
    } else {
        amin :=  (amin < 0 ? -1 : 1) * get_amax(lmaxnum, ldelta);
        amax :=  (amax < 0 ? -1 : 1) * get_amax(lminnum, ldelta);
        af := (abs(amin) - abs(amax)) / lngth;
        as := lngth * abs(amax) / (abs(amin) - abs(amax));
        tti := abs(amax) * minnum * pi/ww * 1.05;
        tta := abs(amin) * maxnum * pi/ww * 1.05;
    };
}

procedure set_ksub(integer lminnum, integer lmaxnum) {
    //set facet ksub (wnum * pi/wd)^4 * bend.modulus/2;
    local maxksub;
    local minksub;

    minnum := lminnum;
    maxnum := lmaxnum;
    if lminnum == lmaxnum then { 
        kf := (maxnum * 2 * pi/(wd * (1+delta)))^4 * bend.modulus/4;
        //amax := get_amx(maxnum, delta);
    } else {
        maxksub := (maxnum * 2 * pi/(wd * (1+delta)))^4 * bend.modulus/4; //remember the K/2 z^2
        minksub := (minnum * 2 * pi/(wd * (1+delta)))^4 * bend.modulus/4;
        kf := (maxksub - minksub)/lngth;
        ks := lngth * minksub / (maxksub - minksub);
        //amax := get_amx(minnum, delta);
    };
    set_amax(lminnum, lmaxnum, delta);
    recalc;
}

procedure set_delta(real dtt) {
    delta := dtt;
    upEndY := wd * (1 + delta);
    ww := upEndY;
    set_amax(minnum, maxnum, delta);
    //amax := get_amx(minnum, delta);
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
    //gbend.modulus := -thicknezz**2 /12/(1 + pr);
    //edgewidth := 1/3;  // Right triangles along edge
    //edgewidth := sqrt(3)/4;  // Equilateral triangles along edge
    //freeedgebend.modulus := bend.modulus/4 * (1 - pr**2) * edgewidth * gridsize;
    set_ksub(minnum, maxnum);
    recalc;
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
        vv.ref_coord[1] :=
            avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0,
                    ref_coord[1]));
        vv.ref_coord[2] :=
            avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0,
                    ref_coord[2]));
        vv.ref_coord[3] :=
            avg(vv.edge ee where old_eid != 0,sum(ee.vertex where old_vid != 0,
                    ref_coord[3]));
        // Set bending energies for those not on an edge:
        if max(vv.edge, border) == 0 then {
            set vv bend;
            set vv gbend;
        };
        //print "set refine";
        
        //else
        //    set vv freeedgebend;
    };
    //exec "foreach vertex vv where old_vid == 0 do { if ! (vv.on_boundary 1 || vv.on_boundary 2) then { if is_defined(\"egap\") then  { if edgegap then { if vv.ref_coord[1] == lngth then set vv egap;} else {set vv egap;};};};}; ";

    // Only swap edges that were newly added or were heights of right triangles

    do {
        flush_counts;
        equiangulate edge where rheight or old_eid==0;
    } while equi_count;

    if edgepr then {
        if is_defined("cenergy") then {
            exec "unset facet cenergy";
            exec "foreach facet ff do { if min(ff.vertex, ref_coord[1]) == 0 || max(ff.vertex, ref_coord[1]) == lngth then set ff cenergy;};";
        };
    };

    set_form_factors;
    // Edge correction needs to be updated:
    set_thickness(thicknezz);
}

fix_bound := {
    foreach edge ee where border and ! (on_boundary 1 || on_boundary 2) do {
        fix ee;
        fix ee.vertex;
    };
}

unfix_bound := {
    foreach edge ee where border and ! (on_boundary 1 || on_boundary 2) do {
        unfix ee;
        unfix ee.vertex where ! (on_boundary 1 || on_boundary 2) || ref_coord[1] == lngth;
    };
}

r :::= {
    set edge divide 1 where length > (1 + 2/sqrt(3))/2 * gridsize;
    refinemarked();
    gridsize := gridsize/2;
    if is_defined("fa") then exec "fa := fa / 4";
    recalc;
    refine_times++;
}

procedure set_bwave(integer lminnum, integer lmaxnum) {
    //local ww;
    //ww := wd * (1 + delta);
    foreach vertex vv where ref_coord[1] == 0 and (! (on_boundary 1 || on_boundary 2)) do {
        vv.z := 2/pi/lminnum * sqrt(- delta * wd^2 * (1 + delta) ) * sin(lminnum * pi / ww * vv.y);
    };
    foreach vertex vv where ref_coord[1] == lngth and (! (on_boundary 1 || on_boundary 2)) do {
        vv.z := 2/pi/lmaxnum * sqrt(- delta * wd^2 * (1 + delta) ) * sin(lmaxnum * pi / ww * vv.y);
    };
}

set_wave := {
    local ww;
    local aa;
    ww := wd * (1 + delta);
    aa := 2/pi/maxnum * sqrt(- delta * wd^2 * (1 + delta) );
    foreach vertex vv where ! (on_boundary 1 || on_boundary 2) do {
        vv.z := aa * sin(maxnum * pi / ww * vv.y);
    };
}

set_cwave := {
    local ww;
    local aa;
    local lnm;
    ww := wd * (1 + delta);
    foreach vertex vv where ! (on_boundary 1 || on_boundary 2) do {
        lnm := ((kf * (vv.x + ks)) * 4 / bend.modulus)^(1/4) * (wd * (1+delta));
        aa := 2/lnm * sqrt(- delta * wd^2 * (1 + delta) );
        vv.z := aa * sin(lnm / ww * vv.y);
    };
}

set_awave := {
    local ww;
    local aa;
    local lnm;
    ww := wd * (1 + delta);
    foreach vertex vv where ! (on_boundary 1 || on_boundary 2) do {
        lnm := ((kf * (vv.x + ks)) * 4 / bend.modulus)^(1/4) * (wd * (1+delta));
        aa := 2/lnm * sqrt(- delta * wd^2 * (1 + delta) );
        vv.z := aa / 2 * (sin(lnm / ww * vv.y) + sin(lnm / ww * (ww - vv.y)));
    };
}

// set_prewavlen := {
//     (maxnum * pi/(wd * (1+delta)))^4 * bend.modulus/4;
// }

//unrefine := {
//    delete edges where old_eid == 0;
    //u; V; u; V;
//    recalc;
//    set_form_factors;
//    set_thickness(thicknezz);
//}

function real ssqrt(real val) {
    if val < 0 then
        return 0;
    return sqrt(val);
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



