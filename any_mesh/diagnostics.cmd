// diagnostics.cmd

gg := { g 100; hessian_seek; }

showx0 := {foreach vertex vv where ref_coord[1]==0 do set vv.edge color red}

colorprofile := {
  local minz;
  local maxz;
  minz := min(vertex, z);
  maxz := max(vertex, z);
  foreach facet ff do set ff color ((min(ff.vertex, z) - minz) * 15 / (maxz - minz) + 0.5) idiv 1 
}

procedure list_z(real xx) {
    local scaleee;
    scaleee := gridsize/sqrt(3)*2/2;
    foreach vertex vv where abs(ref_coord[1] - xx) < scaleee do {
        //if
        printf "%g\t %g\t %g\n",vv.x, vv.y, vv.z;
    };   
}

procedure g_to(real relenergy) {
  local oee;
  oee := total_energy; 
  g; 
  printf "Doing g till %g\n", relenergy;
  if abs(scale) > 0 then {do {oee := total_energy; g3; } while (oee - total_energy)/oee > relenergy};
  //printf "g diff: %g", (oee - total_energy)/oee;
}

procedure hessian_to(real relenergy) {
  local oe;
  oe := total_energy;
  //hessian_seek; 
  printf "Doing hessian till %g\n", relenergy;
  if abs(last_hessian_scale) > 0 then {do {oe := total_energy; hessian_seek} while (oe - total_energy)/oe > relenergy };
}

save_proc := {dump sprintf "%s_h%g_delta%g_r%g_save.dmp", origin_name, thicknezz/in_r, delta, refine_times; }

procedure hgb(real relenergy) { // smart choosing between hessian and g,
  local oe;
  local oee;
  local oe1;
  local oeg;
  local oeh;
  local clc;
  local clcg;
  local clch;
  local gorh; // == 1 for g, ==2 for hessian
  local gtt; // g times, break continuous g
  //exec "metis_factor";
  g;
  gorh := 1;
  gtt := 0;
  if clock < last_stime then {last_stime := clock };

  do {
    oee := total_energy;
    if gorh == 1 then { //last is g, now do hessian, then g
      oe1 := total_energy;
      clc := clock;
      {hessian_seek} 2;
      clch := clock - clc;
      oeh := (oe1 - total_energy)/abs(oe1);
      //gorh := 2;

      //if abs(last_hessian_scale) < 1e-10 then relee := 1e-11;

      oe1 := total_energy;
      clc := clock;
      g 3;
      clcg := clock - clc;
      oeg := (oe1 - total_energy)/abs(oe1);
      gorh := 1;
    } else { //last is hessian, now do g, then hessian
      oe1 := total_energy;
      clc := clock;
      g 3;
      clcg := clock - clc;
      oeg := (oe1 - total_energy)/abs(oe1);
      
      oe1 := total_energy;
      clc := clock;
      {hessian_seek} 2;
      clch := clock - clc;
      oeh := (oe1 - total_energy)/abs(oe1);
      gorh := 2;
    };
    if (oee - total_energy)/abs(oee) > relenergy then {
      do {
        oe := total_energy;
        if oeg/clcg > oeh/clch then 
        {
          oe1 := total_energy;
          clc := clock;
          g 6;
          clcg := clock - clc;
          oeg := (oe1 - total_energy)/abs(oe1);
          gorh := 1;
          gtt := gtt + 6;
          if gtt >= 400 then {{hessian_seek} 6; oeh := oeg * 1000 * clch / clcg; gtt := 0 }
          else {if gtt >= 200 && gtt <= 210  then { 
            oeh := oeg * 1000 * clch / clcg; //gtt := 0 
          };}; // if g for 200 times, then try hessian
        } else {
          oe1 := total_energy;
          clc := clock;
          {hessian_seek} 2;
          clch := clock - clc;
          oeh := (oe1 - total_energy)/abs(oe1);
          gorh := 2;
        }; 
        printf "diff: %g\n", (oe - total_energy)/abs(oe);
        // For saving the process
        if (clock - last_stime) >= save_period then {
          print "Saving the process";
          // dump sprintf "%s_save.dmp",origin_name;
          save_proc;
          last_stime := clock;
        };
      } while (oe - total_energy)/abs(oe) > relenergy;
    };
       
  } while (oee - total_energy)/abs(oee) > relenergy;
}

// converge_to(max allowed relative energy step)
procedure converge_to(real relenergy, integer hnum, integer gnum) { //convert to relative energy, hessian num, g num, saddle num

    if clock idiv (run_time*60) == ( submit_times - ost ) then {
           submit_times++;
           //dump sprintf "%s_s%d.dmp",origin_name,submit_times;
           save_proc;
           //print thicknezz | "sleep 260";
           //quit;
    };

    //local thseek; //the time needed for 1 hessian_seek;
    local oen;
    local en;
    recalc;
    g;
    if clock idiv (run_time*60) == ( submit_times - ost ) then {
           submit_times++;
           //dump sprintf "%s_s%d.dmp",origin_name,submit_times;
           save_proc;
           //print thicknezz | "sleep 260";
           //quit;
    };
    thseek := clock;
    oen := total_energy;
    {hessian_seek} hnum;
    g gnum;
    en := total_energy;
    thseek := clock - thseek;
    
    if ((oen-en)/en) > relenergy then {
        do {
            if clock idiv (run_time*60) == ( submit_times - ost ) then {
                submit_times++;
                //dump sprintf "%s_s%d.dmp",origin_name,submit_times;
                save_proc;
                //print thicknezz | "sleep 200";
                //quit;
                //quit;
                //q;
            };
            thseek := clock;
            oen := total_energy;
            {{hessian_seek} hnum; g gnum; };
            en := total_energy;
            thseek := clock - thseek;
            printf "Difference: %g\n", (oen-en)/abs(en);
            printf "Target:     %g\n", relenergy;
        } while (((oen-en)/en) > relenergy);
    };
}

cc3 := {converge_to(1e-3,1, 3)}
cc6 := {converge_to(1e-6,1, 3)}
cc9 := {converge_to(1e-9,1, 3)}
cc12 := {converge_to(1e-12,1, 3)}
cc0 := {converge_to(0,1, 3)}
ct0 := {converge_to(0,1,0)}

//ss0 := { do { oe := total_energy; {saddle; converge_to(1e-14,1,0)};  } while abs((oe - total_energy)/oe) > 0 }
// saddle to zero
ss0 := { 
  local oe;
  local oee;
  do { 
    oe := total_energy; 
    {saddle; oee := total_energy; g; 
    if (oee - total_energy)/oee > 1e-15 then {do {oee := total_energy; g3; } while (oee - total_energy)/oee > 1e-13}; 
    hgb(1e-14)}; 
    oee := total_energy; g; if (oee - total_energy)/oee > 1e-15 then g2; } while ((oe - total_energy)/oe) > 0 
}
// fix bound then unfix to evolve to zero
/*
ff0 := {
  local oe;
  local oee;
  do {
    oe := total_energy;
    print "fix bound now\n";
    fix_bound;
    //g_to(6e-5);
    print "Doing hessian till 1e-14\n";
    hgb(1e-12);
    oee := total_energy; g; if (oee - total_energy)/oee > 1e-15 then g2;
    print "unfix bound\n";
    unfix_bound;
    //g_to(6e-5);
    print "Doing hessian till 1e-14\n";
    hgb(1e-12);
    oee := total_energy; g; if (oee - total_energy)/oee > 1e-15 then g2;
    printf "current diff: %g", (oe - total_energy)/oe;
  } while ((oe - total_energy)/oe) > 1e-14 
}
*/
// Output the surface - the vertices and facets, their energies, and
// the un-strained configuration

/*
outputsurf := {
    printf "# datafilename: %s\n",datafilename;
    printf "# Surface output: vertices, faces, thickness, gridsize, total_energy, refine_times, poisson_ratio, period or tension\n";
    printf "%g %g %g %0.15g %0.15g %g %g %g\n",
        vertex_count, facet_count, thicknezz, gridsize, total_energy, refine_times, facets[1].poisson_ratio, upEndY;
    printf "# substrate e, stretch energy, bend energy, gbend energy, edge slope e, 0, min_wav_num, max_wav_num\n";
    printf "%0.15g %0.15g %0.15g ", subenergy.value, stretch.value, bend.value;
    if is_defined("gbend") then exec "printf \"%0.15g \", gbend.value;" else printf "0 ";
    if is_defined("righttension") then exec "printf \"%0.15g \", righttension.value" else printf "0 ";
    if is_defined("ese") then exec "printf \"%0.15g \", ese.value " else printf "0 ";
    //if is_defined("shift") then exec "printf \"%0.15g \", shift;" else printf "0 ";
    
    printf "%g %g\n", num0, numl;

    printf "# Vertices: coords, bend, gbend, ref_coords\n";
    foreach vertices do
        printf "%0.15g %0.15g %0.15g %0.15g %0.15g %0.15g %0.15g %0.15g\n",
                x, y, z, bend, 0, ref_coord[1],  ref_coord[2], ref_coord[3];
    printf "# Faces: vertices, stretch, form factors, facet is in periodic b\n";
    foreach facets ff do {
        printf "%g %g %g %0.15g %0.15g %0.15g %0.15g %g\n",
            ff.vertex[1].id, ff.vertex[2].id, ff.vertex[3].id, stretch,
            form_factors[1],  form_factors[2],  form_factors[3], (max(ff.edge, wrap) != 0);
    };
}
*/

outsurfdict := {
    print "# type: dict";
    //printf "# Surface output: vertices, faces, thickness, gridsize, total_energy, refine_times, poisson_ratio, period or tension\n";
    printf "{";
    printf "\"nv\": %d, \"nf\": %d, \"thickness\": %0.15g, \"gridsize\": %0.15g, \"energy\": %0.15g, \"refine_times\": %d, \"pratio\": %g, \"in_r\": %g, \"out_r\": %g, ",
        vertex_count, facet_count, thicknezz, gridsize, total_energy, refine_times, facets[1].poisson_ratio, in_r, out_r;
    printf "\"stretche\": %0.15g, \"bende\": %0.15g, ", stretch.value, bend.value;
    if is_defined("gbend") then exec "printf \"\\\"gbende\\\": %0.15g, \", gbend.value;";
    if is_defined("righttension") then exec "printf \"\\\"tension\\\": %0.15g, \\\"_work\\\": %0.15g, \", righttension.modulus, righttension.value;";
    if is_defined("ese") then exec "printf \"\\\"esee\\\": %0.15g, \", ese.value";
    if is_defined("delta") then exec "printf \"\\\"_delta\\\": %0.15g, \", delta";
    if is_defined("shift") then exec "printf \"\\\"shift\\\": %0.15g, \", shift";
    if is_defined("max_refine_times") then exec "printf \"\\\"max_refine_times\\\": %g, \", max_refine_times";
    //printf "\n# Vertices: coords, bend, gbend, ref_coords\n";
    printf "\"vlist\": [";
    foreach vertices do
        printf "[%0.15g, %0.15g, %0.15g, %0.15g, %0.15g, %0.15g, %0.15g, %0.15g], ",
                x, y, z, bend, 0, ref_coord[1],  ref_coord[2], ref_coord[3];
    printf "], ";
    //printf "\n# Faces: vertices, stretch, form factors, facet is in periodic b\n";
    printf "\"flist\": [";
    foreach facets ff do {
        printf "[%d, %d, %d, %0.15g, %0.15g, %0.15g, %0.15g, %g], ",
            ff.vertex[1].id, ff.vertex[2].id, ff.vertex[3].id, stretch,
            form_factors[1],  form_factors[2],  form_factors[3], (max(ff.edge, wrap) != 0);
    };
    printf "]}\n";
}

/*
fakegbend := { 
    print ((count(vertex where on_constraint 1 and gbend, 1)+count(vertex where on_constraint 2 and gbend,1)-1)*pi + pi/2)*gbend.modulus;
    print total_energy - ((count(vertex where on_constraint 1 and gbend, 1)+count(vertex where on_constraint 2 and gbend,1)-1)*pi + pi/2)*gbend.modulus;
}
*/

/*
outputsurfmirror := {
    printf "# mirror from datafilename: %s\n",datafilename;
    printf "# Surface output: vertices, faces, thickness, gridsize, length, compression, poisson_ratio\n";
    printf "%g %g %g %0.15g %g %g %g 0\n",
        2*vertex_count, 2*facet_count, thicknezz, gridsize, lngth, (upEndY-wd), facets[1].poisson_ratio;
    printf "# Vertices: coords, bend, gbend, ref_coords\n";
    foreach vertices do
        printf "%0.15g %0.15g %0.15g %0.15g %0.15g %0.15g %0.15g %0.15g\n",
                x, y, z, bend, gbend, ref_coord[1],  ref_coord[2], ref_coord[3];
    printf "# Faces: vertices, stretch, form factors, 0\n";
    foreach facets ff do
        printf "%g %g %g %0.15g %0.15g %0.15g %0.15g 0\n",
            ff.vertex[1].id, ff.vertex[2].id, ff.vertex[3].id, 0,
            form_factors[1],  form_factors[2],  form_factors[3];
}
*/

// save the surface to file
//save_surf := { exec sprintf "outputsurf | \"cat > %s_h%g_delta%g_r%g_$(date +%%Y_%%m_%%d_%%H_%%M_%%S).out\"", origin_name,thicknezz/in_r,delta, refine_times }
save_dict := { exec sprintf "outsurfdict | \"cat > %s_h%g_delta%g_r%g_$(date +%%Y_%%m_%%d_%%H_%%M_%%S).out\"", origin_name,thicknezz/in_r,delta, refine_times }


/* used some trick to dump the surface to the {datafilename}_{YYYY_mm_dd_HH_MM_SS} */

dump_surf :=
  { exec sprintf "print \"\" | \"echo dump \\\\\\\"%s_h%g_delta%g_r%g_$(date +\\\"%%Y_%%m_%%d_%%H_%%M_%%S\\\").dmp\\\\\\\" > dump.cmd; sleep 1\"",origin_name,thicknezz/in_r, delta, refine_times;
    read "dump.cmd";
    print "" | "rm dump.cmd";
  }

/*
print "" | "echo dump \"sd_$(date +\"%s\"
echo dump \"sd_$(date +%Y_%m_%d_%H_%M_%S)\""\n"
"echo dump \\\"sd_$(date +\"%s\")\\\"\"\\n\"  "
"\"echo dump \\\\\\\"sd_$(date +\\\"%%s\\\")\\\\\\\"\\\"\\\\n\\\"  \"  "
*/

procedure converge_refine(integer rtm) {   // coverge to 0, after rtm times refine
   if rtm < 0 then {
      quit; 
      q;
   };
   local ini;
   if clock idiv (run_time*60) == ( submit_times - ost ) then {
           submit_times++;
           dump sprintf "%s_s%d.dmp",origin_name,submit_times;
           //print thicknezz | "sleep 260";
           //quit;
   };
   go 10;
   cc0;
   for(ini:=1;ini<=rtm;ini++) {
      //cc0;
      if count(vertex where z!=0, 1) == 0 then: {
        print "Saddle: ";
        saddle;
      };
      r;
      printf "current refine times: %g\n",refine_times;
      cc0;
   };
   if count(vertex where z!=0, 1) == 0 then: {
        print "Saddle: ";
        saddle;
        if count(vertex where z!=0, 1) > 0 then: {
           cc0;
        } else {
	   eigenprobe 0;
	   if eigenneg > 0 then: { 
	      j 0.00001; g 110; cc0;
	   };
	};
   };
   //cc0;
   dump_surf;
   save_dict;
}

function real get_step_t(real ten) {
   local num_ad;
   local ctension;
   local cthick;
   cthick := ten;
   num_ad := 0;
   while (cthick < 1.0) do {
      num_ad++;
      cthick := cthick * 10;
   };

   return 1/(10**num_ad);
}

/*
procedure converge_thick(real fh, real hs, integer rtm) {  //converge to thickness fh with step hs, hs can be positive or negtive
  local ch;
  //local gap;
  //gap = abs(thicknezz -fh);
  local eh;
  if hs > 0 then: {
    if thicknezz >= fh then:{
      printf "error: fh: %g is smaller than thicknezz: %g, but step is positive\n",fh, thicknezz;
      q;
    };
    eh := 1.01 * fh;
    for (ch := thicknezz; ch < eh; ch += hs) {
      set_thickness(ch);
      printf "Current thickness: %g\n",thicknezz;
      converge_refine(rtm - refine_times);
    };
  } else {
    if thicknezz <= fh then: {
       printf "error: fh: %g is larger than thicknezz: %g, but step is negtive\n",fh, thicknezz;
       q;
    };
    eh := 0.99 * fh;
    for (ch := thicknezz; ch > eh; ch += hs) {
      set_thickness(ch);
      printf "Current thickness: %g\n",thicknezz;
      converge_refine(rtm - refine_times);
    };
  };
}
*/

/*
procedure converge_tension(real tensionL, real tstep, integer rtm) {  //converge to tension tensionL with step tstep, tstep can be positive or negtive
  local ctension;
  //ctension := tensionG;
  //local cstep;
  //cstep := get_step_t(tensionG);
  if abs(tstep) > tensionG then:
      printf "Warning: tstep is larger than tensionG, tstep: %g, tensionG: %g",tstep, tensionG;
  local te;
  if tstep > 0 then: {
    if tensionL < tensionG then: {
      print "Error: tensionL is smaller than tensionG";
      q;
    };
    te := (1+0.01)*tensionL;
    for (ctension:=tensionG; ctension <= te; ctension+=tstep) {
      //cstep := get_step_t(ctension);
      //tensionG := ctension;
      //recalc;
      set_tension(ctension);
      printf "current tension: %g\n",tensionG;
      converge_refine(rtm - refine_times);
      
    };
  } else {
    if tensionL > tensionG then: {
      print "Error: tensionL is larger than tensionG";
      q;
    };
    te := (1-0.01)*tensionL;
    for (ctension:=tensionG; ctension >= te; ctension+=tstep) {
      //cstep := get_step_t(ctension);
      //tensionG := ctension;
      //recalc;
      set_tension(ctension);
      printf "current tension: %g\n",tensionG;
      converge_refine(rtm - refine_times);

    };
  };
    
  //printf "ctension: %g, tensionL: %g\n", ctension, tensionL;
  //if ctension == tensionL then:
  //  print "true";
}
*/

function real get_step(real thick) {
   local num_ad;
   local cthick;
   cthick := thick;
   num_ad := 0;
   while (cthick <= 1.0) do {
      num_ad++;
      cthick := cthick * 10;
   };
   
   return 1/(10**num_ad);
}

procedure converge_stage(real e_thick) {
   local c_thick; // current thick
   //local num_ad; // number after decimal
   local c_step;
   if(thicknezz < e_thick) then {
       printf "current thickness %g is larger than endthickness %g",thicknezz,e_thick;
       quit;
   };
   
   if thicknezz >= 1 then
	return;

   c_step := get_step(thicknezz);

   for(c_thick:=thicknezz; c_thick >= e_thick; c_thick-=c_step) {
	c_step := get_step(c_thick);
	set_thickness(c_thick);
	converge_refine(7-refine_times);
   };
}

test_verts := {
   local ct;
   ct := 0;

   foreach facets ff do {
      if ff.edges[1].vertex[1].id == ff.vertex[1].id && ff.edges[1].vertex[2].id == ff.vertex[2].id
         && ff.edges[2].vertex[1].id == ff.vertex[2].id && ff.edges[2].vertex[2].id == ff.vertex[3].id
         && ff.edges[3].vertex[1].id == ff.vertex[3].id && ff.edges[3].vertex[2].id == ff.vertex[1].id then
             ct++;
   };
   print ct;
   print facet_count;
}
