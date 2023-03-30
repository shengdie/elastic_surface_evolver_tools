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

procedure set_delta(real dtt) {
    delta := dtt;
    dratio := 1+delta;
    foreach vertex vv do {
        if max(vv.edge, inborder) > 0 then {
            vv.x := vv.ref_coord[1] * dratio;
            vv.y := vv.ref_coord[2] * dratio;
            set vv fixed;
        };
    };
    if in_r == 0 then {
        c_r := dratio * out_r;
    } else {
        c_r := dratio * in_r;
    };
    
    recalc;
}

procedure set_thickness(real th) {
    local pr;

    thicknezz := th;
    // Assume all facets have the same Poisson ratio
    pr := facets[1].poisson_ratio;
    bend.modulus := thicknezz**2 / 6 / (1 - pr**2);
    if is_defined("gbend") then exec "gbend.modulus := -thicknezz**2 /12/(1 + facets[1].poisson_ratio);";
    //gbend.modulus := -bend.modulus * (1 - pr)/2;
    //gbend.modulus := -thicknezz**2 /12/(1 + pr);
    //edgewidth := 1/3;  // Right triangles along edge
    //edgewidth := sqrt(3)/4;  // Equilateral triangles along edge
    //freeedgebend.modulus := bend.modulus/4 * (1 - pr**2) * edgewidth * gridsize;
    //set_ksub(minnum, maxnum);

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
    local temp_r;

    set vertex old_vid id;
    set edge old_eid id;

    refine edges;
    set edge divide 0;

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
        //print "set refine";
        if max(vv.edge, inborder) > 0 then {
            temp_r := in_r / sqrt(vv.ref_coord[1]^2 + vv.ref_coord[2]^2);
            vv.ref_coord[1] := temp_r * vv.ref_coord[1];
            vv.ref_coord[2] := temp_r * vv.ref_coord[2];
            vv.x := vv.ref_coord[1] * dratio;
            vv.y := vv.ref_coord[2] * dratio;
            set vv fixed;
            // foreach vv.edge.vertex vs where vs.id != vv.id do {
            //     if vs.v_constraint_list[1] == 1 then set vv constraint (vs.v_constraint_list[2] imod 0x100000);
            // };
        } else if max(vv.edge, outborder) > 0 then {
            temp_r := out_r / sqrt(vv.ref_coord[1]^2 + vv.ref_coord[2]^2);
            vv.ref_coord[1] := temp_r * vv.ref_coord[1];
            vv.ref_coord[2] := temp_r * vv.ref_coord[2];
        } else set vv gbend;
    };
    //if is_defined("gbend") then exec "foreach vertex vv where old_vid == 0 do {if max(vv.edge, inborder) == 0 && max(vv.edge, outborder) == 0 then set vv gbend}";

    do {
        flush_counts;
        equiangulate edge where rheight or old_eid==0;
    } while equi_count;

    set_form_factors;
    // Edge correction needs to be updated:
    set_thickness(thicknezz);
}

/*
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
*/

r :::= {
    //set edge divide 1 where length > (1 + 2/sqrt(3))/2 * gridsize;
    refinemarked();
    gridsize := gridsize/2;
    //if is_defined("fa") then exec "fa := fa / 4";
    recalc;
    refine_times++;
}
