import json
import os
import configparser

my_path = os.getcwd()

def parseOptions():
    parser = configparser.ConfigParser()
    setPars= configparser.ConfigParser()
    pars = {}
    parser.read(os.path.join(my_path, "params", "racecar6.ini"))
        
    pars['ggv_file'] = json.loads(parser.get('GENERAL_OPTIONS', 'ggv_file'))
    pars["ax_max_machines_file"] = json.loads(parser.get('GENERAL_OPTIONS', 'ax_max_machines_file'))
    pars["stepsize_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'stepsize_opts'))
    pars["reg_smooth_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'reg_smooth_opts'))
    pars["veh_params"] = json.loads(parser.get('GENERAL_OPTIONS', 'veh_params'))
    pars["vel_calc_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'vel_calc_opts'))
    pars["optim_opts_shortest_path"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'optim_opts_shortest_path'))
    pars["optim_opts_mincurv"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'optim_opts_mincurv'))
    
    pars["curv_calc_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'curv_calc_opts'))
    pars["optim_opts_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'optim_opts_mintime'))
    pars["vehicle_params_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'vehicle_params_mintime'))
    pars["tire_params_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'tire_params_mintime'))
    pars["pwr_params_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'pwr_params_mintime'))
    pars["vehicle_params_mintime"]["wheelbase"] = (pars["vehicle_params_mintime"]["wheelbase_front"]
                                                       + pars["vehicle_params_mintime"]["wheelbase_rear"])
    
    for vals in pars:
        print(vals)
    try:
        choice1=input('Select Option to modify (str): ')
    except:
        print("Just enter something valid here")
        return False
    testBool = False
    try:
        print(len(pars[choice1].split(',')))
        testBool = True
    except:
        testBool = False
            
    if (testBool) :
        print("Current Setting")
        print(choice1 +' = ', pars[choice1])
        changeVal = input('Enter the modified input: ')
        pars[choice1] = changeVal
    else:
        print(choice1, " Menu:")
        print()
        for val in pars[choice1]:
            print (val)
        try:
            choice2=input('Select Option to modify from here (str): ')
            pars[choice1][choice2]
                    
        except:
            print("Just enter something valid here too")
            return False
        print("Current Setting")
        print(choice2 +' = ', pars[choice1][choice2])
        changeVal = input('Enter the modified input: ')
        pars[choice1][choice2] = changeVal
        print(pars[choice1][choice2])

    
    with open(os.path.join(my_path, "params", "racecar5.ini"), 'w') as f:
        setPars.write(f)
    
    f = open(os.path.join(my_path, "params", "racecar5.ini"), 'w')
    f.write("[GENERAL_OPTIONS]\n")
    f.write("ggv_file="+'\"'+pars["ggv_file"]+'\"'+'\n')
    f.write("ax_max_machines_file="+'\"'+pars["ax_max_machines_file"]+'\"'+'\n')
    f.write("stepsize_opts={\"stepsize_prep\": %.2f,\n" %pars["stepsize_opts"]["stepsize_prep"])
    f.write("\"stepsize_reg\": %.2f,\n" %pars["stepsize_opts"]["stepsize_reg"])
    f.write("\"stepsize_interp_after_opt\": %.2f}\n" %pars["stepsize_opts"]["stepsize_interp_after_opt"])
    f.write("reg_smooth_opts={\"k_reg\": %.2f,\n" %pars["reg_smooth_opts"]["k_reg"])
    f.write("\"s_reg\": %.2f}\n" %pars["reg_smooth_opts"]["s_reg"])
    f.write("curv_calc_opts={\"d_preview_curv\": %.2f,\n" %pars["curv_calc_opts"]["d_preview_curv"])
    f.write("\"d_review_curv\": %.2f,\n" %pars["curv_calc_opts"]["d_review_curv"])
    f.write("\"d_preview_head\": %.2f,\n" %pars["curv_calc_opts"]["d_preview_head"])
    f.write("\"d_review_head\": %.2f}\n" %pars["curv_calc_opts"]["d_review_head"])
    f.write("veh_params={\"v_max\": %.2f,\n" %pars["veh_params"]["v_max"])
    f.write("\"length\": %.2f,\n" %pars["veh_params"]["length"])
    f.write("\"width\": %.2f,\n" %pars["veh_params"]["width"])
    f.write("\"mass\": %.2f,\n" %pars["veh_params"]["mass"])
    f.write("\"dragcoeff\": %.2f,\n" %pars["veh_params"]["dragcoeff"])
    f.write("\"curvlim\": %.2f,\n" %pars["veh_params"]["curvlim"])
    f.write("\"g\": %.2f}\n" %pars["veh_params"]["g"])
    f.write("vel_calc_opts={\"dyn_model_exp\": %.2f,\n" %pars["vel_calc_opts"]["dyn_model_exp"])
    f.write("\"dyn_model_exp\": %.2f}\n" %pars["vel_calc_opts"]["dyn_model_exp"])
    f.write("[OPTIMIZATION_OPTIONS]\n")
    f.write("optim_opts_shortest_path={\"width_opt\": %.2f }\n" %int(pars["optim_opts_shortest_path"]["width_opt"]))
    f.write("optim_opts_mincurv={\"width_opt\": %.2f,\n" %int(pars["optim_opts_mincurv"]["width_opt"]))
    f.write("\"iqp_iters_min\": %.2f,\n" %int(pars["optim_opts_mincurv"]["iqp_iters_min"]))
    f.write("\"iqp_curverror_allowed\": %.2f}\n" %pars["optim_opts_mincurv"]["iqp_curverror_allowed"])
    
    f.write("optim_opts_mintime={\"width_opt\": %.2f,\n" %int(pars["optim_opts_mintime"]["width_opt"]))
    f.write("\"penalty_delta\": %.2f,\n" %int(pars["optim_opts_mintime"]["penalty_delta"]))
    f.write("\"penalty_F\": %.2f,\n" %pars["optim_opts_mintime"]["penalty_F"])
    f.write("\"mue\": %.2f,\n" %pars["optim_opts_mintime"]["mue"])
    f.write("\"n_gauss\": %.2f,\n" %pars["optim_opts_mintime"]["n_gauss"])
    f.write("\"dn\": %.2f,\n" %int(pars["optim_opts_mintime"]["dn"]))
    f.write("\"limit_energy\": " + str(pars["optim_opts_mintime"]["limit_energy"]) + ',\n')
    f.write("\"energy_limit\": %.2f,\n" %pars["optim_opts_mintime"]["energy_limit"])
    f.write("\"safe_traj\": " + str(pars["optim_opts_mintime"]["safe_traj"]) + ',\n')
    f.write("\"ax_pos_safe\": " + str(pars["optim_opts_mintime"]["ax_pos_safe"]) + ',\n')
    f.write("\"ax_neg_safe\": " + str(pars["optim_opts_mintime"]["ax_neg_safe"]) + ',\n')
    f.write("\"ay_safe\":" + str(pars["optim_opts_mintime"]["ay_safe"]) + ',\n')
    f.write("\"w_tr_reopt\": %.2f,\n" %pars["optim_opts_mintime"]["w_tr_reopt"])
    f.write("\"w_veh_reopt\": %.2f,\n" %pars["optim_opts_mintime"]["w_veh_reopt"])
    f.write("\"w_add_spl_regr\": %.2f,\n" %pars["optim_opts_mintime"]["w_add_spl_regr"])
    f.write("\"step_non_reg\": %.2f,\n" %pars["optim_opts_mintime"]["step_non_reg"])
    f.write("\"eps_kappa\": %.2f}\n" %pars["optim_opts_mintime"]["eps_kappa"])
    
    f.write("vehicle_params_mintime = {\"wheelbase_front\": %.2f,\n" %pars["vehicle_params_mintime"]["wheelbase_front"])
    f.write("\"wheelbase_rear\": %.2f,\n" %pars["vehicle_params_mintime"]["wheelbase_rear"])
    f.write("\"track_width_front\": %.2f,\n" %pars["vehicle_params_mintime"]["track_width_front"])
    f.write("\"track_width_rear\": %.2f,\n" %pars["vehicle_params_mintime"]["track_width_rear"])
    f.write("\"cog_z\": %.2f,\n" %pars["vehicle_params_mintime"]["cog_z"])
    f.write("\"I_z\": %.2f,\n" %pars["vehicle_params_mintime"]["I_z"])
    f.write("\"liftcoeff_front\": %.2f,\n" %pars["vehicle_params_mintime"]["liftcoeff_front"])
    f.write("\"liftcoeff_rear\": %.2f,\n" %pars["vehicle_params_mintime"]["liftcoeff_rear"])
    f.write("\"k_brake_front\": %.2f,\n" %pars["vehicle_params_mintime"]["k_brake_front"])
    f.write("\"k_drive_front\": %.2f,\n" %pars["vehicle_params_mintime"]["k_drive_front"])
    f.write("\"k_roll\": %.2f,\n" %pars["vehicle_params_mintime"]["k_roll"])
    f.write("\"t_delta\": %.2f,\n" %pars["vehicle_params_mintime"]["t_delta"])
    f.write("\"t_drive\": %.2f,\n" %pars["vehicle_params_mintime"]["t_drive"])
    f.write("\"t_brake\": %.2f,\n" %pars["vehicle_params_mintime"]["t_brake"])
    f.write("\"power_max\": %.2f,\n" %pars["vehicle_params_mintime"]["power_max"])
    f.write("\"f_drive_max\": %.2f,\n" %pars["vehicle_params_mintime"]["f_drive_max"])
    f.write("\"f_brake_max\": %.2f,\n" %pars["vehicle_params_mintime"]["f_brake_max"])
    f.write("\"delta_max\": %.2f}\n" %pars["vehicle_params_mintime"]["delta_max"])
    
    f.write("tire_params_mintime = {\"c_roll\": %.3f,\n" %pars["tire_params_mintime"]["c_roll"])
    f.write("\"f_z0\": %.2f,\n" %pars["tire_params_mintime"]["f_z0"])
    f.write("\"B_front\": %.2f,\n" %pars["tire_params_mintime"]["B_front"])
    f.write("\"C_front\": %.2f,\n" %pars["tire_params_mintime"]["C_front"])
    f.write("\"eps_front\": %.2f,\n" %pars["tire_params_mintime"]["eps_front"])
    f.write("\"E_front\": %.2f,\n" %pars["tire_params_mintime"]["E_front"])
    f.write("\"B_rear\": %.2f,\n" %pars["tire_params_mintime"]["B_rear"])
    f.write("\"C_rear\": %.2f,\n" %pars["tire_params_mintime"]["C_rear"])
    f.write("\"eps_rear\": %.2f,\n" %pars["tire_params_mintime"]["eps_rear"])
    f.write("\"E_rear\": %.2f}\n" %pars["tire_params_mintime"]["E_rear"])
    
    f.write("pwr_params_mintime = {\"pwr_behavior\": " + str(pars["pwr_params_mintime"]["pwr_behavior"]) +'\n')
    f.write("\"simple_loss\": " + str(pars["pwr_params_mintime"]["simple_loss"]) + '\n')
    f.write("\"T_env\": %.2f,\n" %pars["pwr_params_mintime"]["T_env"])
    f.write("\"T_mot_ini\": %.2f,\n" %pars["pwr_params_mintime"]["T_mot_ini"])
    f.write("\"T_batt_ini\": %.2f,\n" %pars["pwr_params_mintime"]["T_batt_ini"])
    f.write("\"T_inv_ini\": %.2f,\n" %pars["pwr_params_mintime"]["T_inv_ini"])
    f.write("\"T_cool_mi_ini\": %.2f,\n" %pars["pwr_params_mintime"]["T_cool_mi_ini"])
    f.write("\"T_cool_b_ini\": %.2f,\n" %pars["pwr_params_mintime"]["T_cool_b_ini"])
    f.write("\"r_wheel\": %.2f,\n" %pars["pwr_params_mintime"]["r_wheel"])
    f.write("\"R_i_sumo\": %.2f,\n" %pars["pwr_params_mintime"]["R_i_sumo"])
    f.write("\"R_i_simple\": %.2f,\n" %pars["pwr_params_mintime"]["R_i_simple"])
    f.write("\"R_i_offset\": %.2f,\n" %pars["pwr_params_mintime"]["R_i_offset"])
    f.write("\"R_i_slope\": %.2f,\n" %pars["pwr_params_mintime"]["R_i_slope"])
    f.write("\"V_OC_simple\": %.2f,\n" %pars["pwr_params_mintime"]["V_OC_simple"])
    f.write("\"SOC_ini\": %.2f,\n" %pars["pwr_params_mintime"]["SOC_ini"])
    f.write("\"C_batt\": %.2f,\n" %pars["pwr_params_mintime"]["C_batt"])
    f.write("\"N_cells_serial\": %.2f,\n" %pars["pwr_params_mintime"]["N_cells_serial"])
    f.write("\"N_cells_parallel\": %.2f,\n" %pars["pwr_params_mintime"]["N_cells_parallel"])
    f.write("\"temp_mot_max\": %.2f,\n" %pars["pwr_params_mintime"]["temp_mot_max"])
    f.write("\"temp_batt_max\": %.2f,\n" %pars["pwr_params_mintime"]["temp_batt_max"])
    f.write("\"temp_inv_max\": %.2f,\n" %pars["pwr_params_mintime"]["temp_inv_max"])
    f.write("\"N_machines\": %.2f,\n" %pars["pwr_params_mintime"]["N_machines"])
    f.write("\"transmission\": %.2f,\n" %pars["pwr_params_mintime"]["transmission"])
    f.write("\"MotorConstant\": %.2f,\n" %pars["pwr_params_mintime"]["MotorConstant"])
    f.write("\"C_therm_machine\": %.2f,\n" %pars["pwr_params_mintime"]["C_therm_machine"])
    f.write("\"C_therm_inv\": %.2f,\n" %pars["pwr_params_mintime"]["C_therm_inv"])
    f.write("\"C_therm_cell\": %.2f,\n" %pars["pwr_params_mintime"]["C_therm_cell"])
    f.write("\"C_TempCopper\": %.2f,\n" %pars["pwr_params_mintime"]["C_TempCopper"])
    f.write("\"m_therm_fluid_mi\": %.2f,\n" %pars["pwr_params_mintime"]["m_therm_fluid_mi"])
    f.write("\"m_therm_fluid_b\": %.2f,\n" %pars["pwr_params_mintime"]["m_therm_fluid_b"])
    f.write("\"R_Phase\": %.2f,\n" %pars["pwr_params_mintime"]["R_Phase"])
    f.write("\"r_rotor_int\": %.2f,\n" %pars["pwr_params_mintime"]["r_rotor_int"])
    f.write("\"r_rotor_ext\": %.2f,\n" %pars["pwr_params_mintime"]["r_rotor_ext"])
    f.write("\"r_stator_int\": %.2f,\n" %pars["pwr_params_mintime"]["r_stator_int"])
    f.write("\"r_stator_ext\": %.2f,\n" %pars["pwr_params_mintime"]["r_stator_ext"])
    f.write("\"l_machine\": %.2f,\n" %pars["pwr_params_mintime"]["l_machine"])
    f.write("\"A_cool_inflate_machine\": %.2f,\n" %pars["pwr_params_mintime"]["A_cool_inflate_machine"])
    f.write("\"A_cool_inv\": %.2f,\n" %pars["pwr_params_mintime"]["A_cool_inv"])
    f.write("\"A_cool_rad\": %.2f,\n" %pars["pwr_params_mintime"]["A_cool_rad"])
    f.write("\"k_iro\": %.2f,\n" %pars["pwr_params_mintime"]["k_iro"])
    f.write("\"h_air\": %.2f,\n" %pars["pwr_params_mintime"]["h_air"])
    f.write("\"h_air_gap\": %.2f,\n" %pars["pwr_params_mintime"]["h_air_gap"])
    f.write("\"h_fluid_mi\": %.2f,\n" %pars["pwr_params_mintime"]["h_fluid_mi"])
    f.write("\"c_heat_fluid\": %.2f,\n" %pars["pwr_params_mintime"]["c_heat_fluid"])
    f.write("\"flow_rate_inv\": %.2f,\n" %pars["pwr_params_mintime"]["flow_rate_inv"])
    f.write("\"flow_rate_rad\": %.2f,\n" %pars["pwr_params_mintime"]["flow_rate_rad"])
    f.write("\"machine_simple_a\": %.2f,\n" %pars["pwr_params_mintime"]["machine_simple_a"])
    f.write("\"machine_simple_b\": %.2f,\n" %pars["pwr_params_mintime"]["machine_simple_b"])
    f.write("\"machine_simple_c\": %.2f,\n" %pars["pwr_params_mintime"]["machine_simple_c"])
    f.write("\"V_ref\": %.2f,\n" %pars["pwr_params_mintime"]["V_ref"])
    f.write("\"I_ref\": %.2f,\n" %pars["pwr_params_mintime"]["I_ref"])
    f.write("\"V_ce_offset\": %.2f,\n" %pars["pwr_params_mintime"]["V_ce_offset"])
    f.write("\"V_ce_slope\": %.2f,\n" %pars["pwr_params_mintime"]["V_ce_slope"])
    f.write("\"E_on\": %.2f,\n" %pars["pwr_params_mintime"]["E_on"])
    f.write("\"E_off\": %.2f,\n" %pars["pwr_params_mintime"]["E_off"])
    f.write("\"E_rr\": %.2f,\n" %pars["pwr_params_mintime"]["E_rr"])
    f.write("\"f_sw\": %.2f,\n" %pars["pwr_params_mintime"]["f_sw"])
    f.write("\"inverter_simple_a\": %.2f,\n" %pars["pwr_params_mintime"]["inverter_simple_a"])
    f.write("\"inverter_simple_b\": %.2f,\n" %pars["pwr_params_mintime"]["inverter_simple_b"])
    f.write("\"inverter_simple_c\": %.2f}\n" %pars["pwr_params_mintime"]["inverter_simple_c"])
    
    f.close()
    
    return True

if __name__=='__main__':
    pass
