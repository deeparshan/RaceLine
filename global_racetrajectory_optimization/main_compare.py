import opt_mintime_traj
import numpy as np
import time
import json
import os
import trajectory_planning_helpers as tph
import copy
import matplotlib.pyplot as plt
import configparser
import pkg_resources
import helper_funcs_glob

"""
Created by:
Alexander Heilmeier

Documentation:
This script has to be executed to generate an optimal trajectory based on a given reference track.
"""
def main_compare():    
    # ----------------------------------------------------------------------------------------------------------------------
    # USER INPUT -----------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # choose vehicle parameter file ----------------------------------------------------------------------------------------
    file_paths1 = {"veh_params_file": "racecar1.ini"}
    file_paths2 = {"veh_params_file": "racecar2.ini"}
    
    # debug and plot options -----------------------------------------------------------------------------------------------
    debug = True                                    # print console messages
    plot_opts = {"mincurv_curv_lin": True,         # plot curv. linearization (original and solution based) (mincurv only)
                 "raceline": True,                  # plot optimized path
                 "imported_bounds": True,          # plot imported bounds (analyze difference to interpolated bounds)
                 "raceline_curv": True,             # plot curvature profile of optimized path
                 "racetraj_vel": True,              # plot velocity profile
                 "racetraj_vel_3d": True,          # plot 3D velocity profile above raceline
                 "racetraj_vel_3d_stepsize": 1.0,   # [m] vertical lines stepsize in 3D velocity profile plot
                 "spline_normals": True,           # plot spline normals to check for crossings
                 "mintime_plots": True}            # plot states, controls, friction coeffs etc. (mintime only)
    
    # select track file (including centerline coordinates + track widths) --------------------------------------------------
    # file_paths["track_name"] = "rounded_rectangle"                              # artificial track
    # file_paths["track_name"] = "handling_track"                                 # artificial track
    #file_paths["track_name"] = "berlin_2018"                                    # Berlin Formula E 2018
    file_paths1["track_name"] = "modena_2019"                                    # Modena 2019
    file_paths2["track_name"] = "modena_2019"
    # set import options ---------------------------------------------------------------------------------------------------
    imp_opts = {"flip_imp_track": False,                # flip imported track to reverse direction
                "set_new_start": False,                 # set new starting point (changes order, not coordinates)
                "new_start": np.array([0.0, -47.0]),    # [x_m, y_m]
                "min_track_width": None,                # [m] minimum enforced track width (set None to deactivate)
                "num_laps": 1}                          # number of laps to be driven (significant with powertrain-option),
                                                        # only relevant in mintime-optimization
    
    # set optimization type ------------------------------------------------------------------------------------------------
    # 'shortest_path'       shortest path optimization
    # 'mincurv'             minimum curvature optimization without iterative call
    # 'mincurv_iqp'         minimum curvature optimization with iterative call
    # 'mintime'             time-optimal trajectory optimization
    opt_type = 'mintime'
    
    # set mintime specific options (mintime only) --------------------------------------------------------------------------
    # tpadata:                      set individual friction map data file if desired (e.g. for varmue maps), else set None,
    #                               e.g. "berlin_2018_varmue08-12_tpadata.json"
    # warm_start:                   [True/False] warm start IPOPT if previous result is available for current track
    # var_friction:                 [-] None, "linear", "gauss" -> set if variable friction coefficients should be used
    #                               either with linear regression or with gaussian basis functions (requires friction map)
    # reopt_mintime_solution:       reoptimization of the mintime solution by min. curv. opt. for improved curv. smoothness
    # recalc_vel_profile_by_tph:    override mintime velocity profile by ggv based calculation (see TPH package)
    
    mintime_opts = {"tpadata": None,
                    "warm_start": False,
                    "var_friction": None,
                    "reopt_mintime_solution": False,
                    "recalc_vel_profile_by_tph": False}
    
    # lap time calculation table -------------------------------------------------------------------------------------------
    lap_time_mat_opts = {"use_lap_time_mat": False,             # calculate a lap time matrix (diff. top speeds and scales)
                         "gg_scale_range": [0.3, 1.0],          # range of gg scales to be covered
                         "gg_scale_stepsize": 0.05,             # step size to be applied
                         "top_speed_range": [100.0, 150.0],     # range of top speeds to be simulated [in km/h]
                         "top_speed_stepsize": 5.0,             # step size to be applied
                         "file": "lap_time_matrix.csv"}         # file name of the lap time matrix (stored in "outputs")
    
    # ----------------------------------------------------------------------------------------------------------------------
    # CHECK USER INPUT -----------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    if opt_type not in ["shortest_path", "mincurv", "mincurv_iqp", "mintime"]:
        raise IOError("Unknown optimization type!")
    
    if opt_type == "mintime" and not mintime_opts["recalc_vel_profile_by_tph"] and lap_time_mat_opts["use_lap_time_mat"]:
        raise IOError("Lap time calculation table should be created but velocity profile recalculation with TPH solver is"
                      " not allowed!")
    
    # ----------------------------------------------------------------------------------------------------------------------
    # CHECK PYTHON DEPENDENCIES --------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # get current path
    file_paths1["module"] = os.path.dirname(os.path.abspath(__file__))
    file_paths2["module"] = os.path.dirname(os.path.abspath(__file__))
    
    # read dependencies from requirements.txt
    requirements_path = os.path.join(file_paths1["module"], 'requirements.txt')
    requirements_path = os.path.join(file_paths2["module"], 'requirements.txt')
    dependencies = []
    
    with open(requirements_path, 'r') as fh:
        line = fh.readline()
    
        while line:
            dependencies.append(line.rstrip())
            line = fh.readline()
    
    # check dependencies
    pkg_resources.require(dependencies)
    
    # ----------------------------------------------------------------------------------------------------------------------
    # INITIALIZATION OF PATHS ----------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # assemble track import path
    file_paths1["track_file"] = os.path.join(file_paths1["module"], "inputs", "tracks", file_paths1["track_name"] + ".csv")
    file_paths2["track_file"] = os.path.join(file_paths2["module"], "inputs", "tracks", file_paths2["track_name"] + ".csv")
    
    
    # assemble friction map import paths
    file_paths1["tpamap"] = os.path.join(file_paths1["module"], "inputs", "frictionmaps",
                                        file_paths1["track_name"] + "_tpamap.csv")
    file_paths2["tpamap"] = os.path.join(file_paths2["module"], "inputs", "frictionmaps",
                                        file_paths2["track_name"] + "_tpamap.csv")
    
    if mintime_opts["tpadata"] is None:
        file_paths1["tpadata"] = os.path.join(file_paths1["module"], "inputs", "frictionmaps",
                                             file_paths1["track_name"] + "_tpadata.json")
        file_paths2["tpadata"] = os.path.join(file_paths2["module"], "inputs", "frictionmaps",
                                             file_paths2["track_name"] + "_tpadata.json")
    else:
        file_paths1["tpadata"] = os.path.join(file_paths1["module"], "inputs", "frictionmaps", mintime_opts["tpadata"])
        file_paths2["tpadata"] = os.path.join(file_paths2["module"], "inputs", "frictionmaps", mintime_opts["tpadata"])
    
    # check if friction map files are existing if the var_friction option was set
    if opt_type == 'mintime' \
            and mintime_opts["var_friction"] is not None \
            and not (os.path.exists(file_paths1["tpadata"]) and os.path.exists(file_paths1["tpamap"])):
    
        mintime_opts["var_friction"] = None
        print("WARNING: var_friction option is not None but friction map data is missing for current track -> Setting"
              " var_friction to None!")
    
    # create outputs folder(s)
    os.makedirs(file_paths1["module"] + "/outputs", exist_ok=True)
    os.makedirs(file_paths2["module"] + "/outputs", exist_ok=True)
    
    if opt_type == 'mintime':
        os.makedirs(file_paths1["module"] + "/outputs/mintime", exist_ok=True)
        os.makedirs(file_paths2["module"] + "/outputs/mintime", exist_ok=True)
    
    # assemble export paths
    file_paths1["mintime_export"] = os.path.join(file_paths1["module"], "outputs", "mintime")
    file_paths1["traj_race_export"] = os.path.join(file_paths1["module"], "outputs", "traj_race_cl.csv")
    # file_paths["traj_ltpl_export"] = os.path.join(file_paths["module"], "outputs", "traj_ltpl_cl.csv")
    file_paths1["lap_time_mat_export"] = os.path.join(file_paths1["module"], "outputs", lap_time_mat_opts["file"])
    
    file_paths2["mintime_export"] = os.path.join(file_paths1["module"], "outputs", "mintime")
    file_paths2["traj_race_export"] = os.path.join(file_paths1["module"], "outputs", "traj_race_cl.csv")
    # file_paths["traj_ltpl_export"] = os.path.join(file_paths["module"], "outputs", "traj_ltpl_cl.csv")
    file_paths2["lap_time_mat_export"] = os.path.join(file_paths1["module"], "outputs", lap_time_mat_opts["file"])
    
    # ----------------------------------------------------------------------------------------------------------------------
    # IMPORT VEHICLE DEPENDENT PARAMETERS ----------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # load vehicle 1 parameter file into a "pars1" dict
    parser = configparser.ConfigParser()
    pars1 = {}
    
    if not parser.read(os.path.join(file_paths1["module"], "params", file_paths1["veh_params_file"])):
        raise ValueError('Specified config file does not exist or is empty!')
    
    pars1["ggv_file"] = json.loads(parser.get('GENERAL_OPTIONS', 'ggv_file'))
    pars1["ax_max_machines_file"] = json.loads(parser.get('GENERAL_OPTIONS', 'ax_max_machines_file'))
    pars1["stepsize_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'stepsize_opts'))
    pars1["reg_smooth_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'reg_smooth_opts'))
    pars1["veh_params"] = json.loads(parser.get('GENERAL_OPTIONS', 'veh_params'))
    pars1["vel_calc_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'vel_calc_opts'))
    
    if opt_type == 'shortest_path':
        pars1["optim_opts"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'optim_opts_shortest_path'))
    
    elif opt_type in ['mincurv', 'mincurv_iqp']:
        pars1["optim_opts"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'optim_opts_mincurv'))
    
    elif opt_type == 'mintime':
        pars1["curv_calc_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'curv_calc_opts'))
        pars1["optim_opts"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'optim_opts_mintime'))
        pars1["vehicle_params_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'vehicle_params_mintime'))
        pars1["tire_params_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'tire_params_mintime'))
        pars1["pwr_params_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'pwr_params_mintime'))
    
        # modification of mintime options/parameters
        pars1["optim_opts"]["var_friction"] = mintime_opts["var_friction"]
        pars1["optim_opts"]["warm_start"] = mintime_opts["warm_start"]
        pars1["vehicle_params_mintime"]["wheelbase"] = (pars1["vehicle_params_mintime"]["wheelbase_front"]
                                                       + pars1["vehicle_params_mintime"]["wheelbase_rear"])
    
    # set import path for ggv diagram and ax_max_machines (if required)
    if not (opt_type == 'mintime' and not mintime_opts["recalc_vel_profile_by_tph"]):
        file_paths1["ggv_file"] = os.path.join(file_paths1["module"], "inputs", "veh_dyn_info", pars1["ggv_file"])
        file_paths1["ax_max_machines_file"] = os.path.join(file_paths1["module"], "inputs", "veh_dyn_info",
                                                          pars1["ax_max_machines_file"])
    
    #
    # load vehicle 2 parameter file into a "pars2" dict
    parser = configparser.ConfigParser()
    pars2 = {}
    
    if not parser.read(os.path.join(file_paths2["module"], "params", file_paths2["veh_params_file"])):
        raise ValueError('Specified config file does not exist or is empty!')
    
    pars2["ggv_file"] = json.loads(parser.get('GENERAL_OPTIONS', 'ggv_file'))
    pars2["ax_max_machines_file"] = json.loads(parser.get('GENERAL_OPTIONS', 'ax_max_machines_file'))
    pars2["stepsize_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'stepsize_opts'))
    pars2["reg_smooth_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'reg_smooth_opts'))
    pars2["veh_params"] = json.loads(parser.get('GENERAL_OPTIONS', 'veh_params'))
    pars2["vel_calc_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'vel_calc_opts'))
    
    if opt_type == 'shortest_path':
        pars2["optim_opts"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'optim_opts_shortest_path'))
    
    elif opt_type in ['mincurv', 'mincurv_iqp']:
        pars2["optim_opts"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'optim_opts_mincurv'))
    
    elif opt_type == 'mintime':
        pars2["curv_calc_opts"] = json.loads(parser.get('GENERAL_OPTIONS', 'curv_calc_opts'))
        pars2["optim_opts"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'optim_opts_mintime'))
        pars2["vehicle_params_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'vehicle_params_mintime'))
        pars2["tire_params_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'tire_params_mintime'))
        pars2["pwr_params_mintime"] = json.loads(parser.get('OPTIMIZATION_OPTIONS', 'pwr_params_mintime'))
    
        # modification of mintime options/parameters
        pars2["optim_opts"]["var_friction"] = mintime_opts["var_friction"]
        pars2["optim_opts"]["warm_start"] = mintime_opts["warm_start"]
        pars2["vehicle_params_mintime"]["wheelbase"] = (pars2["vehicle_params_mintime"]["wheelbase_front"]
                                                       + pars2["vehicle_params_mintime"]["wheelbase_rear"])
    
    # set import path for ggv diagram and ax_max_machines (if required)
    if not (opt_type == 'mintime' and not mintime_opts["recalc_vel_profile_by_tph"]):
        file_paths2["ggv_file"] = os.path.join(file_paths2["module"], "inputs", "veh_dyn_info", pars2["ggv_file"])
        file_paths2["ax_max_machines_file"] = os.path.join(file_paths2["module"], "inputs", "veh_dyn_info",
                                                          pars2["ax_max_machines_file"])
    
    # ----------------------------------------------------------------------------------------------------------------------
    # IMPORT TRACK AND VEHICLE 1 DYNAMICS INFORMATION ------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # save start time
    t_start = time.perf_counter()
    
    # import track
    reftrack_imp1 = helper_funcs_glob.src.import_track.import_track(imp_opts=imp_opts,
                                                                   file_path=file_paths1["track_file"],
                                                                   width_veh=pars1["veh_params"]["width"])
    
    # import ggv and ax_max_machines (if required)
    if not (opt_type == 'mintime' and not mintime_opts["recalc_vel_profile_by_tph"]):
        ggv, ax_max_machines = tph.import_veh_dyn_info.\
            import_veh_dyn_info(ggv_import_path=file_paths1["ggv_file"],
                                ax_max_machines_import_path=file_paths1["ax_max_machines_file"])
    else:
        ggv = None
        ax_max_machines = None
    
    # set ax_pos_safe / ax_neg_safe / ay_safe if required and not set in parameters file
    if opt_type == 'mintime' and pars1["optim_opts"]["safe_traj"] \
            and (pars1["optim_opts"]["ax_pos_safe"] is None
                 or pars1["optim_opts"]["ax_neg_safe"] is None
                 or pars1["optim_opts"]["ay_safe"] is None):
    
        # get ggv if not available
        if ggv is None:
            ggv = tph.import_veh_dyn_info. \
                import_veh_dyn_info(ggv_import_path=file_paths1["ggv_file"],
                                    ax_max_machines_import_path=file_paths1["ax_max_machines_file"])[0]
    
        # limit accelerations
        if pars1["optim_opts"]["ax_pos_safe"] is None:
            pars1["optim_opts"]["ax_pos_safe"] = np.amin(ggv[:, 1])
        if pars1["optim_opts"]["ax_neg_safe"] is None:
            pars1["optim_opts"]["ax_neg_safe"] = -np.amin(ggv[:, 1])
        if pars1["optim_opts"]["ay_safe"] is None:
            pars1["optim_opts"]["ay_safe"] = np.amin(ggv[:, 2])
    
    #
    # ----------------------------------------------------------------------------------------------------------------------
    # IMPORT TRACK AND VEHICLE 2 DYNAMICS INFORMATION ------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # save start time
    t_start = time.perf_counter()
    
    # import track
    reftrack_imp2 = helper_funcs_glob.src.import_track.import_track(imp_opts=imp_opts,
                                                                   file_path=file_paths2["track_file"],
                                                                   width_veh=pars2["veh_params"]["width"])
    
    # import ggv and ax_max_machines (if required)
    if not (opt_type == 'mintime' and not mintime_opts["recalc_vel_profile_by_tph"]):
        ggv, ax_max_machines = tph.import_veh_dyn_info.\
            import_veh_dyn_info(ggv_import_path=file_paths2["ggv_file"],
                                ax_max_machines_import_path=file_paths2["ax_max_machines_file"])
    else:
        ggv = None
        ax_max_machines = None
    
    # set ax_pos_safe / ax_neg_safe / ay_safe if required and not set in parameters file
    if opt_type == 'mintime' and pars2["optim_opts"]["safe_traj"] \
            and (pars2["optim_opts"]["ax_pos_safe"] is None
                 or pars2["optim_opts"]["ax_neg_safe"] is None
                 or pars2["optim_opts"]["ay_safe"] is None):
    
        # get ggv if not available
        if ggv is None:
            ggv = tph.import_veh_dyn_info. \
                import_veh_dyn_info(ggv_import_path=file_paths2["ggv_file"],
                                    ax_max_machines_import_path=file_paths2["ax_max_machines_file"])[0]
    
        # limit accelerations
        if pars2["optim_opts"]["ax_pos_safe"] is None:
            pars2["optim_opts"]["ax_pos_safe"] = np.amin(ggv[:, 1])
        if pars2["optim_opts"]["ax_neg_safe"] is None:
            pars2["optim_opts"]["ax_neg_safe"] = -np.amin(ggv[:, 1])
        if pars2["optim_opts"]["ay_safe"] is None:
            pars2["optim_opts"]["ay_safe"] = np.amin(ggv[:, 2])
    
    # ----------------------------------------------------------------------------------------------------------------------
    # PREPARE REFTRACK 1 -----------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    reftrack_interp1, normvec_normalized_interp1, a_interp1, coeffs_x_interp1, coeffs_y_interp1 = \
        helper_funcs_glob.src.prep_track.prep_track(reftrack_imp=reftrack_imp1,
                                                    reg_smooth_opts=pars1["reg_smooth_opts"],
                                                    stepsize_opts=pars1["stepsize_opts"],
                                                    debug=debug,
                                                    min_width=imp_opts["min_track_width"])
    
    #
    # ----------------------------------------------------------------------------------------------------------------------
    # PREPARE REFTRACK 2 -----------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    reftrack_interp2, normvec_normalized_interp2, a_interp2, coeffs_x_interp2, coeffs_y_interp2 = \
        helper_funcs_glob.src.prep_track.prep_track(reftrack_imp=reftrack_imp2,
                                                    reg_smooth_opts=pars2["reg_smooth_opts"],
                                                    stepsize_opts=pars2["stepsize_opts"],
                                                    debug=debug,
                                                    min_width=imp_opts["min_track_width"])
    
    # ----------------------------------------------------------------------------------------------------------------------
    # CALL OPTIMIZATION 1 ----------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # if reoptimization of mintime solution is used afterwards we have to consider some additional deviation in the first
    # optimization
    if opt_type == 'mintime' and mintime_opts["reopt_mintime_solution"]:
        w_veh_tmp = pars1["optim_opts"]["width_opt"] + (pars1["optim_opts"]["w_tr_reopt"] - pars1["optim_opts"]["w_veh_reopt"])
        w_veh_tmp += pars1["optim_opts"]["w_add_spl_regr"]
        pars_tmp1 = copy.deepcopy(pars1)
        pars_tmp1["optim_opts"]["width_opt"] = w_veh_tmp
    else:
        pars_tmp1 = pars1
        pars_tmp2 = pars2
    
    # call optimization
    if opt_type == 'mincurv':
        alpha_opt1 = tph.opt_min_curv.opt_min_curv(reftrack=reftrack_interp1,
                                                  normvectors=normvec_normalized_interp1,
                                                  A=a_interp1,
                                                  kappa_bound=pars1["veh_params"]["curvlim"],
                                                  w_veh=pars1["optim_opts"]["width_opt"],
                                                  print_debug=debug,
                                                  plot_debug=plot_opts["mincurv_curv_lin"])[0]
    
    elif opt_type == 'mincurv_iqp':
        alpha_opt1, reftrack_interp1, normvec_normalized_interp1 = tph.iqp_handler.\
            iqp_handler(reftrack=reftrack_interp1,
                        normvectors=normvec_normalized_interp1,
                        A=a_interp1,
                        kappa_bound=pars1["veh_params"]["curvlim"],
                        w_veh=pars1["optim_opts"]["width_opt"],
                        print_debug=debug,
                        plot_debug=plot_opts["mincurv_curv_lin"],
                        stepsize_interp=pars1["stepsize_opts"]["stepsize_reg"],
                        iters_min=pars1["optim_opts"]["iqp_iters_min"],
                        curv_error_allowed=pars1["optim_opts"]["iqp_curverror_allowed"])
    
    elif opt_type == 'shortest_path':
        alpha_opt1 = tph.opt_shortest_path.opt_shortest_path(reftrack=reftrack_interp1,
                                                            normvectors=normvec_normalized_interp1,
                                                            w_veh=pars1["optim_opts"]["width_opt"],
                                                            print_debug=debug)
    
    elif opt_type == 'mintime':
        # reftrack_interp, a_interp and normvec_normalized_interp are returned for the case that non-regular sampling was
        # applied
        alpha_opt1, v_opt1, reftrack_interp1, a_interp_tmp1, normvec_normalized_interp1, alpha_opt2, v_opt2, reftrack_interp2, a_interp_tmp2, normvec_normalized_interp2 = opt_mintime_traj.src.opt_compare.\
            opt_compare(reftrack1=reftrack_interp1,
                        coeffs_x1=coeffs_x_interp1,
                        coeffs_y1=coeffs_y_interp1,
                        normvectors1=normvec_normalized_interp1,
                        pars1=pars_tmp1,
                        tpamap_path1=file_paths1["tpamap"],
                        tpadata_path1=file_paths1["tpadata"],
                        export_path1=file_paths1["mintime_export"],
                        reftrack2=reftrack_interp2,
                        coeffs_x2=coeffs_x_interp2,
                        coeffs_y2=coeffs_y_interp2,
                        normvectors2=normvec_normalized_interp2,
                        pars2=pars_tmp2,
                        tpamap_path2=file_paths2["tpamap"],
                        tpadata_path2=file_paths2["tpadata"],
                        export_path2=file_paths2["mintime_export"],
                        print_debug=debug,
                        plot_debug=plot_opts["mintime_plots"])
    
        # replace a_interp if necessary
        if a_interp_tmp1 is not None:
            a_interp1 = a_interp_tmp1
        if a_interp_tmp2 is not None:
            a_interp2 = a_interp_tmp1
    
    else:
        raise ValueError('Unknown optimization type!')
        # alpha_opt = np.zeros(reftrack_interp.shape[0])
    
    # ----------------------------------------------------------------------------------------------------------------------
    # CALL OPTIMIZATION 2 ----------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # if reoptimization of mintime solution is used afterwards we have to consider some additional deviation in the first
    # optimization
        
        
    #if opt_type == 'mintime' and mintime_opts["reopt_mintime_solution"]:
    #    w_veh_tmp = pars2["optim_opts"]["width_opt"] + (pars2["optim_opts"]["w_tr_reopt"] - pars2["optim_opts"]["w_veh_reopt"])
    #    w_veh_tmp += pars2["optim_opts"]["w_add_spl_regr"]
    #    pars_tmp = copy.deepcopy(pars2)
    #    pars_tmp["optim_opts"]["width_opt"] = w_veh_tmp
    #else:
    #    pars_tmp = pars2
    #
    ## call optimization
    #if opt_type == 'mincurv':
    #    alpha_opt2 = tph.opt_min_curv.opt_min_curv(reftrack=reftrack_interp2,
    #                                              normvectors=normvec_normalized_interp2,
    #                                              A=a_interp2,
    #                                              kappa_bound=pars2["veh_params"]["curvlim"],
    #                                              w_veh=pars2["optim_opts"]["width_opt"],
    #                                              print_debug=debug,
    #                                              plot_debug=plot_opts["mincurv_curv_lin"])[0]
    #
    #elif opt_type == 'mincurv_iqp':
    #    alpha_opt2, reftrack_interp2, normvec_normalized_interp2 = tph.iqp_handler.\
    #        iqp_handler(reftrack=reftrack_interp2,
    #                    normvectors=normvec_normalized_interp2,
    #                    A=a_interp2,
    #                    kappa_bound=pars2["veh_params"]["curvlim"],
    #                    w_veh=pars2["optim_opts"]["width_opt"],
    #                    print_debug=debug,
    #                    plot_debug=plot_opts["mincurv_curv_lin"],
    #                    stepsize_interp=pars2["stepsize_opts"]["stepsize_reg"],
    #                    iters_min=pars2["optim_opts"]["iqp_iters_min"],
    #                    curv_error_allowed=pars2["optim_opts"]["iqp_curverror_allowed"])
    #
    #elif opt_type == 'shortest_path':
    #    alpha_opt2 = tph.opt_shortest_path.opt_shortest_path(reftrack=reftrack_interp2,
    #                                                        normvectors=normvec_normalized_interp2,
    #                                                        w_veh=pars2["optim_opts"]["width_opt"],
    #                                                        print_debug=debug)
    #
    #elif opt_type == 'mintime':
    #    # reftrack_interp, a_interp and normvec_normalized_interp are returned for the case that non-regular sampling was
    #    # applied
    #    alpha_opt2, v_opt2, reftrack_interp2, a_interp_tmp2, normvec_normalized_interp2 = opt_mintime_traj.src.opt_mintime.\
    #        opt_mintime(reftrack=reftrack_interp2,
    #                    coeffs_x=coeffs_x_interp2,
    #                    coeffs_y=coeffs_y_interp2,
    #                    normvectors=normvec_normalized_interp2,
    #                    pars=pars_tmp,
    #                    tpamap_path=file_paths2["tpamap"],
    #                    tpadata_path=file_paths2["tpadata"],
    #                    export_path=file_paths2["mintime_export"],
    #                    print_debug=debug,
    #                    plot_debug=plot_opts["mintime_plots"])
    #
    #    # replace a_interp if necessary
    #    if a_interp_tmp2 is not None:
    #        a_interp2 = a_interp_tmp2
    #
    #else:
    #    raise ValueError('Unknown optimization type!')
    #    # alpha_opt = np.zeros(reftrack_interp.shape[0])
    
    #----------------------------------------------------------------------------------------------------------------------
    # INTERPOLATE SPLINES TO SMALL DISTANCES BETWEEN RACELINE POINTS 1 #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    raceline_interp1, a_opt1, coeffs_x_opt1, coeffs_y_opt1, spline_inds_opt_interp1, t_vals_opt_interp1, s_points_opt_interp1,\
        spline_lengths_opt1, el_lengths_opt_interp1 = tph.create_raceline.\
        create_raceline(refline=reftrack_interp1[:, :2],
                        normvectors=normvec_normalized_interp1,
                        alpha=alpha_opt1,
                        stepsize_interp=pars1["stepsize_opts"]["stepsize_interp_after_opt"])
    
    #----------------------------------------------------------------------------------------------------------------------
    # INTERPOLATE SPLINES TO SMALL DISTANCES BETWEEN RACELINE POINTS 2 #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    raceline_interp2, a_opt2, coeffs_x_opt2, coeffs_y_opt2, spline_inds_opt_interp2, t_vals_opt_interp2, s_points_opt_interp2,\
        spline_lengths_opt2, el_lengths_opt_interp2 = tph.create_raceline.\
        create_raceline(refline=reftrack_interp2[:, :2],
                        normvectors=normvec_normalized_interp2,
                        alpha=alpha_opt2,
                        stepsize_interp=pars2["stepsize_opts"]["stepsize_interp_after_opt"])
    
    
    # ----------------------------------------------------------------------------------------------------------------------
    # CALCULATE HEADING AND CURVATURE 1 --------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # calculate heading and curvature (analytically)
    psi_vel_opt1, kappa_opt1 = tph.calc_head_curv_an.\
        calc_head_curv_an(coeffs_x=coeffs_x_opt1,
                          coeffs_y=coeffs_y_opt1,
                          ind_spls=spline_inds_opt_interp1,
                          t_spls=t_vals_opt_interp1)
    
    #----------------------------------------------------------------------------------------------------------------------
    # CALCULATE HEADING AND CURVATURE 2 --------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # calculate heading and curvature (analytically)
    psi_vel_opt2, kappa_opt2 = tph.calc_head_curv_an.\
        calc_head_curv_an(coeffs_x=coeffs_x_opt2,
                          coeffs_y=coeffs_y_opt2,
                          ind_spls=spline_inds_opt_interp2,
                          t_spls=t_vals_opt_interp2)
    
    #----------------------------------------------------------------------------------------------------------------------
    # CALCULATE VELOCITY AND ACCELERATION PROFILE -------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    if opt_type == 'mintime' and not mintime_opts["recalc_vel_profile_by_tph"]:
        # interpolation 1
        s_splines1 = np.cumsum(spline_lengths_opt1)
        s_splines1 = np.insert(s_splines1, 0, 0.0)
        vx_profile_opt1 = np.interp(s_points_opt_interp1, s_splines1[:-1], v_opt1)
    
    else:
        vx_profile_opt1 = tph.calc_vel_profile.\
            calc_vel_profile(ggv=ggv,
                             ax_max_machines=ax_max_machines,
                             v_max=pars1["veh_params"]["v_max"],
                             kappa=kappa_opt1,
                             el_lengths=el_lengths_opt_interp1,
                             closed=True,
                             filt_window=pars1["vel_calc_opts"]["vel_profile_conv_filt_window"],
                             dyn_model_exp=pars1["vel_calc_opts"]["dyn_model_exp"],
                             drag_coeff=pars1["veh_params"]["dragcoeff"],
                             m_veh=pars1["veh_params"]["mass"])
            
    if opt_type == 'mintime' and not mintime_opts["recalc_vel_profile_by_tph"]:
        # interpolation 2
        s_splines2 = np.cumsum(spline_lengths_opt1)
        s_splines2 = np.insert(s_splines2, 0, 0.0)
        vx_profile_opt2 = np.interp(s_points_opt_interp2, s_splines2[:-1], v_opt2)
    
    else:
        vx_profile_opt2 = tph.calc_vel_profile.\
            calc_vel_profile(ggv=ggv,
                             ax_max_machines=ax_max_machines,
                             v_max=pars2["veh_params"]["v_max"],
                             kappa=kappa_opt2,
                             el_lengths=el_lengths_opt_interp2,
                             closed=True,
                             filt_window=pars2["vel_calc_opts"]["vel_profile_conv_filt_window"],
                             dyn_model_exp=pars2["vel_calc_opts"]["dyn_model_exp"],
                             drag_coeff=pars2["veh_params"]["dragcoeff"],
                             m_veh=pars2["veh_params"]["mass"])
    
    # calculate longitudinal acceleration profile 1
    vx_profile_opt_cl1 = np.append(vx_profile_opt1, vx_profile_opt1[0])
    ax_profile_opt1 = tph.calc_ax_profile.calc_ax_profile(vx_profile=vx_profile_opt_cl1,
                                                         el_lengths=el_lengths_opt_interp1,
                                                         eq_length_output=False)
    
    # calculate laptime 1
    t_profile_cl1 = tph.calc_t_profile.calc_t_profile(vx_profile=vx_profile_opt1,
                                                     ax_profile=ax_profile_opt1,
                                                     el_lengths=el_lengths_opt_interp1)
    print("INFO: Estimated laptime 1: %.2fs" % t_profile_cl1[-1])
    
    # calculate longitudinal acceleration profile 2
    vx_profile_opt_cl2 = np.append(vx_profile_opt2, vx_profile_opt2[0])
    ax_profile_opt2 = tph.calc_ax_profile.calc_ax_profile(vx_profile=vx_profile_opt_cl2,
                                                         el_lengths=el_lengths_opt_interp2,
                                                         eq_length_output=False)
    
    # calculate laptime 2
    t_profile_cl2 = tph.calc_t_profile.calc_t_profile(vx_profile=vx_profile_opt2,
                                                     ax_profile=ax_profile_opt2,
                                                     el_lengths=el_lengths_opt_interp2)
    print("INFO: Estimated laptime 2: %.2fs" % t_profile_cl2[-1])
    
    if plot_opts["racetraj_vel"]:
        s_points1 = np.cumsum(el_lengths_opt_interp1[:-1])
        s_points1 = np.insert(s_points1, 0, 0.0)
        
        s_points2 = np.cumsum(el_lengths_opt_interp2[:-1])
        s_points2 = np.insert(s_points2, 0, 0.0)
    
        plt.figure(10)
        plt.plot(s_points1, vx_profile_opt1)
        plt.plot(s_points2, vx_profile_opt2)
        plt.grid()
        plt.xlabel("distance in m")
        plt.legend(["Car 1 vx in m/s", "Car2 vx in m/s"])
        plt.show()
        
        plt.figure(11)
        plt.plot(s_points1, ax_profile_opt1)
        plt.plot(s_points2, ax_profile_opt2)
        plt.grid()
        plt.xlabel("distance in m")
        plt.legend(["Car 1 ax in m/s2", "Car2 ax in m/s2"])
        plt.show()
        
        
        plt.figure(12)
        plt.plot(s_points1, t_profile_cl1[:-1])
        plt.plot(s_points2, t_profile_cl2[:-1])
        plt.grid()
        plt.xlabel("distance in m")
        plt.legend(["Car 1 t in s", "Car 2 t in s"])
        plt.show()
    
    # ----------------------------------------------------------------------------------------------------------------------
    # CALCULATE LAP TIMES (AT DIFFERENT SCALES AND TOP SPEEDS) -RN/ UNCHANGED-------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    if lap_time_mat_opts["use_lap_time_mat"]:
        # simulate lap times
        ggv_scales = np.linspace(lap_time_mat_opts['gg_scale_range'][0],
                                 lap_time_mat_opts['gg_scale_range'][1],
                                 int((lap_time_mat_opts['gg_scale_range'][1] - lap_time_mat_opts['gg_scale_range'][0])
                                     / lap_time_mat_opts['gg_scale_stepsize']) + 1)
        top_speeds = np.linspace(lap_time_mat_opts['top_speed_range'][0] / 3.6,
                                 lap_time_mat_opts['top_speed_range'][1] / 3.6,
                                 int((lap_time_mat_opts['top_speed_range'][1] - lap_time_mat_opts['top_speed_range'][0])
                                     / lap_time_mat_opts['top_speed_stepsize']) + 1)
    
        # setup results matrix
        lap_time_matrix = np.zeros((top_speeds.shape[0] + 1, ggv_scales.shape[0] + 1))
    
        # write parameters in first column and row
        lap_time_matrix[1:, 0] = top_speeds * 3.6
        lap_time_matrix[0, 1:] = ggv_scales
    
        for i, top_speed in enumerate(top_speeds):
            for j, ggv_scale in enumerate(ggv_scales):
                tph.progressbar.progressbar(i*ggv_scales.shape[0] + j,
                                            top_speeds.shape[0] * ggv_scales.shape[0],
                                            prefix="Simulating laptimes ")
    
                ggv_mod = np.copy(ggv)
                ggv_mod[:, 1:] *= ggv_scale
    
                vx_profile_opt = tph.calc_vel_profile.\
                    calc_vel_profile(ggv=ggv_mod,
                                     ax_max_machines=ax_max_machines,
                                     v_max=top_speed,
                                     kappa=kappa_opt1,
                                     el_lengths=el_lengths_opt_interp1,
                                     dyn_model_exp=pars1["vel_calc_opts"]["dyn_model_exp"],
                                     filt_window=pars1["vel_calc_opts"]["vel_profile_conv_filt_window"],
                                     closed=True,
                                     drag_coeff=pars1["veh_params"]["dragcoeff"],
                                     m_veh=pars1["veh_params"]["mass"])
    
                # calculate longitudinal acceleration profile
                vx_profile_opt_cl = np.append(vx_profile_opt, vx_profile_opt[0])
                ax_profile_opt = tph.calc_ax_profile.calc_ax_profile(vx_profile=vx_profile_opt_cl,
                                                                     el_lengths=el_lengths_opt_interp1,
                                                                     eq_length_output=False)
    
                # calculate lap time
                t_profile_cl = tph.calc_t_profile.calc_t_profile(vx_profile=vx_profile_opt,
                                                                 ax_profile=ax_profile_opt,
                                                                 el_lengths=el_lengths_opt_interp1)
    
                # store entry in lap time matrix
                lap_time_matrix[i + 1, j + 1] = t_profile_cl[-1]
    
        # store lap time matrix to file
        np.savetxt(file_paths1["lap_time_mat_export"], lap_time_matrix, delimiter=",", fmt="%.3f")
    
    # ----------------------------------------------------------------------------------------------------------------------
    # DATA POSTPROCESSING 1 --------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # arrange data into one trajectory
    trajectory_opt1 = np.column_stack((s_points_opt_interp1,
                                      raceline_interp1,
                                      psi_vel_opt1,
                                      kappa_opt1,
                                      vx_profile_opt1,
                                      ax_profile_opt1))
    spline_data_opt1 = np.column_stack((spline_lengths_opt1, coeffs_x_opt1, coeffs_y_opt1))
    
    # create a closed race trajectory array
    traj_race_cl1 = np.vstack((trajectory_opt1, trajectory_opt1[0, :]))
    traj_race_cl1[-1, 0] = np.sum(spline_data_opt1[:, 0])  # set correct length
    
    # ----------------------------------------------------------------------------------------------------------------------
    # DATA POSTPROCESSING 2 -----------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # arrange data into one trajectory
    trajectory_opt2 = np.column_stack((s_points_opt_interp2,
                                      raceline_interp2,
                                      psi_vel_opt2,
                                      kappa_opt2,
                                      vx_profile_opt2,
                                      ax_profile_opt2))
    spline_data_opt2 = np.column_stack((spline_lengths_opt2, coeffs_x_opt2, coeffs_y_opt2))
    
    # create a closed race trajectory array
    traj_race_cl2 = np.vstack((trajectory_opt2, trajectory_opt2[0, :]))
    traj_race_cl2[-1, 0] = np.sum(spline_data_opt2[:, 0])  # set correct length
    
    # print end time
    print("INFO: Runtime from import to final trajectory was %.2fs" % (time.perf_counter() - t_start))
    
    # ----------------------------------------------------------------------------------------------------------------------
    # CHECK TRAJECTORY 1 --------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    bound11, bound12 = helper_funcs_glob.src.check_traj.\
        check_traj(reftrack=reftrack_interp1,
                   reftrack_normvec_normalized=normvec_normalized_interp1,
                   length_veh=pars1["veh_params"]["length"],
                   width_veh=pars1["veh_params"]["width"],
                   debug=debug,
                   trajectory=trajectory_opt1,
                   ggv=ggv,
                   ax_max_machines=ax_max_machines,
                   v_max=pars1["veh_params"]["v_max"],
                   curvlim=pars1["veh_params"]["curvlim"],
                   mass_veh=pars1["veh_params"]["mass"],
                   dragcoeff=pars1["veh_params"]["dragcoeff"])
    
    # ----------------------------------------------------------------------------------------------------------------------
    # CHECK TRAJECTORY 2---------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    bound21, bound22 = helper_funcs_glob.src.check_traj.\
        check_traj(reftrack=reftrack_interp2,
                   reftrack_normvec_normalized=normvec_normalized_interp2,
                   length_veh=pars2["veh_params"]["length"],
                   width_veh=pars2["veh_params"]["width"],
                   debug=debug,
                   trajectory=trajectory_opt2,
                   ggv=ggv,
                   ax_max_machines=ax_max_machines,
                   v_max=pars2["veh_params"]["v_max"],
                   curvlim=pars2["veh_params"]["curvlim"],
                   mass_veh=pars2["veh_params"]["mass"],
                   dragcoeff=pars2["veh_params"]["dragcoeff"])
    
    # ----------------------------------------------------------------------------------------------------------------------
    # EXPORT 1------------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # export race trajectory  to CSV
    if "traj_race_export" in file_paths1.keys():
        helper_funcs_glob.src.export_traj_race.export_traj_race(file_paths=file_paths1,
                                                                traj_race=traj_race_cl1)
    
    # if requested, export trajectory including map information (via normal vectors) to CSV
    if "traj_ltpl_export" in file_paths1.keys():
        helper_funcs_glob.src.export_traj_ltpl.export_traj_ltpl(file_paths=file_paths1,
                                                                spline_lengths_opt=spline_lengths_opt1,
                                                                trajectory_opt=trajectory_opt1,
                                                                reftrack=reftrack_interp1,
                                                                normvec_normalized=normvec_normalized_interp1,
                                                                alpha_opt=alpha_opt1)
    
    # ----------------------------------------------------------------------------------------------------------------------
    # EXPORT 2-------------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # export race trajectory  to CSV
    if "traj_race_export" in file_paths2.keys():
        helper_funcs_glob.src.export_traj_race.export_traj_race(file_paths=file_paths2,
                                                                traj_race=traj_race_cl2)
    
    # if requested, export trajectory including map information (via normal vectors) to CSV
    if "traj_ltpl_export" in file_paths2.keys():
        helper_funcs_glob.src.export_traj_ltpl.export_traj_ltpl(file_paths=file_paths2,
                                                                spline_lengths_opt=spline_lengths_opt2,
                                                                trajectory_opt=trajectory_opt2,
                                                                reftrack=reftrack_interp2,
                                                                normvec_normalized=normvec_normalized_interp2,
                                                                alpha_opt=alpha_opt2)
    
    print("INFO: Finished export of trajectory:", time.strftime("%H:%M:%S"))
    
    
    # ----------------------------------------------------------------------------------------------------------------------
    # PLOT RESULTS --- HERE ----------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    
    # get bound of imported map (for reference in final plot) 1
    bound11_imp = None
    bound12_imp = None
    
    if plot_opts["imported_bounds"]:
        # try to extract four times as many points as in the interpolated version (in order to hold more details)
        n_skip1 = max(int(reftrack_imp1.shape[0] / (bound11.shape[0] * 4)), 1)
    
        _, _, _, normvec_imp1 = tph.calc_splines.calc_splines(path=np.vstack((reftrack_imp1[::n_skip1, 0:2],
                                                                             reftrack_imp1[0, 0:2])))
    
        bound11_imp = reftrack_imp1[::n_skip1, :2] + normvec_imp1 * np.expand_dims(reftrack_imp1[::n_skip1, 2], 1)
        bound12_imp = reftrack_imp1[::n_skip1, :2] - normvec_imp1 * np.expand_dims(reftrack_imp1[::n_skip1, 3], 1)
    
    # get bound of imported map (for reference in final plot) 2
    bound21_imp = None
    bound22_imp = None
    
    if plot_opts["imported_bounds"]:
        # try to extract four times as many points as in the interpolated version (in order to hold more details)
        n_skip2 = max(int(reftrack_imp2.shape[0] / (bound21.shape[0] * 4)), 1)
    
        _, _, _, normvec_imp2 = tph.calc_splines.calc_splines(path=np.vstack((reftrack_imp2[::n_skip2, 0:2],
                                                                             reftrack_imp2[0, 0:2])))
    
        bound21_imp = reftrack_imp2[::n_skip2, :2] + normvec_imp2 * np.expand_dims(reftrack_imp2[::n_skip2, 2], 1)
        bound22_imp = reftrack_imp2[::n_skip2, :2] - normvec_imp2 * np.expand_dims(reftrack_imp2[::n_skip2, 3], 1)
    
    ###
    c1_fast = False
    c2_fast = False
    if (t_profile_cl1[-1]<t_profile_cl2[-1]):
        print("Car 1 has a better time, and it takes the following trajectory to do so")
        c1_fast = True
    elif (t_profile_cl2[-1]<t_profile_cl1[-1]):
        print("Car 2 has a better time, and it takes the following trajectory to do so")
        c2_fast = True
    else:
        print("They both take the same time, so here are both the outs")
        c1_fast = True
        c2_fast = True
    ####
    # plot results
    if (c1_fast):
        helper_funcs_glob.src.result_plots.result_plots(plot_opts=plot_opts,
                                                    width_veh_opt=pars1["optim_opts"]["width_opt"],
                                                    width_veh_real=pars1["veh_params"]["width"],
                                                    refline=reftrack_interp1[:, :2],
                                                    bound1_imp=bound11_imp,
                                                    bound2_imp=bound12_imp,
                                                    bound1_interp=bound11,
                                                    bound2_interp=bound12,
                                                    trajectory=trajectory_opt1)
    if (c2_fast):
        helper_funcs_glob.src.result_plots.result_plots(plot_opts=plot_opts,
                                                    width_veh_opt=pars2["optim_opts"]["width_opt"],
                                                    width_veh_real=pars2["veh_params"]["width"],
                                                    refline=reftrack_interp2[:, :2],
                                                    bound1_imp=bound21_imp,
                                                    bound2_imp=bound22_imp,
                                                    bound1_interp=bound21,
                                                    bound2_interp=bound22,
                                                    trajectory=trajectory_opt2)
    print("INFO: Your comparative simulation is done")
        

if __name__=='__main__':
	pass