import numpy as np
import matplotlib.pyplot as plt


def result_plots_mintime_compare(pars1: dict,
                         reftrack1: np.ndarray,
                         s1: np.ndarray,
                         t1: np.ndarray,
                         x1: np.ndarray,
                         u1: np.ndarray,
                         ax1: np.ndarray,
                         ay1: np.ndarray,
                         atot1: np.ndarray,
                         tf1: np.ndarray,
                         ec1: np.ndarray,
                         pars2: dict,
                         reftrack2: np.ndarray,
                         s2: np.ndarray,
                         t2: np.ndarray,
                         x2: np.ndarray,
                         u2: np.ndarray,
                         ax2: np.ndarray,
                         ay2: np.ndarray,
                         atot2: np.ndarray,
                         tf2: np.ndarray,
                         ec2: np.ndarray,
                         pwr2: dict = None,
                         pwr1: dict = None) -> None:

    """
    Created by:
    Fabian Christ

    Extended by:
    Thomas Herrmann (thomas.herrmann@tum.de)

    Documentation:
    This function plots several figures containing relevant trajectory information after trajectory optimization.

    Inputs:
    pars:       parameters dictionary
    reftrack:   contains the information of the reftrack -> [x, y, w_tr_right, w_tr_left]
    s:          contains the curvi-linear distance along the trajectory
    t:          contains the time along the trajectory
    x:          contains all state variables along the trajectory
    u:          contains all control variables along the trajectory
    ax:         contains the longitudinal acceleration along the trajectory
    ay:         contains the lateral acceleration along the trajectory
    atot:       contains the total acceleration along the trajectory
    tf:         contains all tire forces along the trajectory
    ec:         contains the used energy along the trajectory
    """

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT OPTIONS -----------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    plt.rcParams['axes.labelsize'] = 10.0
    plt.rcParams['axes.titlesize'] = 11.0
    plt.rcParams['legend.fontsize'] = 10.0
    plt.rcParams['figure.figsize'] = 25 / 2.54, 20 / 2.54

    plot_opts = {"v_a_t": True,
                 "general": True,
                 "lateral_distance": True,
                 "power": True,
                 "kamm_circle": True,
                 "tire_forces": True,
                 "tire_forces_longitudinal": True,
                 "tire_forces_dynamic": True,
                 "energy_consumption": True,
                 "pwr_states": True,
                 "pwr_soc": True,
                 "pwr_losses": True}

    # ------------------------------------------------------------------------------------------------------------------
    # EXTRACT PLOT DATA 1----------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------
    
    # state variables
    v1 = x1[:, 0]
    beta1 = x1[:, 1]
    omega_z1 = x1[:, 2]
    n1 = x1[:, 3]
    xi1 = x1[:, 4]
    if pars1["pwr_params_mintime"]["pwr_behavior"]:
        temp_mot1 = x1[:, 5]
        temp_batt1 = x1[:, 6]
        temp_inv1 = x1[:, 7]
        temp_radiators_cool_mi1 = x1[:, 8]
        temp_radiators_cool_b1 = x1[:, 9]
        soc_batt1 = x1[:, 10]

    # control variables
    delta1 = np.append(u1[:, 0], u1[0, 0])
    f_drive1 = np.append(u1[:, 1], u1[0, 1])
    f_brake1 = np.append(u1[:, 2], u1[0, 2])
    gamma_y1 = np.append(u1[:, 3], u1[0, 3])
    
    # tire forces
    tf_x_fl1 = tf1[:, 0]
    tf_y_fl1 = tf1[:, 1]
    tf_z_fl1 = tf1[:, 2]
    tf_x_fr1 = tf1[:, 3]
    tf_y_fr1 = tf1[:, 4]
    tf_z_fr1 = tf1[:, 5] 
    tf_x_rl1 = tf1[:, 6]
    tf_y_rl1 = tf1[:, 7]
    tf_z_rl1 = tf1[:, 8]
    tf_x_rr1 = tf1[:, 9]
    tf_y_rr1 = tf1[:, 10]
    tf_z_rr1 = tf1[:, 11]

    # parameters
    g1 = pars1["veh_params"]["g"]
    veh1 = pars1["vehicle_params_mintime"]
    tire1 = pars1["tire_params_mintime"]
    
    # ------------------------------------------------------------------------------------------------------------------
    # EXTRACT PLOT DATA 2----------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------
    
    # state variables
    v2 = x2[:, 0]
    beta2 = x2[:, 1]
    omega_z2 = x2[:, 2]
    n2 = x2[:, 3]
    xi2 = x2[:, 4]
    if pars2["pwr_params_mintime"]["pwr_behavior"]:
        temp_mot2 = x2[:, 5]
        temp_batt2 = x2[:, 6]
        temp_inv2 = x2[:, 7]
        temp_radiators_cool_mi2 = x2[:, 8]
        temp_radiators_cool_b2 = x2[:, 9]
        soc_batt2 = x2[:, 10]

    # control variables
    delta2 = np.append(u2[:, 0], u2[0, 0])
    f_drive2 = np.append(u2[:, 1], u2[0, 1])
    f_brake2 = np.append(u2[:, 2], u2[0, 2])
    gamma_y2 = np.append(u2[:, 3], u2[0, 3])
    
    # tire forces
    tf_x_fl2 = tf2[:, 0]
    tf_y_fl2 = tf2[:, 1]
    tf_z_fl2 = tf2[:, 2]
    tf_x_fr2 = tf2[:, 3]
    tf_y_fr2 = tf2[:, 4]
    tf_z_fr2 = tf2[:, 5] 
    tf_x_rl2 = tf2[:, 6]
    tf_y_rl2 = tf2[:, 7]
    tf_z_rl2 = tf2[:, 8]
    tf_x_rr2 = tf2[:, 9]
    tf_y_rr2 = tf2[:, 10]
    tf_z_rr2 = tf2[:, 11]

    # parameters
    g2 = pars2["veh_params"]["g"]
    veh2 = pars2["vehicle_params_mintime"]
    tire2 = pars2["tire_params_mintime"]


    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: VELOCITY + LONGITUDINAL ACCELERATION + LATERAL ACCELERATION + TOTAL ACCELERATION + TIME --------------------
    # ------------------------------------------------------------------------------------------------------------------

    if plot_opts["v_a_t"]:

        plt.figure(1)
        plt.clf()
        data1, = plt.plot(s1, v1, label=r'$\it{v}$' + ' in ' + r'$\it{\frac{m}{s}}$')
        data2, = plt.plot(s1, ax1, label=r'$\it{a_x}$' + ' in ' + r'$\it{\frac{m}{s^2}}$')
        data3, = plt.plot(s1, ay1, label=r'$\it{a_y}$' + ' in ' + r'$\it{\frac{m}{s^2}}$')
        data4, = plt.plot(s1, atot1, label=r'$\it{a_{tot}}$' + ' in ' + r'$\it{\frac{m}{s^2}}$')
        data5, = plt.plot(s1, t1, label=r'$\it{t}$' + ' in ' + r'$\it{s}$')
        
        data11, = plt.plot(s2, v2, label=r'$\it{v}$' + ' in ' + r'$\it{\frac{m}{s}}$')
        data22, = plt.plot(s2, ax2, label=r'$\it{a_x}$' + ' in ' + r'$\it{\frac{m}{s^2}}$')
        data33, = plt.plot(s2, ay2, label=r'$\it{a_y}$' + ' in ' + r'$\it{\frac{m}{s^2}}$')
        data44, = plt.plot(s2, atot2, label=r'$\it{a_{tot}}$' + ' in ' + r'$\it{\frac{m}{s^2}}$')
        data55, = plt.plot(s2, t2, label=r'$\it{t}$' + ' in ' + r'$\it{s}$')
        
        plt.grid()
        plt.ylim(bottom=-15)
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        first_legend = plt.legend(handles=[data1, data2, data3, data4, data5], title="Car1", loc=2)
        
        plt.gca().add_artist(first_legend)
        
        plt.legend(handles=[data11, data22, data33, data44, data55], title = "Car 2")
        
        #plt.legend([r'$\it{v}$' + ' in ' + r'$\it{\frac{m}{s}}$',
        #            r'$\it{a_x}$' + ' in ' + r'$\it{\frac{m}{s^2}}$',
        #            r'$\it{a_y}$' + ' in ' + r'$\it{\frac{m}{s^2}}$',
        #            r'$\it{a_{tot}}$' + ' in ' + r'$\it{\frac{m}{s^2}}$',
        #            r'$\it{t}$' + ' in ' + r'$\it{s}$'])
        
        
        
        plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: SIDE SLIP ANGLE + YAW RATE + RELATIVE ANGLE TO TANGENT ON REFLINE + STEERING ANGLE -------------------------
    # ------------------------------------------------------------------------------------------------------------------

    if plot_opts["general"]:

        plt.figure(2)
        plt.clf()
        plt.subplot(221)
        plt.plot(s1, beta1 * 180 / np.pi)
        plt.plot(s2, beta2 * 180 / np.pi)
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel('side slip angle ' + r'$\beta$' + ' in ' + r'$\it{°}$')
        plt.grid()
        plt.subplot(222)
        plt.plot(s1, omega_z1 * 180 / np.pi)
        plt.plot(s2, omega_z2 * 180 / np.pi)
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel('yaw rate ' + r'$\omega_{z}$' + ' in ' + r'$\it{\frac{°}{s}}$')
        plt.grid()
        plt.subplot(223)
        plt.plot(s1, xi1 * 180 / np.pi)
        plt.plot(s2, xi2 * 180 / np.pi)
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel('relative angle to tangent on reference line ' + r'$\xi$' + ' in ' + r'$\it{°}$')
        plt.grid()
        plt.subplot(224)
        plt.step(s1, delta1 * 180 / np.pi, where='post')
        plt.step(s2, delta2 * 180 / np.pi, where='post')
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel('steering angle ' + r'$\delta$' + ' in ' + r'$\it{°}$')
        plt.grid()
        plt.legend(['Car 1', 'Car 2'])
        plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: LATERAL DISTANCE TO REFERENCE LINE + ROAD BOUNDARIES 1----------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    if plot_opts["lateral_distance"]:

        plt.figure(3)
        plt.clf()
        plt.plot(s1, n1, color = 'green')
        plt.plot(s2, n2, color = 'red')
        reftrack_cl = np.vstack((reftrack1, reftrack1[0, :]))
        plt.plot(s1, reftrack_cl[:, 3], color='black')
        plt.plot(s1, reftrack_cl[:, 3] - pars1["optim_opts"]["width_opt"] / 2, color='grey')
        plt.plot(s1, -reftrack_cl[:, 2], color='black')
        plt.plot(s1, -reftrack_cl[:, 2] + pars1["optim_opts"]["width_opt"] / 2, color='grey')
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel('lateral distance to reference line ' + r'$\it{n}$' + ' in ' + r'$\it{m}$')
        plt.legend(['raceline - Car 1', 'raceline - Car 2', 'road boundaries', 'road boundaries - safety margin'], ncol=1, loc=4)
        plt.grid()
        plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: KAMM's CIRCLE ----------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    if plot_opts["kamm_circle"]:

        plt.figure(5)
        plt.clf()
        plt.suptitle("Kamm's Circle")
        plt.subplot(221)
        circle1 = plt.Circle((0, 0), 1, fill=False)
        fig = plt.gcf()
        ax = fig.gca()
        ax.add_artist(circle1)
        plt.plot(tf_y_fl1 / (tf_z_fl1 * pars1["optim_opts"]["mue"]),
                 tf_x_fl1 / (tf_z_fl1 * pars1["optim_opts"]["mue"]), '^:', color = 'green')
        plt.plot(tf_y_fl2 / (tf_z_fl2 * pars2["optim_opts"]["mue"]),
                 tf_x_fl2 / (tf_z_fl2 * pars2["optim_opts"]["mue"]), '^:', color = 'blue')
        plt.xlim(-1.2, 1.2)
        plt.ylim(-1.2, 1.2)
        plt.xlabel(r'$\it{\frac{F_{y}}{F_{ymax}}}$')
        plt.ylabel(r'$\it{\frac{F_{x}}{F_{xmax}}}$')
        plt.axis('equal')
        plt.grid()

        plt.subplot(222)
        circle1 = plt.Circle((0, 0), 1, fill=False)
        fig = plt.gcf()
        ax = fig.gca()
        ax.add_artist(circle1)
        plt.plot(tf_y_fr1 / (tf_z_fr1 * pars1["optim_opts"]["mue"]),
                 tf_x_fr1 / (tf_z_fr1 * pars1["optim_opts"]["mue"]), '^:', color = 'green')
        plt.plot(tf_y_fr2 / (tf_z_fr2 * pars2["optim_opts"]["mue"]),
                 tf_x_fr2 / (tf_z_fr2 * pars2["optim_opts"]["mue"]), '^:', color = 'blue')
        plt.xlim(-1.2, 1.2)
        plt.ylim(-1.2, 1.2)
        plt.xlabel(r'$\it{\frac{F_{y}}{F_{ymax}}}$')
        plt.ylabel(r'$\it{\frac{F_{x}}{F_{xmax}}}$')
        plt.axis('equal')
        plt.grid()

        plt.subplot(223)
        circle1 = plt.Circle((0, 0), 1, fill=False)
        fig = plt.gcf()
        ax = fig.gca()
        ax.add_artist(circle1)
        plt.plot(tf_y_rl1 / (tf_z_rl1 * pars1["optim_opts"]["mue"]),
                 tf_x_rl1 / (tf_z_rl1 * pars1["optim_opts"]["mue"]), '^:',color='green')
        plt.plot(tf_y_rl2 / (tf_z_rl2 * pars2["optim_opts"]["mue"]),
                 tf_x_rl2 / (tf_z_rl2 * pars2["optim_opts"]["mue"]), '^:', color='blue')
        plt.xlim(-1.2, 1.2)
        plt.ylim(-1.2, 1.2)
        plt.xlabel(r'$\it{\frac{F_{y}}{F_{ymax}}}$')
        plt.ylabel(r'$\it{\frac{F_{x}}{F_{xmax}}}$')
        plt.axis('equal')
        plt.grid()

        plt.subplot(224)
        circle1 = plt.Circle((0, 0), 1, fill=False)
        fig = plt.gcf()
        ax = fig.gca()
        ax.add_artist(circle1)
        plt.plot(tf_y_rr1 / (tf_z_rr1 * pars1["optim_opts"]["mue"]),
                 tf_x_rr1 / (tf_z_rr1 * pars1["optim_opts"]["mue"]), '^:')
        plt.plot(tf_y_rr2 / (tf_z_rr2 * pars2["optim_opts"]["mue"]),
                 tf_x_rr2 / (tf_z_rr2 * pars2["optim_opts"]["mue"]), '^:', color='blue')
        plt.xlim(-1.2, 1.2)
        plt.ylim(-1.2, 1.2)
        plt.xlabel(r'$\it{\frac{F_{y}}{F_{ymax}}}$')
        plt.ylabel(r'$\it{\frac{F_{x}}{F_{xmax}}}$')
        plt.axis('equal')
        plt.grid()
        plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: TIRE FORCES (LONGITUDINAL + LATERAL + NORMAL) --------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    if plot_opts["tire_forces"]:

        plt.figure(6)
        plt.clf()
        plt.suptitle("Tire Forces")
        plt.subplot(221)
        plt.plot(s1, tf_x_fl1)
        plt.plot(s1, tf_y_fl1)
        plt.plot(s1, tf_z_fl1)
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F_{i}}$' + ' in ' + r'$\it{N}$')
        plt.legend([r'$\it{F_{x}}$', r'$\it{F_{y}}$', r'$\it{F_{z}}$'], ncol=3, loc=4)
        plt.grid()

        plt.subplot(222)
        plt.plot(s1, tf_x_fr1)
        plt.plot(s1, tf_y_fr1)
        plt.plot(s1, tf_z_fr1)
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F_{i}}$' + ' in ' + r'$\it{N}$')
        plt.legend([r'$\it{F_{x}}$', r'$\it{F_{y}}$', r'$\it{F_{z}}$'], ncol=3, loc=4)
        plt.grid()

        plt.subplot(223)
        plt.plot(s1, tf_x_rl1)
        plt.plot(s1, tf_y_rl1)
        plt.plot(s1, tf_z_rl1)
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F_{i}}$' + ' in ' + r'$\it{N}$')
        plt.legend([r'$\it{F_{x}}$', r'$\it{F_{y}}$', r'$\it{F_{z}}$'], ncol=3, loc=4)
        plt.grid()

        plt.subplot(224)
        plt.plot(s1, tf_x_rr1)
        plt.plot(s1, tf_y_rr1)
        plt.plot(s1, tf_z_rr1)
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F_{i}}$' + ' in ' + r'$\it{N}$')
        plt.legend([r'$\it{F_{x}}$', r'$\it{F_{y}}$', r'$\it{F_{z}}$'], ncol=3, loc=4)
        plt.grid()
        plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: TIRE FORCES (LONGITUDINAL) ---------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    if plot_opts["tire_forces_longitudinal"]:
        plt.figure(7)
        plt.step(s1, f_drive1 / 1000, where="post")
        plt.step(s1, f_brake1 / 1000, where='post')
        plt.step(s1, (f_drive1 + f_brake1) / 1000, where='post')
        plt.plot(s1, veh1["power_max"] / (v1 * 1000), linewidth=0.5)
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F}$' + ' in ' + r'$\it{kN}$')
        plt.legend([r'$\it{F_{drive}}$', r'$\it{F_{brake}}$',
                    r'$\it{F_{drive}}$' + " + " + r'$\it{F_{brake}}$',
                    r'$\it{F_{P_{max}}}$'], ncol=1, loc=4)
        plt.grid()
        plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: DYNAMIC WHEEL LOAD TRANSFER 1-----------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    if plot_opts["tire_forces_dynamic"]:

        f_xroll1 = tire1["c_roll"] * pars1["veh_params"]["mass"] * g1
        f_xdrag1 = pars1["veh_params"]["dragcoeff"] * v1 ** 2

        f_zlift_fl1 = 0.5 * veh1["liftcoeff_front"] * v1 ** 2
        f_zlift_fr1 = 0.5 * veh1["liftcoeff_front"] * v1 ** 2
        f_zlift_rl1 = 0.5 * veh1["liftcoeff_rear"] * v1 ** 2
        f_zlift_rr1 = 0.5 * veh1["liftcoeff_rear"] * v1 ** 2

        f_zlong_fl1 = -0.5 * veh1["cog_z"] / veh1["wheelbase"] * (f_drive1 + f_brake1 - f_xroll1 - f_xdrag1)
        f_zlong_fr1 = -0.5 * veh1["cog_z"] / veh1["wheelbase"] * (f_drive1 + f_brake1 - f_xroll1 - f_xdrag1)
        f_zlong_rl1 = 0.5 * veh1["cog_z"] / veh1["wheelbase"] * (f_drive1 + f_drive1 - f_xroll1 - f_xdrag1)
        f_zlong_rr1 = 0.5 * veh1["cog_z"] / veh1["wheelbase"] * (f_drive1 + f_drive1 - f_xroll1 - f_xdrag1)

        f_zlat_fl1 = - veh1["k_roll"] * gamma_y1
        f_zlat_fr1 = veh1["k_roll"] * gamma_y1
        f_zlat_rl1 = - (1 - veh1["k_roll"]) * gamma_y1
        f_zlat_rr1 = (1 - veh1["k_roll"]) * gamma_y1

        plt.figure(8)
        plt.suptitle("Dynamic Wheel Load")
        plt.subplot(221)
        plt.plot(s1, f_zlift_fl1)
        plt.plot(s1, f_zlong_fl1)
        plt.plot(s1, f_zlat_fl1)
        plt.plot(s1, f_zlift_fl1 + f_zlong_fl1 + f_zlat_fl1, color='black')
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F_{i}}$' + ' in ' + r'$\it{N}$')

        plt.grid()
        plt.subplot(222)
        plt.plot(s1, f_zlift_fr1)
        plt.plot(s1, f_zlong_fr1)
        plt.plot(s1, f_zlat_fr1)
        plt.plot(s1, f_zlift_fr1 + f_zlong_fr1 + f_zlat_fr1, color='black')
        plt.xlabel(r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F_{i}}$' + ' in ' + r'$\it{N}$')

        plt.grid()
        plt.subplot(223)
        plt.plot(s1, f_zlift_rl1)
        plt.plot(s1, f_zlong_rl1)
        plt.plot(s1, f_zlat_rl1)
        plt.plot(s1, f_zlift_rl1 + f_zlong_rl1 + f_zlat_rl1, color='black')
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F_{i}}$' + ' in ' + r'$\it{N}$')

        plt.grid()
        plt.subplot(224)
        plt.plot(s1, f_zlift_rr1)
        plt.plot(s1, f_zlong_rr1)
        plt.plot(s1, f_zlat_rr1)
        plt.plot(s1, f_zlift_rr1 + f_zlong_rr1 + f_zlat_rr1, color='black')
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F_{i}}$' + ' in ' + r'$\it{N}$')
        plt.legend([r'$\it{F_{lift}}$', r'$\it{F_{dyn,long}}$', r'$\it{F_{dyn,lat}}$',
                    r'$\it{F_{lift}}$' + ' + ' + r'$\it{F_{dyn,long}}$' + ' + ' + r'$\it{F_{dyn,lat}}$'], ncol=2, loc=4)
        plt.grid()
        plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: DYNAMIC WHEEL LOAD TRANSFER 2-----------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    if plot_opts["tire_forces_dynamic"]:

        f_xroll2 = tire2["c_roll"] * pars1["veh_params"]["mass"] * g2
        f_xdrag2 = pars2["veh_params"]["dragcoeff"] * v2 ** 2

        f_zlift_fl2 = 0.5 * veh2["liftcoeff_front"] * v2 ** 2
        f_zlift_fr2 = 0.5 * veh2["liftcoeff_front"] * v2 ** 2
        f_zlift_rl2 = 0.5 * veh2["liftcoeff_rear"] * v2 ** 2
        f_zlift_rr2 = 0.5 * veh2["liftcoeff_rear"] * v2 ** 2

        f_zlong_fl2 = -0.5 * veh2["cog_z"] / veh2["wheelbase"] * (f_drive2 + f_brake2 - f_xroll2 - f_xdrag2)
        f_zlong_fr2 = -0.5 * veh2["cog_z"] / veh2["wheelbase"] * (f_drive2 + f_brake2 - f_xroll2 - f_xdrag2)
        f_zlong_rl2 = 0.5 * veh2["cog_z"] / veh2["wheelbase"] * (f_drive2 + f_drive2 - f_xroll2 - f_xdrag2)
        f_zlong_rr2 = 0.5 * veh2["cog_z"] / veh2["wheelbase"] * (f_drive2 + f_drive2 - f_xroll2 - f_xdrag2)

        f_zlat_fl2 = - veh2["k_roll"] * gamma_y2
        f_zlat_fr2 = veh2["k_roll"] * gamma_y2
        f_zlat_rl2 = - (1 - veh2["k_roll"]) * gamma_y2
        f_zlat_rr2 = (1 - veh2["k_roll"]) * gamma_y2

        plt.figure(8)
        plt.suptitle("Dynamic Wheel Load")
        plt.subplot(221)
        plt.plot(s2, f_zlift_fl2)
        plt.plot(s2, f_zlong_fl2)
        plt.plot(s2, f_zlat_fl2)
        plt.plot(s2, f_zlift_fl2 + f_zlong_fl2 + f_zlat_fl2, color='black')
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F_{i}}$' + ' in ' + r'$\it{N}$')

        plt.grid()
        plt.subplot(222)
        plt.plot(s2, f_zlift_fr2)
        plt.plot(s2, f_zlong_fr2)
        plt.plot(s2, f_zlat_fr2)
        plt.plot(s2, f_zlift_fr2 + f_zlong_fr2 + f_zlat_fr2, color='black')
        plt.xlabel(r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F_{i}}$' + ' in ' + r'$\it{N}$')

        plt.grid()
        plt.subplot(223)
        plt.plot(s2, f_zlift_rl2)
        plt.plot(s2, f_zlong_rl2)
        plt.plot(s2, f_zlat_rl2)
        plt.plot(s2, f_zlift_rl2 + f_zlong_rl2 + f_zlat_rl2, color='black')
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F_{i}}$' + ' in ' + r'$\it{N}$')

        plt.grid()
        plt.subplot(224)
        plt.plot(s2, f_zlift_rr2)
        plt.plot(s2, f_zlong_rr2)
        plt.plot(s2, f_zlat_rr2)
        plt.plot(s2, f_zlift_rr2 + f_zlong_rr2 + f_zlat_rr2, color='black')
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel(r'$\it{F_{i}}$' + ' in ' + r'$\it{N}$')
        plt.legend([r'$\it{F_{lift}}$', r'$\it{F_{dyn,long}}$', r'$\it{F_{dyn,lat}}$',
                    r'$\it{F_{lift}}$' + ' + ' + r'$\it{F_{dyn,long}}$' + ' + ' + r'$\it{F_{dyn,lat}}$'], ncol=2, loc=4)
        plt.grid()
        plt.show()
    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: ENERGY CONSUMPTION -----------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    if plot_opts["energy_consumption"]:

        plt.figure(9)
        plt.clf()
        plt.plot(s1, ec1)
        plt.plot(s2, ec2)
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel('energy consumption ' + r'$\it{ec}$' + ' in ' + r'$\it{Wh}$')
        plt.legend(['Car 1', 'Car2'])
        plt.grid()
        plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: POWER ------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    if plot_opts["power"]:
        plt.figure(4)
        plt.clf()
        plt.plot(s1, v1 * (f_drive1 + f_brake1) / 1000.0)
        plt.plot(s2, v2 * (f_drive2 + f_brake2) / 1000.0)
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$m$')
        plt.ylabel('power ' + r'$\it{P}$' + ' in ' + r'$kW$')
        plt.grid()
        plt.legend(['Car1', 'Car2'])
        if pwr1 is not None:
            plt.plot(s1[:-1], pwr1["batt"].p_loss_total + pwr1["batt"].p_out_batt)
            plt.legend([r'$\it{P_{wheel}}$', r'$\it{P_{system}}$'])
        plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: POWERTRAIN TEMPERATURES ------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    if pars1["pwr_params_mintime"]["pwr_behavior"] and plot_opts["pwr_states"]:

        plt.figure(10)
        plt.plot(s1, temp_mot1)
        plt.plot(s1, temp_batt1)
        plt.plot(s1, temp_inv1)
        plt.plot(s1, temp_radiators_cool_mi1)
        plt.plot(s1, temp_radiators_cool_b1)

        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel('component temperatures ' + r'$\it{T}$' + ' in ' + r'°C')
        plt.legend([r'$\it{T_\mathrm{Machine}}$', r'$\it{T_\mathrm{Battery}}$', r'$\it{T_\mathrm{Inverter}}$',
                    r'$\it{T_\mathrm{Fluid_{MI}}}$', r'$\it{T_\mathrm{Fluid_B}}$'])
        plt.grid()
        plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: SOC BATTERY ------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    if pars1["pwr_params_mintime"]["pwr_behavior"] and plot_opts["pwr_soc"]:

        plt.figure(11)
        plt.plot(s1, soc_batt1)
        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.ylabel('SOC battery [1 - 0]')
        plt.grid()
        plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # PLOT: POWER LOSSES -----------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    if pars1["pwr_params_mintime"]["pwr_behavior"] and plot_opts["pwr_losses"]:

        if pars1["pwr_params_mintime"]["simple_loss"]:
            plt.figure(12)
            plt.plot(s1[:-1], pwr1["machine"].p_loss_total)
            plt.plot(s1[:-1], pwr1["inverter"].p_loss_total)
            plt.plot(s1[:-1], pwr1["batt"].p_loss_total)
            plt.legend([r'$\it{P_\mathrm{loss,machine}}$', r'$\it{P_\mathrm{loss,inverter}}$',
                        r'$\it{P_\mathrm{loss,battery}}$'])
            plt.ylabel('Power loss ' + r'$\it{P_\mathrm{loss}}$' + ' in ' + r'kW')
        else:
            plt.figure(12)
            plt.subplot(311)
            plt.plot(s1[:-1], pwr1["machine"].p_loss_total)
            plt.plot(s1[:-1], pwr1["machine"].p_loss_copper)
            plt.plot(s1[:-1], pwr1["machine"].p_loss_stator_iron)
            plt.plot(s1[:-1], pwr1["machine"].p_loss_rotor)
            plt.ylabel('Power loss single machine\n' + r'$\it{P_\mathrm{loss}}$' + ' in ' + r'kW')
            plt.legend([r'$\it{P_\mathrm{loss,total}}$', r'$\it{P_\mathrm{loss,copper}}$',
                        r'$\it{P_\mathrm{loss,statorIron}}$', r'$\it{P_\mathrm{loss,rotor}}$'])
            plt.grid()
            plt.subplot(312)
            plt.plot(s1[:-1], pwr1["inverter"].p_loss_total)
            plt.plot(s1[:-1], pwr1["inverter"].p_loss_switch)
            plt.plot(s1[:-1], pwr1["inverter"].p_loss_cond)
            plt.legend([r'$\it{P_\mathrm{loss,total}}$', r'$\it{P_\mathrm{loss,switching}}$',
                        r'$\it{P_\mathrm{loss,conducting}}$'])
            plt.ylabel('Power loss single inverter\n' + r'$\it{P_\mathrm{loss}}$' + ' in ' + r'kW')
            plt.grid()
            plt.subplot(313)
            plt.plot(s1[:-1], pwr1["batt"].p_loss_total)
            plt.ylabel('Power loss battery\n' + r'$\it{P_\mathrm{loss}}$' + ' in ' + r'kW')

        plt.xlabel('distance ' + r'$\it{s}$' + ' in ' + r'$\it{m}$')
        plt.grid()
        plt.show()

# testing --------------------------------------------------------------------------------------------------------------
    if __name__ == "__main__":
        pass
