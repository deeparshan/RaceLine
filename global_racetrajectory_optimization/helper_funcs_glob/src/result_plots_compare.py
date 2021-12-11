import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import trajectory_planning_helpers
import os

my_path = os.getcwd()


def result_plots_compare(plot_opts1: dict,
                 width_veh_opt1: float,
                 width_veh_real1: float,
                 refline1: np.ndarray,
                 bound1_imp1: np.ndarray,
                 bound2_imp1: np.ndarray,
                 bound1_interp1: np.ndarray,
                 bound2_interp1: np.ndarray,
                 trajectory1: np.ndarray,
                 plot_opts2: dict,
                 width_veh_opt2: float,
                 width_veh_real2: float,
                 refline2: np.ndarray,
                 bound1_imp2: np.ndarray,
                 bound2_imp2: np.ndarray,
                 bound1_interp2: np.ndarray,
                 bound2_interp2: np.ndarray,
                 trajectory2: np.ndarray) -> None:
    """
    Created by:
    Alexander Heilmeier

    Documentation:
    This function plots several figures containing relevant trajectory information after trajectory optimization.

    Inputs:
    plot_opts:      dict containing the information which figures to plot
    width_veh_opt:  vehicle width used during optimization in m
    width_veh_real: real vehicle width in m
    refline:        contains the reference line coordinates [x_m, y_m]
    bound1_imp:     first track boundary (as imported) (mostly right) [x_m, y_m]
    bound2_imp:     second track boundary (as imported) (mostly left) [x_m, y_m]
    bound1_interp:  first track boundary (interpolated) (mostly right) [x_m, y_m]
    bound2_interp:  second track boundary (interpolated) (mostly left) [x_m, y_m]
    trajectory:     trajectory data [s_m, x_m, y_m, psi_rad, kappa_radpm, vx_mps, ax_mps2]
    """

    if plot_opts1["raceline"]:
        # calculate vehicle boundary points (including safety margin in vehicle width)
        normvec_normalized_opt1 = trajectory_planning_helpers.calc_normal_vectors.\
            calc_normal_vectors(trajectory1[:, 3])
        normvec_normalized_opt2 = trajectory_planning_helpers.calc_normal_vectors.\
            calc_normal_vectors(trajectory2[:, 3])

        veh_bound1_virt1 = trajectory1[:, 1:3] + normvec_normalized_opt1 * width_veh_opt1 / 2
        veh_bound2_virt1 = trajectory1[:, 1:3] - normvec_normalized_opt1 * width_veh_opt1 / 2
        
        veh_bound1_virt2 = trajectory2[:, 1:3] + normvec_normalized_opt2 * width_veh_opt2 / 2
        veh_bound2_virt2 = trajectory2[:, 1:3] - normvec_normalized_opt2 * width_veh_opt2 / 2

        veh_bound1_real1 = trajectory1[:, 1:3] + normvec_normalized_opt1 * width_veh_real1 / 2
        veh_bound2_real1 = trajectory1[:, 1:3] - normvec_normalized_opt1 * width_veh_real1 / 2
        
        veh_bound1_real2 = trajectory2[:, 1:3] + normvec_normalized_opt2 * width_veh_real2 / 2
        veh_bound2_real2 = trajectory2[:, 1:3] - normvec_normalized_opt2 * width_veh_real2 / 2

        point1_arrow1 = refline1[0]
        point2_arrow1 = refline1[3]
        vec_arrow1 = point2_arrow1 - point1_arrow1

        #point1_arrow2 = refline2[0]
        #point2_arrow2 = refline2[3]
        #vec_arrow2 = point2_arrow2 - point1_arrow2
        # plot track including optimized path
        plt.figure()
        plt.plot(refline1[:, 0], refline1[:, 1], "k--", linewidth=0.7)
        plt.plot(veh_bound1_virt1[1:, 0], veh_bound1_virt1[:, 1], "b", linewidth=0.5)
        plt.plot(veh_bound2_virt1[:, 0], veh_bound2_virt1[:, 1], "b", linewidth=0.5)
        plt.plot(veh_bound1_real1[:, 0], veh_bound1_real1[:, 1], "c", linewidth=0.5)
        plt.plot(veh_bound2_real1[:, 0], veh_bound2_real1[:, 1], "c", linewidth=0.5)
        plt.plot(bound1_interp1[:, 0], bound1_interp1[:, 1], "k-", linewidth=0.7)
        plt.plot(bound2_interp1[:, 0], bound2_interp1[:, 1], "k-", linewidth=0.7)
        plt.plot(trajectory1[:, 1], trajectory1[:, 2], "r-", linewidth=0.7)
        
        plt.plot(veh_bound1_virt2[:, 0], veh_bound1_virt2[:, 1], "b", linewidth=0.5)
        plt.plot(veh_bound2_virt2[:, 0], veh_bound2_virt2[:, 1], "b", linewidth=0.5)
        plt.plot(veh_bound1_real2[:, 0], veh_bound1_real2[:, 1], "c", linewidth=0.5)
        plt.plot(veh_bound2_real2[:, 0], veh_bound2_real2[:, 1], "c", linewidth=0.5)
        plt.plot(bound1_interp2[:, 0], bound1_interp2[:, 1], "k-", linewidth=0.7)
        plt.plot(bound2_interp2[:, 0], bound2_interp2[:, 1], "k-", linewidth=0.7)
        plt.plot(trajectory2[:, 1], trajectory2[:, 2], "r-", linewidth=0.7)

        if plot_opts1["imported_bounds"] and bound1_imp1 is not None and bound2_imp1 is not None:
            plt.plot(bound1_imp1[:, 0], bound1_imp1[:, 1], "y-", linewidth=0.7)
            plt.plot(bound2_imp1[:, 0], bound2_imp1[:, 1], "y-", linewidth=0.7)

        plt.grid()
        ax = plt.gca()
        ax.arrow(point1_arrow1[0], point1_arrow1[1], vec_arrow1[0], vec_arrow1[1],
                 head_width=7.0, head_length=7.0, fc='g', ec='g')
        ax.set_aspect("equal", "datalim")
        plt.xlabel("east in m")
        plt.ylabel("north in m")
        plt.savefig(my_path + '/outputs/Optimized_Path.png')
        plt.show()


# testing --------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    pass
