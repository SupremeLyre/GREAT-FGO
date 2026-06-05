import os
from pathlib import Path
from file_reader import read_position_data_enu, read_velocity_data, read_attitude_data
from statistics import calculate_all_metrics, save_comprehensive_statistics, print_statistics
from plot import plot_position_errors, plot_velocity_errors, plot_attitude_errors


def main():
    # Parameter settings
    # file_data = r'../sample_data/FGO_20250928/result/SEPT-RTK.fgo'
    # file_ref = r'../sample_data/FGO_20250928/ref/groundtruth_0928_GNSS.txt'
    # output_dir = r'./results/20250928'
    file_data = r'../sample_data/MSF_20201027/result/campus-01/CAMP-MSF.ins'
    file_ref = r'../sample_data/MSF_20201027/groundtruth/campus01-groundtruth.txt'
    output_dir = r'./results/20201027'

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    file_basename = os.path.splitext(os.path.basename(file_data))[0]

    # Initialize variables
    pos_errors_x, pos_errors_y, pos_errors_z = None, None, None
    vel_errors_x, vel_errors_y, vel_errors_z = None, None, None
    att_errors_x, att_errors_y, att_errors_z = None, None, None

    print("Processing position data (ENU coordinates)...")
    try:
        # Process position data with ENU conversion
        matched_time_pos, de, dn, du = read_position_data_enu(file_data, file_ref)
        pos_errors_x, pos_errors_y, pos_errors_z = de, dn, du
        pos_metrics = calculate_all_metrics(de, dn, du, "Position")
        plot_position_errors(matched_time_pos, de, dn, du, pos_metrics, file_basename, output_dir)
        print_statistics(pos_metrics, "Position")

        # # save enu data
        enu_data_path = os.path.join(output_dir, f'{file_basename}_enu_data.txt')
        with open(enu_data_path, 'w') as f:
            f.write("Time(s), East(m), North(m), Up(m)\n")
            for t, e, n, u in zip(matched_time_pos, de, dn, du):
                f.write(f"{t}, {e:.6f}, {n:.6f}, {u:.6f}\n")
        print(f"ENU coordinate data saved to: {enu_data_path}")

    except Exception as e:
        print(f"Error processing position data: {e}")

    # # Alternative：process vel data（for RTK/INS）
    print("\nProcessing velocity data...")
    try:
        # Process velocity data (if available)
        matched_time_vel, dvx, dvy, dvz = read_velocity_data(file_data, file_ref)
        vel_errors_x, vel_errors_y, vel_errors_z = dvx, dvy, dvz
        vel_metrics = calculate_all_metrics(dvx, dvy, dvz, "Velocity")
        plot_velocity_errors(matched_time_vel, dvx, dvy, dvz, vel_metrics, file_basename, output_dir)
        print_statistics(vel_metrics, "Velocity")
    except Exception as e:
        print(f"Velocity data not available or error: {e}")
    #
    # # Alternative：process att data（for RTK/INS）
    print("\nProcessing attitude data...")
    try:
        # Process attitude data (if available)
        matched_time_att, dpitch, droll, dyaw = read_attitude_data(file_data, file_ref)
        att_errors_x, att_errors_y, att_errors_z = dpitch, droll, dyaw
        att_metrics = calculate_all_metrics(dpitch, droll, dyaw, "Attitude")
        plot_attitude_errors(matched_time_att, dpitch, droll, dyaw, att_metrics, file_basename, output_dir)
        print_statistics(att_metrics, "Attitude")
    except Exception as e:
        print(f"Attitude data not available or error: {e}")

    # Save comprehensive statistics
    print("\nSaving comprehensive statistics...")
    try:
        save_comprehensive_statistics(
            pos_errors_x=pos_errors_x,
            pos_errors_y=pos_errors_y,
            pos_errors_z=pos_errors_z,
            vel_errors_x=vel_errors_x,
            vel_errors_y=vel_errors_y,
            vel_errors_z=vel_errors_z,
            att_errors_x=att_errors_x,
            att_errors_y=att_errors_y,
            att_errors_z=att_errors_z,
            file_data=file_data,
            file_ref=file_ref,
            output_path=os.path.join(output_dir, f'{file_basename}_statistics.txt')
        )
        print(f"Comprehensive statistics saved successfully!")
    except Exception as e:
        print(f"Error saving comprehensive statistics: {e}")

    print(f"\nAll processing completed! Results saved to: {output_dir}")


if __name__ == "__main__":
    main()