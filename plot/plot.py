import matplotlib.pyplot as plt
import os


def plot_errors(matched_time, errors_x, errors_y, errors_z, metrics,
                file_basename, output_dir, data_type="Position",
                colors=None, labels=None, ylabel_suffix="",
                ylim=(-1, 1), threshold=0.1):
    """
    General plotting function for ENU coordinates
    """
    if colors is None:
        colors = ['red', 'green', 'blue']
    if labels is None:
        # 默认使用ENU标签
        labels = ['East', 'North', 'Up']

    fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    errors = [errors_x, errors_y, errors_z]

    for i, (ax, error, color, label) in enumerate(zip(axs, errors, colors, labels)):
        m = metrics[['x', 'y', 'z'][i]]
        ax.scatter(matched_time, error, color=color, s=10, label=label, alpha=0.7)
        ax.axhline(threshold, color='gray', linestyle='--', linewidth=1, alpha=0.5,
                   label=f'±{threshold}{ylabel_suffix} Threshold')
        ax.axhline(-threshold, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        ax.axhline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.8)
        ax.set_ylabel(f'{label} Error {ylabel_suffix}')
        ax.set_ylim(ylim)
        ax.grid(True, alpha=0.3)
        ax.legend()

        # 设置标题
        title = f'{label} {data_type} Errors (MAE: {m["mae"]:.3f}{ylabel_suffix}, RMSE: {m["rmse"]:.3f}{ylabel_suffix})'
        ax.set_title(title)

    axs[2].set_xlabel('Time (s)')
    plt.suptitle(f'{data_type} Errors: {file_basename}', fontsize=14)
    plt.tight_layout()

    # 保存图片
    output_path = os.path.join(output_dir, f'{file_basename}_{data_type.lower()}.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()

    return output_path


def plot_position_errors(matched_time, de, dn, du, metrics, file_basename, output_dir):
    """绘制ENU位置误差图"""
    return plot_errors(matched_time, de, dn, du, metrics, file_basename, output_dir,
                       data_type="Position", colors=['red', 'green', 'blue'],
                       labels=['East', 'North', 'Up'], ylabel_suffix="(m)",
                       ylim=(-10, 10), threshold=0.5)


def plot_velocity_errors(matched_time, dvx, dvy, dvz, metrics, file_basename, output_dir):
    """绘制速度误差图"""
    return plot_errors(matched_time, dvx, dvy, dvz, metrics, file_basename, output_dir,
                       data_type="Velocity", colors=['red', 'green', 'blue'],
                       labels=['Vx', 'Vy', 'Vz'], ylabel_suffix="(m/s)",
                       ylim=(-0.2, 0.2), threshold=0.05)


def plot_attitude_errors(matched_time, dpitch, droll, dyaw, metrics, file_basename, output_dir):
    """绘制姿态误差图"""
    return plot_errors(matched_time, dpitch, droll, dyaw, metrics, file_basename, output_dir,
                       data_type="Attitude", colors=['red', 'green', 'blue'],
                       labels=[ 'Yaw', 'Pitch', 'Roll'], ylabel_suffix="(°)",
                       ylim=(-3, 3), threshold=1)