# OutageStatistic - 导航中断误差统计与绘图

本仓库用于统计导航中断误差以及绘制误差曲线图和阴影表示的中断区间，主要程序为OutageStatistic.m

OutageStatistic.m: 用于统计GLINS输出结果中断期间误差

GenerateOutages.m: 用于生成中断文件，导出文件也可用于GINS

MisalignedAngleValid.m: 通过绘制误差曲线验证GINS导出到目标IMU的安装角和杆臂是否正确

smooth_angle.m: 用于航向角插值之前进行平滑，避免相邻历元在0°和360°附近跳变

warp_angle.m: 将角度归化至[-180°,180°]



