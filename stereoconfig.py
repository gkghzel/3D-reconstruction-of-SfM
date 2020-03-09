import numpy as np
import cv2

# 双目相机参数
class stereoCameral(object):
    def __init__(self):

        # 左相机内参数
        self.cam_matrix_left = np.array([[530.9002, 0, 137.63037], [0, 581.00362, 162.32884], [0, 0, 1]])
        # 右相机内参数
        self.cam_matrix_right = np.array([[524.84413, 0, 217.17358], [0, 577.11024, 150.76379], [0, 0, 1]])

        # 左右相机畸变系数:[k1, k2, p1, p2, k3]
        self.distortion_l = np.array([-0.2575, 0.6231, 0.0366, -0.0108, 0])
        self.distortion_r = np.array([-0.25745, 0.62307, 0.0366, -0.01082, 0])

        # 旋转矩阵
        om = np.array([-0.009954307627539,-0.042566273591172, 0.011454074192639])
        self.R = cv2.Rodrigues(om)[0]  # 使用Rodrigues变换将om变换为R
        # 平移矩阵
        self.T = np.array([-5.49238, 0.04267, -0.39886])