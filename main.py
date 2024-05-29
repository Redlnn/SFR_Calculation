import math
from dataclasses import dataclass
from pathlib import Path
import time
from typing import Literal

import cv2
import numpy as np

ROI_WIDTH = 100
ROI_HEIGHT = 50

orig_img = cv2.imread((Path('edge_samples') / '1.bmp').as_posix())
start_time = time.time()
orig_img_gray = cv2.cvtColor(orig_img, cv2.COLOR_BGR2GRAY)
img = np.copy(orig_img)
gray_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

# 图片尺寸
width = img.shape[1]
height = img.shape[0]
diagonal_length = math.sqrt(width**2 + height**2)

# 计算中心坐标
center_x = width // 2
center_y = height // 2

cos_height = (height**2 + diagonal_length**2 - width**2) / (
    2 * height * diagonal_length
)
cos_width = (width**2 + diagonal_length**2 - height**2) / (2 * width * diagonal_length)


@dataclass
class RoiImage:
    img: cv2.typing.MatLike
    roi: tuple[tuple[int, int], tuple[int, int]]
    direction: (
        Literal['C'] | Literal['UL'] | Literal['UR'] | Literal['BL'] | Literal['BR']
    )
    field: float


# cv2.rectangle(
#     img,
#     (center_x - 5, center_y - 5),
#     (center_x + 5, center_y + 5),
#     (0, 0, 255),
#     10,
# )
# cv2.rectangle(
#     img,
#     (center_x - 100, center_y - 100),
#     (center_x + 100, center_y + 100),
#     (0, 0, 255),
#     10,
# )
rois: list[RoiImage] = [
    RoiImage(
        img=orig_img[center_y - 100 : center_y + 100, center_x - 100 : center_x + 100],
        roi=(
            (center_x - 100, center_y - 100),
            (center_x + 100, center_y + 100),
        ),
        direction='C',
        field=0,
    )
]
for field in (0.3, 0.5, 0.75, 0.9):
    # 计算偏移向量
    offset = field * (diagonal_length / 2)
    offset_x = int(offset * cos_width)
    offset_y = int(offset * cos_height)

    # 图像中间左上角坐标
    # TODO: 超出图像的话强行压回来
    pos_ul = 'UL', (center_x - offset_x, center_y - offset_y)
    pos_ur = 'UR', (center_x + offset_x, center_y - offset_y)
    pos_bl = 'BL', (center_x - offset_x, center_y + offset_y)
    pos_br = 'BR', (center_x + offset_x, center_y + offset_y)

    for pos in (pos_ul, pos_ur, pos_bl, pos_br):
        # cv2.rectangle(
        #     img,
        #     (pos[1][0] - 5, pos[1][1] - 5),
        #     (pos[1][0] + 5, pos[1][1] + 5),
        #     (0, 0, 255),
        #     10,
        # )
        # cv2.rectangle(
        #     img,
        #     (pos[1][0] - 100, pos[1][1] - 100),
        #     (pos[1][0] + 100, pos[1][1] + 100),
        #     (0, 0, 255),
        #     10,
        # )
        rois.append(
            RoiImage(
                img=orig_img[
                    pos[1][1] - 100 : pos[1][1] + 100, pos[1][0] - 100 : pos[1][0] + 100
                ],
                roi=(
                    (pos[1][0] - 100, pos[1][1] - 100),
                    (pos[1][0] + 100, pos[1][1] + 100),
                ),
                direction=pos[0],
                field=field,
            )
        )

for i in rois:

    '''
    shi-Tomas算法是对Harris角点检测算法的改进,一般会比Harris算法得到更好的角点.
    Harris算法的角点响应函数是将矩阵M的行列式值与M的迹相减,利用差值判断是否为角点,
    后来shi和Tomas提出改进的方法是,若矩阵M的两个特征值中较小的一个大于阈值,则认为
    他是角点,即:
        R=min(a1,a2)

    API:
        corners=cv.goodFeaturesToTrack(image, maxcorners, qualityLevel, minDistance)
        参数:
            image: 输入的灰度图像
            maxCorners: 获取角点数的数目
            qualityLevel: 该参数指出最低可接受的角点质量水平,在0~1之间
            minDistance: 角点之间的最小欧氏距离,避免得到相邻特征点
        返回:
            corners:搜索到的角点,在这里所有低于质量水平的角点被排除,然后把合格的角点按照质量排序,
            然后将质量较好的角点附近(小于最小欧氏距离)的角点删除,最后找到maxCorners个角点返回

    '''
    corners = cv2.goodFeaturesToTrack(
        cv2.cvtColor(i.img, cv2.COLOR_BGR2GRAY), 1000, 0.5, 10
    )
    # TODO: 占用时间长

    x, y = corners[0].ravel()
    x, y = int(x), int(y)

    # 防止ROI选出ROI识别区域
    # 竖直方向超出
    if int(i.roi[0][1] + y - ROI_WIDTH - 20) > i.roi[0][1]:
        cv2.rectangle(
            img,
            (
                i.roi[0][0] + x - ROI_HEIGHT // 3,
                i.roi[0][1] + y - ROI_WIDTH - 20,
            ),
            (
                i.roi[0][0] + x + ROI_HEIGHT - ROI_HEIGHT // 3,
                i.roi[0][1] + y - 20,
            ),
            (0, 0, 255),
            2,
        )
    else:
        cv2.rectangle(
            img,
            (
                i.roi[0][0] + x - ROI_HEIGHT + ROI_HEIGHT // 3,
                i.roi[0][1] + y + 20,
            ),
            (
                i.roi[0][0] + x + ROI_HEIGHT // 3,
                i.roi[0][1] + y + ROI_WIDTH + 20,
            ),
            (0, 255, 0),
            2,
        )

    # 水平方向超出
    if (i.roi[0][0] + x - ROI_WIDTH - 20) > i.roi[0][0]:
        cv2.rectangle(
            img,
            (
                i.roi[0][0] + x - ROI_WIDTH - 20,
                i.roi[0][1] + y - ROI_HEIGHT + ROI_HEIGHT // 3,
            ),
            (
                i.roi[0][0] + x - 20,
                i.roi[0][1] + y + ROI_HEIGHT // 3,
            ),
            (0, 0, 255),
            2,
        )
    else:
        cv2.rectangle(
            img,
            (
                i.roi[0][0] + x + 20,
                i.roi[0][1] + y - ROI_HEIGHT // 3,
            ),
            (
                i.roi[0][0] + x + ROI_WIDTH + 20,
                i.roi[0][1] + y + ROI_HEIGHT - ROI_HEIGHT // 3,
            ),
            (0, 255, 0),
            2,
        )

print(time.time() - start_time)

cv2.namedWindow("result", 0)
cv2.resizeWindow("result", 1600, 1115)  # 设置窗口大小
cv2.imshow('result', img)
# cv2.imwrite('test.png', rois[0].img)
cv2.waitKey()
